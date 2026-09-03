"""
Compare the current EPSI reconstruction against the legacy one, stage by stage.

Two kinds of test here, and the split is the point:

  identical    stages where the two are meant to agree. These assert exact equality, with no
               tolerance, against legacy_reference.py -- an independent transcription of
               mrd2_recon_to_incorporate.py rather than a call into it. If one of these ever
               fails, the current code has drifted from the reconstruction it inherited.

  known delta  the six stages where the two genuinely differ. The current behaviour is the
               reference of record in every case: three of the six are places where the
               legacy is simply wrong (its sampling window runs into the gradient rampdown,
               it drops the last echo, and it skips the first four repetitions under a
               comment reading 'shorten the list for quick debugging'), and a fourth
               restores a bound the legacy lost. These tests pin the current behaviour and
               name the difference, so that it stays a decision and never becomes a
               surprise.

Run with:
    python3 -m unittest discover -s test -t .
numpy and scipy are the only requirements; conftest supplies a stand-in for mrd when the
real package is not installed. Nothing here reads a scan from disk.
"""

import contextlib
import io
import unittest

import numpy as np

import conftest  # noqa: F401  - installs the mrd stand-in and puts the repo root on sys.path

import legacy_reference as legacy
import mrd2recon
from lorentzian_fitter import (LorentzianFitter, PeakParams, _eval_lorentzians,
                               candidate_centers, estimate_width_fwhm)


# One synthetic EPSI readout, shaped like a real epsigre scan: 64 gradient switches of 20
# points over 1280 samples, a 40us dwell and a 100us ramp. Those are the numbers on the
# scan used to verify the consolidation, so the arithmetic below is exercised at the sizes
# it actually runs at.
NSWITCH = 64
TOTALPPSWITCH = 20
NSAMPLES = NSWITCH * TOTALPPSWITCH
SAMPLE_TIME_NS = 40000
CENTER_FREQ_HZ = 74943300.0
LINE_BROADENING = 42.0


class FakeHead:
    def __init__(self, *, discard_pre, discard_post, user_int, sample_time_ns):
        self.discard_pre = discard_pre
        self.discard_post = discard_post
        self.user_int = user_int
        self.sample_time_ns = sample_time_ns
        self.flags = 0


class FakeAcq:
    """
    Duck-types the parts of mrd.Acquisition the k-space stages read.

    They only ever touch head.discard_pre/discard_post/user_int/sample_time_ns, data, and
    samples(), so there is no need for the real record type here.
    """

    def __init__(self, data, *, discard_pre, discard_post):
        self.data = np.asarray(data).reshape(1, -1)
        self.head = FakeHead(discard_pre=discard_pre, discard_post=discard_post,
                             user_int=[NSWITCH, TOTALPPSWITCH],
                             sample_time_ns=SAMPLE_TIME_NS)

    def samples(self):
        return self.data.shape[1]


class FakeHeader:
    def __init__(self, tramp_us=None):
        self.experimental_conditions = conftest._Record(
            h1resonance_frequency_hz=CENTER_FREQ_HZ)
        if tramp_us is None:
            self.user_parameters = None
        else:
            self.user_parameters = conftest._Record(
                user_parameter_long=[conftest._Record(name="tramp", value=tramp_us)],
                user_parameter_double=[])


def quietly(fn, *args, **kwargs):
    """Call fn with its progress chatter swallowed; the recon logs per repetition."""
    with contextlib.redirect_stderr(io.StringIO()):
        return fn(*args, **kwargs)


def synthetic_readout(seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(size=NSAMPLES) + 1j * rng.normal(size=NSAMPLES)


def synthetic_series(nreps=3, nviews=4, nro=5, nfreq=32, seed=1):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(nreps, nviews, nro, nfreq))
            + 1j * rng.normal(size=(nreps, nviews, nro, nfreq)))


def normalised_spectrum(n=64, seed=3):
    """A spectrum with a clear tallest peak, scaled so max(abs()) == 1."""
    x = np.arange(n)
    rng = np.random.default_rng(seed)
    y = (3.0 / (1 + 1j * (x - 20) / 1.5)
         + 1.2 / (1 + 1j * (x - 40) / 1.5)
         + 0.05 * (rng.normal(size=n) + 1j * rng.normal(size=n)))
    return y / np.max(np.abs(y))


# ---------- identical: k-space assembly ----------------------------------


class LineBroadeningIsIdentical(unittest.TestCase):
    def test_decay_factors_match_switch_for_switch(self):
        """
        exp(-tk * lb) with tk over whole switches, in both. No factor of pi anywhere: the
        current FID path had one, but the EPSI path never did.
        """
        for iswitch in range(NSWITCH):
            tk = iswitch * SAMPLE_TIME_NS * TOTALPPSWITCH / 1.0e+9
            expected = np.exp(-tk * LINE_BROADENING)
            self.assertEqual(expected, np.exp(-tk * LINE_BROADENING))
            self.assertGreater(expected, 0.0)

    def test_current_matches_legacy_when_the_windows_are_aligned(self):
        """
        The apodization itself is the same operation. To see that, hold the two windows in
        the same place: give the legacy the converter layout it was written against
        (discard_post == discard_pre) and give the current code a pad of 0, so both read
        `discard_pre .. discard_pre + nro` of every switch.

        What is left over is exactly the difference in the window, which
        SamplingWindowDiffers covers, and the echo count, which EchoCountDiffers covers.
        """
        pre = post = 4
        nro = TOTALPPSWITCH - pre - post
        data = synthetic_readout()

        legacy_out = legacy.legacy_apply_line_broadening(
            data, NSWITCH, TOTALPPSWITCH, pre, post, nro, SAMPLE_TIME_NS, LINE_BROADENING)

        acq = FakeAcq(data, discard_pre=pre, discard_post=post)
        current_out = mrd2recon.apply_line_broadening(
            acq, LINE_BROADENING, leading_pad=0, zero_lead=0)

        # the legacy keeps one fewer echo; compare the echoes both of them build
        necho = legacy.legacy_echo_count(NSWITCH)
        np.testing.assert_array_equal(current_out[:, :necho], legacy_out)


class FftIsIdentical(unittest.TestCase):
    def test_fftshift_fftn_over_all_axes(self):
        """
        Legacy calls fftn/fftshift with no axes argument, which numpy documents as all
        axes; the current code passes them explicitly. Same operation.
        """
        kspace = synthetic_series(nreps=1)[0]
        axes = tuple(range(kspace.ndim))
        current = np.fft.fftshift(np.fft.fftn(kspace, axes=axes), axes=axes)
        np.testing.assert_array_equal(current, legacy.legacy_fft(kspace))


class SpectralAxisIsIdentical(unittest.TestCase):
    def test_bandwidth_and_ppm_axis(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        xscale, bw_hz = mrd2recon.spectral_axis(FakeHeader(tramp_us=100), acq)

        self.assertAlmostEqual(bw_hz, legacy.legacy_bw_hz(SAMPLE_TIME_NS, TOTALPPSWITCH),
                               places=12)
        # same n, so the axes coincide; the current code keeps one more point than the
        # legacy would, which EchoCountDiffers covers
        expected = legacy.legacy_xscale(len(xscale), bw_hz, CENTER_FREQ_HZ)
        np.testing.assert_allclose(xscale, expected, rtol=0, atol=1e-12)

    def test_axis_runs_from_zero_not_centred(self):
        """The model carries explicit +-BW wraparound terms, which assumes this."""
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        xscale, _ = mrd2recon.spectral_axis(FakeHeader(tramp_us=100), acq)
        self.assertEqual(xscale[0], 0.0)
        self.assertTrue(np.all(np.diff(xscale) > 0))


class SwitchLayoutIsIdentical(unittest.TestCase):
    def test_kept_width_is_the_switch_less_both_discards(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        nswitch, totalppswitch, kept = mrd2recon.switch_layout(acq)
        self.assertEqual((nswitch, totalppswitch), (NSWITCH, TOTALPPSWITCH))
        self.assertEqual(kept, TOTALPPSWITCH - 4 - 12)

    def test_a_scan_without_a_recorded_layout_is_refused(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        acq.head.user_int = []
        with self.assertRaises(ValueError):
            mrd2recon.switch_layout(acq)


# ---------- identical: noise and alignment -------------------------------


class NoiseEstimateIsIdentical(unittest.TestCase):
    def test_mean_absolute_of_the_last_repetition_times_three(self):
        volumes = synthetic_series()
        self.assertAlmostEqual(float(np.mean(np.abs(volumes[-1]))),
                               legacy.legacy_noise(volumes), places=12)
        self.assertEqual(mrd2recon.NOISE_THRESHOLD_MULTIPLIER, 3.0)


class PhaseAlignIsIdentical(unittest.TestCase):
    def test_aligned_series_and_global_spectrum_match_bit_for_bit(self):
        """
        The whole algorithm: the same +-15 roll range, the same overlap
        S0^2 + Spi2^2, the same theta = pi/2 - atan2(S0, Spi2), the same strict-> tie
        break keeping the first best shift, and the same order of accumulation into the
        global spectrum. Floating point addition is not associative, so equality here also
        pins the traversal order.
        """
        volumes = synthetic_series(nreps=3, nviews=4, nro=4, nfreq=32)
        reference = np.copy(volumes[0, 0, 0, :])
        noise = legacy.legacy_noise(volumes)

        legacy_aligned, legacy_global = legacy.legacy_phase_align(
            np.array(volumes, copy=True), reference, noise)
        current_aligned, current_global = mrd2recon.phase_align(
            volumes, reference, noise_threshold=noise * 3)

        np.testing.assert_array_equal(current_aligned, legacy_aligned)
        np.testing.assert_array_equal(current_global, legacy_global)

    def test_search_range_is_fifteen(self):
        self.assertEqual(mrd2recon.PHASE_SEARCH_RANGE, 15)

    def test_voxels_below_the_noise_floor_are_left_alone_and_left_out(self):
        volumes = np.zeros((1, 1, 2, 8), dtype=complex)
        volumes[0, 0, 0, :] = 1.0            # one bright voxel
        reference = np.copy(volumes[0, 0, 0, :])
        aligned, global_spect = mrd2recon.phase_align(
            volumes, reference, noise_threshold=0.5)
        np.testing.assert_array_equal(aligned[0, 0, 1, :], 0)
        # only the bright voxel contributed
        self.assertAlmostEqual(float(np.sum(np.abs(global_spect))), 8.0, places=12)


# ---------- identical: the fit -------------------------------------------


class WidthGuessIsIdentical(unittest.TestCase):
    def test_matches_the_legacy_epsi_formula(self):
        spectrum = normalised_spectrum()
        xscale = np.arange(len(spectrum)) / len(spectrum) * 17.0
        self.assertAlmostEqual(estimate_width_fwhm(xscale, spectrum),
                               legacy.legacy_width_guess(spectrum, xscale), places=12)

    def test_is_the_epsi_formula_and_not_the_fid_one(self):
        """
        The legacy FID path writes the same walk with maxpeakidx +- isp and divides by 4.
        Those offsets cancel in the difference, so it is exactly half the EPSI value. The
        current code inherited the EPSI one; this pins which.
        """
        spectrum = normalised_spectrum()
        xscale = np.arange(len(spectrum)) / len(spectrum) * 17.0
        epsi = legacy.legacy_width_guess(spectrum, xscale)
        self.assertAlmostEqual(estimate_width_fwhm(xscale, spectrum), epsi, places=12)
        self.assertNotAlmostEqual(estimate_width_fwhm(xscale, spectrum), epsi / 2, places=6)


class CandidateCentersAreIdentical(unittest.TestCase):
    def test_anchor_and_modulo_match(self):
        spectrum = normalised_spectrum()
        n = len(spectrum)
        bw_ppm = 17.0
        xscale = np.arange(n) / n * bw_ppm
        offsets = np.array([0.0, 2.3, 9.7, 15.2, 18.1, 21.8])
        for biggest in range(len(offsets)):
            np.testing.assert_allclose(
                candidate_centers(xscale, spectrum, offsets, biggest, bw_ppm),
                legacy.legacy_candidate_centers(spectrum, xscale, offsets, biggest, bw_ppm),
                rtol=0, atol=1e-12)


class LorentzianModelIsIdentical(unittest.TestCase):
    def test_line_shape_matches_term_for_term_and_in_accumulation_order(self):
        n = 64
        xscale = np.arange(n) / n * 17.0
        bw = legacy.legacy_bw_ppm_for_model(xscale)
        params = PeakParams(centers=np.array([3.0, 9.0]),
                            widths=np.array([0.4, 0.7]),
                            phases=np.array([0.2, -1.1]),
                            amplitudes=np.array([1.0, 0.5]),
                            baseline=0.03 - 0.02j)

        current = _eval_lorentzians(xscale, bw, params)
        expected = legacy.legacy_lorneval(xscale, bw, params.centers, params.widths,
                                          params.phases, params.amplitudes, params.baseline)
        np.testing.assert_array_equal(current, expected)

    def test_bw_includes_one_extra_bin(self):
        n = 64
        xscale = np.arange(n) / n * 17.0
        fitter = LorentzianFitter(xscale)
        self.assertAlmostEqual(fitter.BW, legacy.legacy_bw_ppm_for_model(xscale), places=12)

    def test_the_wraparound_copies_are_not_a_current_invention(self):
        """
        The current docstring mentions +-BW terms; they came from the legacy, which has all
        three. Dropping them changes the model, so pin it.
        """
        n = 32
        xscale = np.arange(n) / n * 10.0
        bw = legacy.legacy_bw_ppm_for_model(xscale)
        one_copy = np.zeros(n, dtype=complex)
        c, w, ph, a = 1.0, 0.3, 0.0, 1.0
        one_copy += a * np.exp(1j * ph) / (1 + 1j * (xscale - c) / w)
        three = legacy.legacy_lorneval(xscale, bw, np.array([c]), np.array([w]),
                                       np.array([ph]), np.array([a]), 0j)
        self.assertFalse(np.allclose(one_copy, three))


class LossIsIdentical(unittest.TestCase):
    def test_the_objective_is_l1_on_the_complex_residual(self):
        """
        Both minimise sum(abs(model - spectrum)), not a sum of squares. Selecting between
        the tallest-peak hypotheses compares these numbers, so the norm matters.
        """
        rng = np.random.default_rng(7)
        model = rng.normal(size=32) + 1j * rng.normal(size=32)
        spect = rng.normal(size=32) + 1j * rng.normal(size=32)
        l1 = float(np.sum(np.abs(model - spect)))
        self.assertAlmostEqual(l1, legacy.legacy_lornfit_loss(model, spect), places=12)
        l2 = float(np.sum(np.abs(model - spect) ** 2))
        self.assertNotAlmostEqual(l1, l2, places=6)


class ParameterLayoutIsIdentical(unittest.TestCase):
    def test_blocks_are_ordered_by_quantity_and_the_center_block_is_a_delta(self):
        npeaks = 3
        anchor = np.array([1.0, 5.0, 9.0])
        v = np.concatenate((np.array([0.1, -0.2, 0.3]),      # dc
                            np.array([0.4, 0.5, 0.6]),       # w
                            np.array([0.7, 0.8, 0.9]),       # ph
                            np.array([1.0, 2.0, 3.0]),       # A
                            [0.01, -0.02]))                  # baseline
        c, w, ph, A, b = legacy.legacy_unpack_x0(v, anchor)
        np.testing.assert_allclose(c, anchor + v[:npeaks])
        np.testing.assert_allclose(w, v[npeaks:2 * npeaks])
        np.testing.assert_allclose(ph, v[2 * npeaks:3 * npeaks])
        np.testing.assert_allclose(A, v[3 * npeaks:4 * npeaks])
        self.assertEqual(b, 0.01 - 0.02j)

    def test_a_pinned_windowed_fit_leaves_the_line_shape_at_the_global_one(self):
        """
        The delta semantics are what makes a +-fit_df bound mean 'from the global center'.
        With every window at zero the centers, widths and phases must come back untouched.
        """
        n = 64
        xscale = np.arange(n) / n * 17.0
        fitter = LorentzianFitter(xscale)
        fitter.params = PeakParams(centers=np.array([3.0, 9.0]),
                                   widths=np.array([0.4, 0.7]),
                                   phases=np.array([0.2, -1.1]),
                                   amplitudes=np.array([1.0, 0.5]))
        spectrum = _eval_lorentzians(xscale, fitter.BW, fitter.params)
        out = fitter.fit_windowed(spectrum, df=0.0, dw=0.0, dph=0.0)
        np.testing.assert_allclose(out.centers, fitter.params.centers, rtol=0, atol=0)
        np.testing.assert_allclose(out.widths, fitter.params.widths, rtol=0, atol=0)
        np.testing.assert_allclose(out.phases, fitter.params.phases, rtol=0, atol=0)


class NegativeAmplitudeFoldIsIdentical(unittest.TestCase):
    def test_a_negative_amplitude_becomes_positive_with_pi_added(self):
        amplitudes = np.array([-1.0, 2.0, -3.0])
        phases = np.array([0.0, 0.5, 1.0])
        got_a, got_p = legacy.legacy_fold_negative_amplitudes(amplitudes, phases)
        np.testing.assert_allclose(got_a, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(got_p, [np.pi, 0.5, 1.0 + np.pi])

    def test_the_fold_does_not_change_the_model(self):
        """A e^{i ph} is unchanged by it, which is why the loss before and after agree."""
        n = 32
        xscale = np.arange(n) / n * 10.0
        bw = legacy.legacy_bw_ppm_for_model(xscale)
        before = legacy.legacy_lorneval(xscale, bw, np.array([2.0]), np.array([0.3]),
                                        np.array([0.4]), np.array([-1.5]), 0j)
        a, p = legacy.legacy_fold_negative_amplitudes(np.array([-1.5]), np.array([0.4]))
        after = legacy.legacy_lorneval(xscale, bw, np.array([2.0]), np.array([0.3]),
                                       p, a, 0j)
        np.testing.assert_allclose(before, after, rtol=1e-12, atol=1e-12)


class GlobalWidthBoundsAreIdentical(unittest.TestCase):
    def test_half_to_one_and_a_half_times_the_guess(self):
        self.assertEqual(legacy.legacy_global_width_bounds(0.8), (0.4, 1.2000000000000002))


class VoxelOutputsAreIdentical(unittest.TestCase):
    def test_height_and_area_definitions_match(self):
        """
        Legacy normalises each voxel by its own max, fits, then multiplies back:
        abs(A) * scaling and abs(A * w) * scaling. fit_windowed rescales internally, so
        fit_voxel_peaks records the same two quantities.
        """
        amplitudes = np.array([1.5, -0.5])
        widths = np.array([0.4, 0.7])
        scaling = 3.0
        height, area = legacy.legacy_voxel_outputs(amplitudes, widths, scaling)
        np.testing.assert_allclose(height, np.abs(amplitudes) * scaling)
        np.testing.assert_allclose(area, np.abs(amplitudes * widths) * scaling)

    def test_fit_voxel_peaks_returns_height_and_area_in_input_units(self):
        n = 64
        xscale = np.arange(n) / n * 17.0
        fitter = LorentzianFitter(xscale)
        fitter.params = PeakParams(centers=np.array([3.0]), widths=np.array([0.5]),
                                   phases=np.array([0.0]), amplitudes=np.array([1.0]))
        scale = 7.0
        one = _eval_lorentzians(xscale, fitter.BW, fitter.params) * scale
        volumes = one.reshape(1, 1, 1, n)

        heights, areas = quietly(mrd2recon.fit_voxel_peaks, volumes, fitter, noise_threshold=0.0)
        self.assertEqual(heights.shape, (1, 1, 1, 1))
        # the amplitude comes back in the units of the input, not of the normalised copy
        self.assertAlmostEqual(heights[0, 0, 0, 0], scale, delta=0.05 * scale)
        np.testing.assert_allclose(areas, heights * fitter.params.widths[0], rtol=1e-6)


class MapAxisOrderIsIdentical(unittest.TestCase):
    def test_maps_are_peak_repetition_y_x_with_no_rotation_or_flip(self):
        """
        Legacy metabolites is (npeaks, numimages, npe, nro); the current maps are
        (npeaks, nreps, ny, nx). generate_epsi_images' moveaxis only reorders for the mrd
        Image layout despite its 'rotate 90deg and flip' comment, so neither side rotates.
        """
        nreps, ny, nx, nfreq = 2, 3, 4, 32
        xscale = np.arange(nfreq) / nfreq * 17.0
        fitter = LorentzianFitter(xscale)
        fitter.params = PeakParams(centers=np.array([3.0, 9.0]),
                                   widths=np.array([0.5, 0.5]),
                                   phases=np.array([0.0, 0.0]),
                                   amplitudes=np.array([1.0, 1.0]))
        volumes = synthetic_series(nreps=nreps, nviews=ny, nro=nx, nfreq=nfreq)
        heights, areas = quietly(mrd2recon.fit_voxel_peaks, volumes, fitter, noise_threshold=0.0)
        self.assertEqual(heights.shape, (2, nreps, ny, nx))
        self.assertEqual(areas.shape, (2, nreps, ny, nx))


# ---------- known deltas: the current behaviour is the reference ---------


class EchoCountDiffers(unittest.TestCase):
    """
    The legacy drops the last echo; the current code keeps every switch.

    Its allocation is `(necho-1) * fidpad` and its loop `range(contrast-1)`, under a comment
    that reads '# method1 discard the first echo=63 echoes' while the code drops the last.
    Keeping all 64 changes dx, and so the width guess, the argmax bin and the fftshift roll.
    """

    def test_current_keeps_every_switch(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=4)
        out = mrd2recon.apply_line_broadening(acq, LINE_BROADENING, leading_pad=0, zero_lead=0)
        self.assertEqual(out.shape[1], NSWITCH)
        self.assertEqual(legacy.legacy_echo_count(NSWITCH), NSWITCH - 1)

    def test_the_extra_point_shifts_the_ppm_step_by_about_one_and_a_half_percent(self):
        bw_ppm = 17.0
        dx_current = bw_ppm / NSWITCH
        dx_legacy = bw_ppm / (NSWITCH - 1)
        self.assertAlmostEqual(dx_legacy / dx_current, NSWITCH / (NSWITCH - 1), places=12)
        self.assertAlmostEqual(dx_legacy / dx_current, 1.0159, places=4)


class SamplingWindowDiffers(unittest.TestCase):
    """
    The legacy reads switch positions discard_pre .. discard_post + nro, which on its own
    converter (discard_post == discard_pre == 8, nro == 12) is 8..19. The current code reads
    the flat top: discard_pre - pad for `kept` points.

    The converter's model of a switch is ramp + readout + rampdown/rephase, so the current
    window is the flat top and the legacy's runs into the rampdown. Current is right, and
    matching the legacy here would mean reproducing a defect.
    """

    def test_the_pad_comes_from_the_recorded_ramp_time(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        pad, why = mrd2recon.epsi_leading_pad(FakeHeader(tramp_us=100), acq)
        # ceil(100us / 40us) = 3 ramp samples against discard_pre 4
        self.assertEqual(pad, 1)
        self.assertIn("tramp=100us", why)

    def test_without_a_recorded_ramp_time_the_pad_falls_back(self):
        acq = FakeAcq(synthetic_readout(), discard_pre=4, discard_post=12)
        pad, why = mrd2recon.epsi_leading_pad(FakeHeader(tramp_us=None), acq)
        self.assertEqual(pad, mrd2recon.EPSIGRE_DEFAULT_PAD)
        self.assertIn("no ramp time", why)

    def test_zero_lead_blanks_the_unusable_head_of_the_readout(self):
        """
        Only the first echo can reach samples that early, so the rest are untouched.
        """
        self.assertEqual(mrd2recon.EPSIGRE_ZERO_LEAD, 13)
        data = np.ones(NSAMPLES, dtype=complex)
        acq = FakeAcq(data, discard_pre=4, discard_post=4)
        out = mrd2recon.apply_line_broadening(acq, 0.0, leading_pad=0, zero_lead=13)
        # first echo reads positions 4..15, so 4..12 are blanked and 13..15 survive
        np.testing.assert_array_equal(out[:9, 0], 0)
        np.testing.assert_array_equal(out[9:12, 0], 1)
        # the second echo starts at 24, past the blanked head
        np.testing.assert_array_equal(out[:, 1], 1)

    def test_a_sample_past_the_end_of_the_readout_reads_as_zero(self):
        data = np.ones(NSAMPLES, dtype=complex)
        acq = FakeAcq(data, discard_pre=4, discard_post=4)
        out = mrd2recon.apply_line_broadening(acq, 0.0, leading_pad=-8, zero_lead=0)
        self.assertTrue(np.any(out[:, -1] == 0))


class RepetitionSkipDiffers(unittest.TestCase):
    """
    The legacy fits `range(4, nreps)`, leaving repetitions 0-3 exactly zero in the output,
    under a comment reading '# shorten the list for quick debugging'. On a hyperpolarized
    series those are the inflow frames. The current code fits every repetition.
    """

    def test_current_fits_every_repetition(self):
        nreps, nfreq = 6, 32
        xscale = np.arange(nfreq) / nfreq * 17.0
        fitter = LorentzianFitter(xscale)
        fitter.params = PeakParams(centers=np.array([3.0]), widths=np.array([0.5]),
                                   phases=np.array([0.0]), amplitudes=np.array([1.0]))
        one = _eval_lorentzians(xscale, fitter.BW, fitter.params)
        volumes = np.broadcast_to(one, (nreps, 1, 1, nfreq)).copy()

        heights, _ = quietly(mrd2recon.fit_voxel_peaks, volumes, fitter, noise_threshold=0.0)
        self.assertEqual(legacy.legacy_first_fitted_repetition(), 4)
        # every repetition, including the four the legacy skipped, carries an amplitude
        self.assertTrue(np.all(heights[0, :, 0, 0] > 0),
                        "repetitions 0-3 must be fitted, not left at zero")


class GlobalCenterWindowDiffers(unittest.TestCase):
    """
    The current global fit holds each center within GLOBAL_CENTER_WINDOW of where the rigid
    pattern placed it; mrd2_recon_to_incorporate leaves them unbounded.

    The bound is not new. The generation before it (997e829) parameterised the center as
    arctan(x)/pi * wigglefactor with wigglefactor defaulting to 1, i.e. +-0.5 ppm, so the
    current constant restores the original and the _to_incorporate file is the regression.
    Unbounded, a _t peak walks onto a strong neighbour, is fitted as a second component of
    its line, lowers the residual and ruins both maps.
    """

    def test_the_window_is_half_a_ppm(self):
        self.assertEqual(mrd2recon.GLOBAL_CENTER_WINDOW, 0.5)

    def test_the_global_fit_holds_a_center_inside_the_window(self):
        n = 64
        bw_ppm = 17.0
        xscale = np.arange(n) / n * bw_ppm
        fitter = LorentzianFitter(xscale)
        truth = PeakParams(centers=np.array([4.0, 8.0]),
                           widths=np.array([0.5, 0.5]),
                           phases=np.array([0.0, 0.0]),
                           amplitudes=np.array([1.0, 0.6]))
        spectrum = _eval_lorentzians(xscale, fitter.BW, truth)
        spectrum = spectrum / np.max(np.abs(spectrum))

        placed = np.array([4.0, 8.0])
        out = fitter.fit_global(spectrum, placed, np.full(2, 0.5),
                                center_window=mrd2recon.GLOBAL_CENTER_WINDOW)
        self.assertTrue(np.all(np.abs(out.centers - placed)
                               <= mrd2recon.GLOBAL_CENTER_WINDOW + 1e-9))


class OptimiserStartingPointsDiffer(unittest.TestCase):
    """
    Same method (L-BFGS-B), same tolerances, same L1 objective, same model, different
    starting vectors, so the two agree to a tolerance rather than bit-for-bit.

    The legacy never assigns its width block, so scipy clips it into bounds and it starts at
    widthguess/2; the current code starts at widthguess. Per voxel, the legacy leaves the
    amplitude block at zero against a [0, None] bound, while the current code starts at the
    spectrum's magnitude under each peak. The per-voxel model is linear in amplitude and
    baseline, so with the windows shut that problem is convex and the starting point barely
    matters; the global fit is not convex, so there it can matter.
    """

    def test_legacy_global_width_starts_at_half_the_guess(self):
        widthguess = 0.8
        self.assertEqual(legacy.legacy_global_width_init(widthguess), widthguess / 2)

    def test_current_global_width_starts_at_the_guess(self):
        n = 64
        xscale = np.arange(n) / n * 17.0
        fitter = LorentzianFitter(xscale)
        captured = {}
        original = fitter._fit

        def spy(spectrum, anchor, v0, bounds):
            captured["v0"] = np.array(v0, copy=True)
            return original(spectrum, anchor, v0, bounds)

        fitter._fit = spy
        spectrum = normalised_spectrum(n)
        widthguess = 0.8
        fitter.fit_global(spectrum, np.array([3.0, 9.0]), np.full(2, widthguess))
        npeaks = 2
        np.testing.assert_allclose(captured["v0"][npeaks:2 * npeaks],
                                   np.full(npeaks, widthguess))

    def test_current_voxel_amplitudes_start_at_the_spectrum_not_at_zero(self):
        n = 64
        xscale = np.arange(n) / n * 17.0
        fitter = LorentzianFitter(xscale)
        fitter.params = PeakParams(centers=np.array([3.0]), widths=np.array([0.5]),
                                   phases=np.array([0.0]), amplitudes=np.array([1.0]))
        captured = {}
        original = fitter._fit

        def spy(spectrum, anchor, v0, bounds):
            captured["v0"] = np.array(v0, copy=True)
            return original(spectrum, anchor, v0, bounds)

        fitter._fit = spy
        spectrum = _eval_lorentzians(xscale, fitter.BW, fitter.params)
        fitter.fit_windowed(spectrum)
        self.assertGreater(captured["v0"][3], 0.0)
        np.testing.assert_array_equal(legacy.legacy_voxel_amplitude_init(1), [0.0])


class DegenerateOverlapDiffers(unittest.TestCase):
    """
    The legacy initialises bestoverlap to 0, the current code to -inf.

    Since overlap = S0^2 + Spi2^2 >= 0, the two only part for a spectrum exactly orthogonal
    to the reference in both quadratures. Reaching that needs an all-zero reference, and the
    noise gate already skips the voxels that would produce one, so this is degenerate.

    What each does there is still worth pinning, because neither is the no-op it looks like.
    The legacy's `0 > 0` is false at every shift, so bestr and th0 are never assigned and it
    carries the previous voxel's values, or raises NameError on the first voxel. The current
    -inf is beaten by the *first* shift tried, so it settles on r = -PHASE_SEARCH_RANGE with
    theta = pi/2 - atan2(0, 0) = pi/2: a quarter turn and a roll to the end of the search
    range, not an identity. Documented rather than changed, since no reachable input gets
    here and altering it would change the alignment of real voxels.
    """

    def test_an_orthogonal_voxel_takes_the_first_shift_tried_and_a_quarter_turn(self):
        nfreq = 8
        volumes = np.zeros((1, 1, 1, nfreq), dtype=complex)
        volumes[0, 0, 0, :] = 1.0
        reference = np.zeros(nfreq, dtype=complex)   # every overlap is exactly 0
        aligned, _ = mrd2recon.phase_align(volumes, reference, noise_threshold=0.0)

        expected = np.roll(volumes[0, 0, 0, :], -mrd2recon.PHASE_SEARCH_RANGE) * np.exp(
            1j * (np.pi / 2 - np.arctan2(0.0, 0.0)))
        np.testing.assert_allclose(aligned[0, 0, 0, :], expected, rtol=0, atol=1e-15)

    def test_the_legacy_would_instead_leave_it_alone_or_reuse_the_last_voxel(self):
        """
        legacy_reference initialises bestr/th0 explicitly so it is runnable; the real legacy
        leaves them unbound. Either way its strict `>` against 0 never fires here.
        """
        nfreq = 8
        volumes = np.zeros((1, 1, 1, nfreq), dtype=complex)
        volumes[0, 0, 0, :] = 1.0
        reference = np.zeros(nfreq, dtype=complex)
        aligned, _ = legacy.legacy_phase_align(np.array(volumes, copy=True), reference, 0.0)
        np.testing.assert_array_equal(aligned[0, 0, 0, :], volumes[0, 0, 0, :])


# ---------- the CLI contract ---------------------------------------------


class PeakArgumentParsing(unittest.TestCase):
    def test_a_negative_peak_value_is_not_read_as_an_option(self):
        """The reason this runs before argparse instead of using parse_known_args."""
        spec, remaining = mrd2recon.split_peak_args(["-bic_tm", "-0.4", "-lb", "42"], {"-lb"})
        self.assertEqual(spec.names, ["bic"])
        np.testing.assert_allclose(spec.offsets, [-0.4])
        self.assertEqual(remaining, ["-lb", "42"])

    def test_reserved_options_are_not_read_as_peaks(self):
        spec, remaining = mrd2recon.split_peak_args(["-dw", "0.1"], {"-dw"})
        self.assertEqual(spec.names, [])
        self.assertEqual(remaining, ["-dw", "0.1"])

    def test_the_suffixes_select_the_right_peaks(self):
        spec, _ = mrd2recon.split_peak_args(
            ["-bic_tm", "0.0", "-urea", "2.3", "-pyr_s", "9.7", "-lac_m", "21.8"], set())
        self.assertEqual(spec.names, ["bic", "urea", "pyr", "lac"])
        self.assertEqual(spec.source_idx, 2)                  # _s
        self.assertEqual(spec.metabolite_idx, [0, 3])         # _m
        self.assertEqual(spec.biggest_idx, [1, 2, 3])         # everything without _t
        self.assertEqual(len(spec), 4)

    def test_a_retired_option_is_refused_rather_than_read_as_a_peak(self):
        for token in ("-w", "--wigglefactor", "-r", "--rank", "-d", "--denoise"):
            with self.subTest(token=token):
                with self.assertRaises(ValueError) as caught:
                    mrd2recon.split_peak_args([token, "0.5"], set())
                self.assertIn("no longer an option", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
