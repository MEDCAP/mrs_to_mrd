"""
The legacy EPSI reconstruction's arithmetic, transcribed stage by stage.

Transcribed from mrd2_recon_to_incorporate.py's epsi_recon (L151-479) and the fitter it
uses, lorn_to_incorporate.py. Line numbers below refer to those two files.

Why a transcription rather than an import. The legacy file cannot be run against anything
this repo produces today, for four independent reasons:

  - it imports acqtypes, which is not on this branch, and lorn, which exists here only as
    lorn_to_incorporate.py
  - it reads acq.data as (samples, coils); MRStomrd2 has written (coils, samples) since
    68a5773, and the mrd schema declares that order
  - it reads the centre frequency from acq.head.acquisition_center_frequency, which nothing
    in this repo writes any more, so its whole ppm axis would divide by zero
  - it calls plt.show() throughout

So the two are not one implementation over one format; they are two implementations over
two formats. Transcribing the arithmetic is the only way to compare the parts that are
meant to be the same, and it has the side benefit of being readable as a record of what
the legacy actually did.

Everything here is deliberately a literal transcription, including the parts that are
wrong. Where the current code differs on purpose, test_epsi_parity.py names the difference
and asserts it. Nothing in this file should be "cleaned up" to agree with the current
code: its whole value is being an independent witness.

numpy only, so this runs anywhere.
"""

import numpy as np


# ---------- k-space assembly ---------------------------------------------


def legacy_echo_count(nswitch: int) -> int:
    """
    How many echoes the legacy keeps: all but the last.

    epsi_recon L174 allocates `(necho-1) * fidpad` and L191 loops
    `for iecho in range(a.head.idx.contrast-1)`, so the final echo is never filled.
    The comment above it reads '# method1 discard the first echo=63 echoes', but the code
    drops the last one.
    """
    return nswitch - 1


def legacy_switch_window(discard_pre: int, discard_post: int, nro: int) -> slice:
    """
    Which samples of a switch the legacy reads, relative to the start of that switch.

    epsi_recon L192:
        a.data[(0 + iecho*totalppswitch + discard_pre):(0 + iecho*totalppswitch + discard_post + nro), 0]

    Note it opens at discard_pre and closes at discard_post + nro, mixing the two. That is
    only self-consistent because the converter it was written against set
    discard_post = discard_pre (68a5773 L65). Under today's converter, which sets
    discard_post = 3 * discard_pre, the slice is the wrong length.
    """
    return slice(discard_pre, discard_post + nro)


def legacy_apply_line_broadening(data: np.ndarray,
                                 nswitch: int,
                                 totalppswitch: int,
                                 discard_pre: int,
                                 discard_post: int,
                                 nro: int,
                                 sample_time_ns: float,
                                 lb: float) -> np.ndarray:
    """
    One readout into (kept points, echoes), apodized. epsi_recon L191-195.

    Args:
        - data: the readout as (samples,), i.e. already unwrapped from the legacy's
          (samples, coils) layout
    """
    necho = legacy_echo_count(nswitch)
    out = np.zeros((nro, necho), dtype='complex')
    window = legacy_switch_window(discard_pre, discard_post, nro)
    for iecho in range(necho):
        base = iecho * totalppswitch
        out[:, iecho] = data[base + window.start:base + window.stop]
        # decay is over switches, not over the points inside one
        tk = iecho * sample_time_ns * totalppswitch / 1.0e+9
        out[:, iecho] *= np.exp(-tk * lb)
    return out


def legacy_fft(kspace: np.ndarray) -> np.ndarray:
    """epsi_recon L219: np.fft.fftshift(np.fft.fftn(kspace)), all axes by numpy default."""
    return np.fft.fftshift(np.fft.fftn(kspace))


def legacy_bw_hz(sample_time_ns: float, totalppswitch: int) -> float:
    """epsi_recon L265: BW = 1 / sampletime / totalppswitch, with sampletime in seconds."""
    sampletime = sample_time_ns / 1.0e+9
    return 1.0 / sampletime / totalppswitch


def legacy_xscale(n: int, bw_hz: float, centerfreq_hz: float) -> np.ndarray:
    """epsi_recon L266: arange(n)/n * BW/centerfreq*1e6, running from 0 rather than centred."""
    return np.array(range(n)) / n * bw_hz / centerfreq_hz * 1.0e+6


# ---------- noise and alignment ------------------------------------------


def legacy_noise(imgset: np.ndarray) -> float:
    """epsi_recon L238: the last repetition of the series has decayed, so it measures noise."""
    return float(np.mean(np.abs(imgset[-1, :, :, :])))


def legacy_phase_align(imgset: np.ndarray,
                       maxspect: np.ndarray,
                       noise: float,
                       search_range: int = 15):
    """
    epsi_recon L241-263. Roll and rotate every voxel onto the reference, and sum.

    Transcribed literally, including `bestoverlap = 0` rather than -inf: if every shift
    scored exactly zero the legacy would carry the previous voxel's bestr/th0, or raise
    NameError on the first voxel. Since overlap = S0^2 + Spi2^2 >= 0, that needs a spectrum
    exactly orthogonal to the reference in both quadratures, which the noise gate already
    skips. Kept as-is so the difference is visible rather than papered over.

    Mutates imgset in place, as the legacy does, and returns the global spectrum.
    """
    globalspect = np.zeros(imgset.shape[3], dtype='complex')
    for ide in range(imgset.shape[0]):
        for j in range(imgset.shape[1]):
            for k in range(imgset.shape[2]):
                thisspect = imgset[ide, j, k, :]
                if np.max(np.abs(thisspect)) < noise * 3:
                    continue
                bestoverlap = 0
                bestr = 0
                th0 = 0.0
                for r in range(-search_range, search_range + 1):
                    thisrollspect = np.roll(thisspect, r)
                    S0 = np.sum(np.real(thisrollspect * np.conj(maxspect)))
                    Spi2 = np.sum(np.real(thisrollspect * 1j * np.conj(maxspect)))
                    overlap = S0**2 + Spi2**2
                    if overlap > bestoverlap:
                        bestr = r
                        bestoverlap = overlap
                        th0 = np.pi / 2 - np.arctan2(S0, Spi2)
                imgset[ide, j, k, :] = np.roll(imgset[ide, j, k, :], bestr) * np.exp(1j * th0)
                globalspect += imgset[ide, j, k, :]
    return imgset, globalspect


# ---------- the fit ------------------------------------------------------


def legacy_width_guess(globalspect_norm: np.ndarray, xscale: np.ndarray) -> float:
    """
    epsi_recon L270-278. FWHM of the largest peak, walking outward with wraparound.

    Note the EPSI path divides the span by 2. spectra_recon (the FID path, L575-578) writes
    the same walk with `maxpeakidx +- isp` and divides by 4; those offsets cancel in the
    difference, so the FID version is exactly half this one. This is the EPSI one.
    """
    maxpeakidx = np.argmax(np.abs(globalspect_norm))
    n = len(globalspect_norm)
    leftidx = -1
    rightidx = -1
    for isp in range(n):
        if np.abs(globalspect_norm[(maxpeakidx - isp) % n]) < 0.5 and leftidx == -1:
            leftidx = -isp
        if np.abs(globalspect_norm[(maxpeakidx + isp) % n]) < 0.5 and rightidx == -1:
            rightidx = isp
    return (rightidx - leftidx) * (xscale[1] - xscale[0]) / 2


def legacy_candidate_centers(globalspect_norm: np.ndarray,
                             xscale: np.ndarray,
                             peakoffsets: np.ndarray,
                             biggest: int,
                             bw_ppm: float) -> np.ndarray:
    """
    epsi_recon L333-334. Anchor the rigid offset pattern at the tallest point, mod the width.
    """
    return (xscale[np.argmax(np.abs(globalspect_norm))]
            - (peakoffsets - peakoffsets[biggest])) % bw_ppm


def legacy_bw_ppm_for_model(xscale: np.ndarray) -> float:
    """lorn_to_incorporate.py L62: BW = max(x) - min(x) + (x[1] - x[0])."""
    return max(xscale) - min(xscale) + (xscale[1] - xscale[0])


def legacy_lorneval(xscale: np.ndarray,
                    bw: float,
                    centers: np.ndarray,
                    widths: np.ndarray,
                    phases: np.ndarray,
                    amplitudes: np.ndarray,
                    baseline: complex) -> np.ndarray:
    """
    lorn_to_incorporate.py L83-98. The line shape, with its +-BW wraparound copies.

    Accumulation order is preserved: for each peak, the term at c, then c-BW, then c+BW.
    Floating point addition is not associative, so that order is part of the result.
    """
    y = np.zeros(len(xscale), dtype='complex') + baseline
    for j in range(len(centers)):
        y += amplitudes[j] * np.exp(1j * phases[j]) / (1 + 1j * (xscale - centers[j]) / widths[j])
        y += amplitudes[j] * np.exp(1j * phases[j]) / (1 + 1j * (xscale - centers[j] - bw) / widths[j])
        y += amplitudes[j] * np.exp(1j * phases[j]) / (1 + 1j * (xscale - centers[j] + bw) / widths[j])
    return y


def legacy_lornfit_loss(model: np.ndarray, spect: np.ndarray) -> float:
    """lorn_to_incorporate.py L101-104: sum(abs(y - spect)), an L1 norm on the complex residual."""
    return float(np.sum(abs(model - spect)))


def legacy_unpack_x0(x0: np.ndarray, anchor_centers: np.ndarray):
    """
    lorn_to_incorporate.py L16-32. The optimiser vector, blocked by quantity.

        [dc_0..dc_n-1 | w_0..w_n-1 | ph_0..ph_n-1 | A_0..A_n-1 | re_b | im_b]

    The center block is a *delta* from the stored centers (`c = centers + x0[:npeaks]`),
    which is why the per-voxel bound of +-fit_df bounds the departure from the global
    center rather than from zero.
    """
    npeaks = int(len(x0) / 4)
    c = anchor_centers + x0[:npeaks]
    w = x0[npeaks:2 * npeaks]
    ph = x0[2 * npeaks:3 * npeaks]
    A = x0[3 * npeaks:4 * npeaks]
    b = x0[4 * npeaks] + 1j * x0[4 * npeaks + 1]
    return c, w, ph, A, b


def legacy_fold_negative_amplitudes(amplitudes: np.ndarray, phases: np.ndarray):
    """epsi_recon L343-346: a negative amplitude becomes positive with pi added to its phase."""
    amplitudes = np.array(amplitudes, dtype=float)
    phases = np.array(phases, dtype=float)
    neg = amplitudes < 0
    amplitudes[neg] *= -1
    phases[neg] += np.pi
    return amplitudes, phases


def legacy_global_width_bounds(widthguess: float):
    """epsi_recon L339: bnds[npeaks + ip] = [widthguess / 2, widthguess * 1.5]."""
    return (widthguess / 2, widthguess * 1.5)


def legacy_global_width_init(widthguess: float) -> float:
    """
    What the legacy global fit actually starts its widths at.

    epsi_recon L329 allocates `x0 = np.zeros((4 * npeaks) + 2)` and then assigns only the
    phase block (L338) and the amplitude block (L337). The width block is never written, so
    it starts at 0 and scipy's L-BFGS-B clips it into bounds, landing at widthguess / 2.

    lornputpeakparams(centers, ones*widthguess, ...) at L340 does set the lorn module's
    `widths` global, but lorneval reads its width from x0, not from that global -- the
    global is only read by lor1fit/lor1eval, the amplitudes-only model the EPSI path never
    uses. So the width guess never reaches the starting vector.
    """
    return legacy_global_width_bounds(widthguess)[0]


def legacy_voxel_amplitude_init(npeaks: int) -> np.ndarray:
    """
    epsi_recon L399-406: x0 = zeros(4n+2) with only the phase block assigned, and the
    amplitude block bounded [0, None]. So every voxel amplitude starts at its lower bound.
    """
    return np.zeros(npeaks)


def legacy_first_fitted_repetition() -> int:
    """
    epsi_recon L383: `for ide in range(4, hpimgset.shape[0])`.

    Repetitions 0-3 are never fitted and stay exactly zero in the output. The comment two
    lines above is '# shorten the list for quick debugging', so this is debug residue, and
    on a hyperpolarized series those are the inflow frames.
    """
    return 4


def legacy_voxel_outputs(amplitudes: np.ndarray, widths: np.ndarray, scaling: float):
    """
    epsi_recon L415-416. Peak height and peak area, both back in the input's units.

        metabolites [:, ide, j, k] = np.abs(A) * scaling
        metabolites2[:, ide, j, k] = np.abs(A * w) * scaling

    metabolites2 is computed and plotted but never written to the stream or the .mat.
    """
    return np.abs(amplitudes) * scaling, np.abs(amplitudes * widths) * scaling
