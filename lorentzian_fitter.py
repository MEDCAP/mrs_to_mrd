import numpy as np
from scipy.optimize import minimize, Bounds
from dataclasses import dataclass


@dataclass
class PeakParams:
    centers: np.ndarray
    widths: np.ndarray
    phases: np.ndarray
    amplitudes: np.ndarray
    baseline: complex = 0j
    loss: float = np.inf    # residual of the fit that produced these, for hypothesis selection


def estimate_width_fwhm(xscale: np.ndarray, spectrum_norm: np.ndarray) -> float:
    """Estimate a peak width from the FWHM of the largest peak in a normalised spectrum.

    Walks outward from the maximum in both directions until the magnitude drops below
    half, wrapping around the spectrum, and halves the span.
    Args:
        - xscale: the frequency axis, evenly spaced
        - spectrum_norm: complex spectrum scaled so max(abs()) == 1
    """
    npts = len(spectrum_norm)
    maxidx = int(np.argmax(np.abs(spectrum_norm)))
    leftidx = -1
    rightidx = -1
    for isp in range(npts):
        if np.abs(spectrum_norm[(maxidx - isp) % npts]) < 0.5 and leftidx == -1:
            leftidx = -isp
        if np.abs(spectrum_norm[(maxidx + isp) % npts]) < 0.5 and rightidx == -1:
            rightidx = isp
    if leftidx == -1:
        leftidx = -npts // 2
    if rightidx == -1:
        rightidx = npts // 2
    return float((rightidx - leftidx) * (xscale[1] - xscale[0]) / 2)


def candidate_centers(xscale: np.ndarray,
                      spectrum_norm: np.ndarray,
                      offsets: np.ndarray,
                      biggest_idx: int,
                      bw_ppm: float) -> np.ndarray:
    """Place the peak centers on the assumption that peak `biggest_idx` is the largest one.

    The peak offsets are a rigid pattern; the only unknown is which of them sits under the
    tallest point in the spectrum. Anchoring the pattern there and wrapping modulo the
    spectral width gives one hypothesis per candidate.
    """
    anchor = xscale[int(np.argmax(np.abs(spectrum_norm)))]
    return (anchor - (np.asarray(offsets, dtype=float) - offsets[biggest_idx])) % bw_ppm


def _eval_lorentzians(xscale: np.ndarray, BW: float, params: PeakParams) -> np.ndarray:
    y = np.full(len(xscale), params.baseline, dtype='complex')
    for j in range(len(params.centers)):
        scaled = params.amplitudes[j] * np.exp(1j * params.phases[j])
        for shift in (0.0, -BW, BW):
            y += scaled / (1 + 1j * (xscale - params.centers[j] + shift) / params.widths[j])
    return y


class LorentzianFitter:
    """Fits Lorentzian peaks to MRS spectra.

    Global fit: optimises centers, widths, phases, amplitudes and baseline.
    Amplitude fit: fixed peak shape from global fit, optimises amplitudes only.

    Usage:
        fitter = LorentzianFitter(xscale, wigglefactor)
        params = fitter.fit_global(spectrum, centers_init, widths_init)
        fitted_spectrum = fitter.eval()
        amplitudes = fitter.fit_amplitudes(voxel_spectrum)
    """

    def __init__(self, xscale: np.ndarray, wigglefactor: float = 1.0):
        self.xscale = np.asarray(xscale, dtype=float)
        self.BW = float(self.xscale[-1] - self.xscale[0] + (self.xscale[1] - self.xscale[0]))
        self.wigglefactor = wigglefactor
        self.params: PeakParams | None = None

    def eval(self) -> np.ndarray:
        """Evaluate the fitted spectrum at xscale points."""
        return _eval_lorentzians(self.xscale, self.BW, self.params)

    def fit_global(self,
                   spectrum: np.ndarray,
                   centers_init: np.ndarray,
                   widths_init: np.ndarray,
                   width_bounds: tuple = None) -> PeakParams:
        """Fit all Lorentzian parameters (centers, widths, phases, amplitudes, baseline).

        The optimizer vector has 4 values per peak followed by 2 baseline values:
          [center_t, width_t, phase, amplitude,  <-- peak 0
           center_t, width_t, phase, amplitude,  <-- peak 1
           ...
           re_baseline, im_baseline]

        center_t and width_t are arctan-transformed to constrain them near their
        initial values. Stores the result in self.params and returns it.

        width_bounds, if given, is an absolute (lo, hi) width range in xscale units;
        it is mapped through the same arctan transform and passed to the optimizer.
        """
        npeaks = len(centers_init)
        c0 = np.asarray(centers_init, dtype=float)
        w0 = np.asarray(widths_init, dtype=float)
        wf = self.wigglefactor

        def _pack_peak(j: int, center: float, width: float, phase: float, amplitude: float) -> np.ndarray:
            return np.array([
                np.tan((center - c0[j]) / wf * np.pi),
                np.tan((width / w0[j] - 1) * np.pi / 1.8),
                phase,
                amplitude,
            ])

        def _unpack_peak(j: int, v4: np.ndarray) -> tuple[float, float, float, float]:
            center    = c0[j] + np.arctan(v4[0]) / np.pi * wf
            width     = w0[j] * (1 + np.arctan(v4[1]) * 1.8 / np.pi)
            phase     = v4[2]
            amplitude = v4[3]
            return center, width, phase, amplitude

        def _pack(params: PeakParams) -> np.ndarray:
            peak_vecs = [_pack_peak(j, params.centers[j], params.widths[j],
                                    params.phases[j], params.amplitudes[j])
                         for j in range(npeaks)]
            return np.concatenate([*peak_vecs, [np.real(params.baseline), np.imag(params.baseline)]])

        def _unpack(v: np.ndarray) -> PeakParams:
            per_peak = [_unpack_peak(j, v[4*j:4*j+4]) for j in range(npeaks)]
            return PeakParams(
                centers=np.array([p[0] for p in per_peak]),
                widths=np.array([p[1] for p in per_peak]),
                phases=np.array([p[2] for p in per_peak]),
                amplitudes=np.array([p[3] for p in per_peak]),
                baseline=v[4*npeaks] + 1j * v[4*npeaks+1],
            )

        def _residual(v: np.ndarray) -> float:
            return float(np.sum(np.abs(_eval_lorentzians(self.xscale, self.BW, _unpack(v)) - spectrum)))

        phases_init = np.array([np.angle(spectrum[np.argmin(np.abs(self.xscale - c))]) for c in c0])
        amps_init = np.array([np.abs(spectrum[np.argmin(np.abs(self.xscale - c))]) for c in c0])
        v0 = _pack(PeakParams(c0, w0, phases_init, amps_init, 0j))

        bounds = None
        if width_bounds is not None:
            lo, hi = width_bounds
            bounds = [(None, None)] * (4 * npeaks + 2)
            for j in range(npeaks):
                # the transform saturates at w0*(1 +- 0.9); clip so tan() stays finite
                t_lo, t_hi = (np.tan(np.clip(w / w0[j] - 1, -0.89, 0.89) * np.pi / 1.8)
                              for w in (lo, hi))
                bounds[4 * j + 1] = (t_lo, t_hi)

        result = minimize(_residual, v0, bounds=bounds)
        params = _unpack(result.x)
        neg = params.amplitudes < 0
        params.amplitudes[neg] *= -1
        params.phases[neg] += np.pi
        params.loss = float(result.fun)

        self.params = params
        return params

    def fit_amplitudes(self, spectrum: np.ndarray) -> np.ndarray:
        """Fit only amplitudes with fixed peak shape (centers/widths/phases) from fit_global.

        Normalises internally so the result is in the same units as the input spectrum.
        Returns an array of amplitudes, one per peak.
        """
        npeaks = len(self.params.centers)
        scaling = np.max(np.abs(spectrum))
        if scaling == 0:
            return np.zeros(npeaks)
        spectrum_norm = spectrum / scaling

        centers = self.params.centers
        widths = self.params.widths
        phases = self.params.phases

        def _residual(v: np.ndarray) -> float:
            p = PeakParams(centers, widths, phases, v[:npeaks], v[npeaks] + 1j * v[npeaks+1])
            return float(np.sum(np.abs(_eval_lorentzians(self.xscale, self.BW, p) - spectrum_norm)))

        amp_init = np.array([np.abs(spectrum_norm[np.argmin(np.abs(self.xscale - c))]) for c in centers])
        v0 = np.concatenate((amp_init, [0.0, 0.0]))
        bounds = Bounds(
            np.concatenate((np.zeros(npeaks), [-0.1, -0.1])),
            np.concatenate((amp_init * 1.5, [0.1, 0.1])),
        )
        return minimize(_residual, v0, bounds=bounds).x[:npeaks] * scaling
