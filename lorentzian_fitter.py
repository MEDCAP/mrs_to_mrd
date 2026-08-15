import numpy as np
from scipy.optimize import minimize
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

    Global fit: optimises centers, widths, phases, amplitudes and baseline over a whole
    spectrum. Windowed fit: refits one voxel with the line shape held at, or within a stated
    distance of, the global result.

    Every parameter is optimised in its own units - a center as a ppm offset from where it was
    placed, a width as an absolute ppm value - and is constrained by an explicit box bound
    rather than by a transform folded into the model.

    Usage:
        fitter = LorentzianFitter(xscale)
        params = fitter.fit_global(spectrum, centers_init, widths_init, width_bounds)
        fitted_spectrum = fitter.eval()
        voxel = fitter.fit_windowed(voxel_spectrum, df=0.1, dw=0.0, dph=0.0)
    """

    def __init__(self, xscale: np.ndarray):
        self.xscale = np.asarray(xscale, dtype=float)
        self.BW = float(self.xscale[-1] - self.xscale[0] + (self.xscale[1] - self.xscale[0]))
        self.params: PeakParams | None = None

    def eval(self) -> np.ndarray:
        """Evaluate the fitted spectrum at xscale points."""
        return _eval_lorentzians(self.xscale, self.BW, self.params)

    def _fit(self,
             spectrum: np.ndarray,
             anchor: np.ndarray,
             v0: np.ndarray,
             bounds: list) -> PeakParams:
        """Minimise the Lorentzian model against a spectrum.

        The optimizer vector holds one block per quantity rather than one block per peak:

            [dc_0..dc_n-1 | w_0..w_n-1 | ph_0..ph_n-1 | A_0..A_n-1 | re_b | im_b]

        where dc is a ppm offset from `anchor` and w is an absolute width in ppm.
        Args:
            - anchor: the centers the dc block is measured from
            - v0: the initial vector, in that layout
            - bounds: one (lo, hi) per entry of v0, either of which may be None
        Returns:
            - the fitted PeakParams, with negative amplitudes folded into their phase
        """
        npeaks = len(anchor)
        anchor = np.asarray(anchor, dtype=float)

        def unpack(v: np.ndarray) -> PeakParams:
            return PeakParams(
                centers=anchor + v[:npeaks],
                widths=np.array(v[npeaks:2 * npeaks]),
                phases=np.array(v[2 * npeaks:3 * npeaks]),
                amplitudes=np.array(v[3 * npeaks:4 * npeaks]),
                baseline=v[4 * npeaks] + 1j * v[4 * npeaks + 1],
            )

        # A parameter whose bounds have collapsed onto a point cannot move, and leaving it in the
        # vector only costs the optimizer a finite difference per iteration to rediscover that.
        # With the voxel windows at their default of zero that is three parameters in every four
        v = np.asarray(v0, dtype=float).copy()
        free = np.ones(len(v), dtype=bool)
        for i, (lo, hi) in enumerate(bounds):
            if lo is not None and hi is not None and lo == hi:
                v[i] = lo
                free[i] = False

        def scatter(v_free: np.ndarray) -> np.ndarray:
            full = v.copy()
            full[free] = v_free
            return full

        def residual(v_free: np.ndarray) -> float:
            model = _eval_lorentzians(self.xscale, self.BW, unpack(scatter(v_free)))
            return float(np.sum(np.abs(model - spectrum)))

        result = minimize(residual, v[free],
                          bounds=[b for b, keep in zip(bounds, free) if keep])

        params = unpack(scatter(result.x))
        neg = params.amplitudes < 0
        params.amplitudes[neg] *= -1
        params.phases[neg] += np.pi
        params.loss = float(result.fun)
        return params

    def fit_global(self,
                   spectrum: np.ndarray,
                   centers_init: np.ndarray,
                   widths_init: np.ndarray,
                   width_bounds: tuple = None) -> PeakParams:
        """Fit all Lorentzian parameters (centers, widths, phases, amplitudes, baseline).

        Centers are free to move anywhere; only the widths are constrained, which is what
        holds the fit together when several peaks are close enough to trade signal.
        Stores the result in self.params and returns it.
        Args:
            - centers_init: where each peak is thought to be, in ppm
            - widths_init: the width guess each peak starts from, in ppm
            - width_bounds: an absolute (lo, hi) width range in ppm, applied to every peak
        """
        npeaks = len(centers_init)
        c0 = np.asarray(centers_init, dtype=float)
        w0 = np.asarray(widths_init, dtype=float)

        nearest = [int(np.argmin(np.abs(self.xscale - c))) for c in c0]
        phases_init = np.array([np.angle(spectrum[i]) for i in nearest])
        amps_init = np.array([np.abs(spectrum[i]) for i in nearest])
        v0 = np.concatenate((np.zeros(npeaks), w0, phases_init, amps_init, [0.0, 0.0]))

        bounds = [(None, None)] * (4 * npeaks + 2)
        for j in range(npeaks):
            # a width is an absolute value here, so an unbounded optimizer can walk one through
            # zero and the model diverges. The default is the range the arctan parameterization
            # this replaces used to enforce implicitly
            bounds[npeaks + j] = (width_bounds if width_bounds is not None
                                  else (0.1 * w0[j], 1.9 * w0[j]))

        self.params = self._fit(spectrum, c0, v0, bounds)
        return self.params

    def fit_windowed(self,
                     spectrum: np.ndarray,
                     *,
                     df: float = 0.0,
                     dw: float = 0.0,
                     dph: float = 0.0) -> PeakParams:
        """Refit one spectrum with its line shape held near the global fit.

        Every peak keeps the center, width and phase fit_global settled on, and may depart from
        each by at most the matching window. At the default of zero windows the line shape is
        pinned outright and only the amplitudes and a complex baseline are free, which is the
        least a voxel too noisy to support a full fit needs to still yield an amplitude.

        Normalises internally, so amplitudes come back in the units of the input spectrum.
        self.params is left alone: it is the anchor every voxel is fitted against.
        Args:
            - df: how far a center may move from the global fit, in ppm
            - dw: how far a width may move from the global fit, in ppm
            - dph: how far a phase may move from the global fit, in radians
        Returns:
            - this spectrum's own PeakParams, not stored on the fitter
        """
        anchor = np.asarray(self.params.centers, dtype=float)
        widths = np.asarray(self.params.widths, dtype=float)
        phases = np.asarray(self.params.phases, dtype=float)
        npeaks = len(anchor)

        scaling = float(np.max(np.abs(spectrum)))
        if scaling == 0.0:
            return PeakParams(centers=anchor.copy(), widths=widths.copy(),
                              phases=phases.copy(), amplitudes=np.zeros(npeaks),
                              baseline=0j, loss=0.0)
        spectrum_norm = spectrum / scaling

        amps_init = np.array([np.abs(spectrum_norm[int(np.argmin(np.abs(self.xscale - c)))])
                              for c in anchor])
        v0 = np.concatenate((np.zeros(npeaks), widths, phases, amps_init, [0.0, 0.0]))

        bounds = [(None, None)] * (4 * npeaks + 2)
        for j in range(npeaks):
            bounds[j] = (-df, df)
            bounds[npeaks + j] = (widths[j] - dw, widths[j] + dw)
            bounds[2 * npeaks + j] = (phases[j] - dph, phases[j] + dph)
            bounds[3 * npeaks + j] = (0.0, None)

        params = self._fit(spectrum_norm, anchor, v0, bounds)
        params.amplitudes = params.amplitudes * scaling
        params.baseline = params.baseline * scaling
        return params
