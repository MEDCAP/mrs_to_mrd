"""
Reconstruct and denoise FID spectral time-series data in mrd2 format, with an
interactive tool to compare spectra with/without truncated-SVD denoising and at
different ranks.

Denoising follows Francischello et al., NMR Biomed. 2021;34:e4285: the complex
time-domain FIDs of the acquisition time series are stacked as the columns of a
matrix, and a low-rank approximation (truncated SVD) recovers the signal subspace
while discarding the noise subspace.

Usage:
    # interactive: sliders for rank and time point, raw vs denoised overlay
    python fid_recon.py -i raw_ndarray.mrd2

    # static comparison of raw + several ranks for one time point
    python fid_recon.py -i raw_ndarray.mrd2 --ranks 2 4 9 --repetition 40

    # save the interactive/static figure to a file (headless-friendly)
    python fid_recon.py -i raw_ndarray.mrd2 --save compare.png
"""

import argparse
import sys
from pathlib import Path
from typing import BinaryIO, Union, Iterable, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

import mrd


def svd_decompose(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Economy complex SVD of the FID matrix (computed once, reused for every rank)."""
    return np.linalg.svd(M, full_matrices=False)


def gavish_donoho_rank(S: np.ndarray, shape: Tuple[int, int]) -> int:
    """
    Gavish-Donoho optimal hard threshold for a matrix with unknown noise level:
    thresh = omega(beta) * median(S), beta = min(m, n) / max(m, n).
    """
    m, n = shape
    beta = min(m, n) / max(m, n)
    omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43
    thresh = omega * np.median(S)
    return max(1, int(np.count_nonzero(S > thresh)))


def low_rank(U: np.ndarray, S: np.ndarray, Vh: np.ndarray, rank: int) -> np.ndarray:
    """Reconstruct the rank-`rank` approximation M_hat = U[:, :r] diag(S[:r]) Vh[:r]."""
    rank = int(np.clip(rank, 1, len(S)))
    return (U[:, :rank] * S[:rank]) @ Vh[:rank]


def denoise_svd(M: np.ndarray, rank: int = None) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Truncated-SVD (low-rank) denoising of a complex signal matrix M whose columns are
    the complex time-domain FIDs of a spectral time series.

    Args:
        - M: complex ndarray of shape (nsamples, nspectra); each column is a FID
        - rank: number of singular values to retain. If None, estimated via the
          Gavish-Donoho optimal hard threshold.
    Returns:
        - (M_hat, S, rank): denoised matrix, the singular values, and the rank used
    """
    U, S, Vh = svd_decompose(M)
    if rank is None:
        rank = gavish_donoho_rank(S, M.shape)
    rank = int(np.clip(rank, 1, len(S)))
    return low_rank(U, S, Vh, rank), S, rank


def build_fid_matrix(head: mrd.Header,
                     acq_stream: Iterable[mrd.Acquisition]) -> Tuple[np.ndarray, List[int], int]:
    """
    Stack the complex FIDs into a matrix, placing each FID at the column given by its
    acquisition header index (`idx.repetition`) rather than by arrival order, so the
    result is robust to out-of-order/interleaved streams. Navigation/phantom
    acquisitions are skipped. Unfilled columns are compacted out before returning.

    Assumes single-channel MRS data: acq.data has shape (coils=1, samples).

    Returns:
        - M: complex128 matrix of shape (nsamples, nreps_present)
        - present_reps: the repetition index of each column (in matrix column order)
        - sample_time_ns: dwell/sample time carried on the acquisitions
    """
    acqs = [a for a in acq_stream
            if not (a.head.flags & int(mrd.AcquisitionFlags.IS_NAVIGATION_DATA))]
    if not acqs:
        raise ValueError("No (non-navigation) acquisitions found in the stream")

    nsamples = acqs[0].samples()
    rep_limit = head.encoding[0].encoding_limits.repetition
    max_rep = max(a.head.idx.repetition for a in acqs)
    nreps = (rep_limit.maximum + 1) if rep_limit is not None else (max_rep + 1)

    M = np.zeros((nsamples, nreps), dtype=np.complex128)
    filled = np.zeros(nreps, dtype=bool)
    sample_time_ns = acqs[0].head.sample_time_ns
    for a in acqs:
        if a.samples() != nsamples:
            raise ValueError(f"Acquisition has {a.samples()} samples, expected {nsamples}")
        j = a.head.idx.repetition
        M[:, j] = np.asarray(a.data).squeeze().astype(np.complex128)
        filled[j] = True

    if not filled.all():
        missing = np.where(~filled)[0]
        print(f"Warning: {len(missing)} repetition column(s) not filled ({missing.tolist()}); "
              f"compacting to {int(filled.sum())} present columns", file=sys.stderr)
    present_reps = np.where(filled)[0].tolist()
    return M[:, filled], present_reps, sample_time_ns


def frequency_axis(nsamples: int, sample_time_ns: float, carrier_hz: float):
    """Return (xaxis, label) in ppm if a carrier frequency is available, else in Hz."""
    dwell_s = sample_time_ns / 1e9
    freq_hz = np.fft.fftshift(np.fft.fftfreq(nsamples, d=dwell_s))
    carrier_mhz = (carrier_hz or 0) / 1e6
    if carrier_mhz > 0:
        return freq_hz / carrier_mhz, "Frequency (ppm)"
    return freq_hz, "Frequency (Hz)"


def fids_to_spectra(M: np.ndarray, sample_time_ns: float, line_broadening: float) -> np.ndarray:
    """Apodize each FID column with exp(-pi*LB*t) and FFT to spectra (fftshift-centered)."""
    nsamples = M.shape[0]
    t = np.arange(nsamples) * (sample_time_ns / 1e9)
    apod = np.exp(-np.pi * line_broadening * t)[:, None]
    return np.fft.fftshift(np.fft.fft(M * apod, axis=0), axes=0)


def _project(spectrum: np.ndarray, mode: str) -> np.ndarray:
    if mode == "real":
        return np.real(spectrum)
    if mode == "imag":
        return np.imag(spectrum)
    return np.abs(spectrum)


def compare_interactive(M: np.ndarray, present_reps: List[int], sample_time_ns: float,
                        carrier_hz: float, line_broadening: float, init_rank: int = None,
                        init_rep: int = 0, save: Path = None):
    """
    Interactive comparison: sliders for time point (repetition) and rank, overlaying the
    raw spectrum against the SVD-denoised spectrum, with a singular-value panel showing
    where the current rank truncates.
    """
    U, S, Vh = svd_decompose(M)
    auto_rank = gavish_donoho_rank(S, M.shape)
    max_rank = len(S)
    rank = auto_rank if init_rank is None else int(np.clip(init_rank, 1, max_rank))

    ncols = M.shape[1]
    init_rep = int(np.clip(init_rep, 0, ncols - 1))
    xaxis, xlabel = frequency_axis(M.shape[0], sample_time_ns, carrier_hz)
    ppm = xlabel.endswith("(ppm)")

    raw_spectra = fids_to_spectra(M, sample_time_ns, line_broadening)

    def denoised_spectra(r):
        return fids_to_spectra(low_rank(U, S, Vh, r), sample_time_ns, line_broadening)

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_axes([0.08, 0.34, 0.66, 0.58])   # spectra overlay
    ax_sv = fig.add_axes([0.80, 0.34, 0.17, 0.58])  # singular values

    state = {"rank": rank, "col": init_rep, "mode": "magnitude"}
    den = denoised_spectra(state["rank"])

    (line_raw,) = ax.plot(xaxis, _project(raw_spectra[:, init_rep], state["mode"]),
                          color="0.6", lw=1.0, label="raw")
    (line_den,) = ax.plot(xaxis, _project(den[:, init_rep], state["mode"]),
                          color="C1", lw=1.4, label="denoised")
    if ppm:
        ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Amplitude")
    ax.legend(loc="upper right")

    # singular-value panel
    ax_sv.semilogy(np.arange(1, max_rank + 1), S, "o-", ms=3, color="0.4")
    sv_cut = ax_sv.axvline(state["rank"] + 0.5, color="C3", lw=1.2)
    ax_sv.set_title("singular values", fontsize=9)
    ax_sv.set_xlabel("index")

    def title():
        rep = present_reps[state["col"]]
        ax.set_title(f"rep {rep}  (col {state['col']+1}/{ncols})   "
                     f"rank = {state['rank']} (auto = {auto_rank})   LB = {line_broadening} Hz")

    title()

    # widgets
    ax_rep = fig.add_axes([0.08, 0.18, 0.66, 0.03])
    ax_rank = fig.add_axes([0.08, 0.12, 0.66, 0.03])
    ax_mode = fig.add_axes([0.80, 0.08, 0.17, 0.16])
    s_rep = Slider(ax_rep, "time point", 0, ncols - 1, valinit=init_rep, valstep=1)
    s_rank = Slider(ax_rank, "rank", 1, max_rank, valinit=state["rank"], valstep=1)
    r_mode = RadioButtons(ax_mode, ("magnitude", "real", "imag"), active=0)

    def redraw(rescale=True):
        col, mode = state["col"], state["mode"]
        line_raw.set_ydata(_project(raw_spectra[:, col], mode))
        line_den.set_ydata(_project(den[:, col], mode))
        title()
        if rescale:
            ax.relim()
            ax.autoscale_view(scalex=False)
        fig.canvas.draw_idle()

    def on_rep(v):
        state["col"] = int(v)
        redraw()

    def on_rank(v):
        nonlocal den
        state["rank"] = int(v)
        den = denoised_spectra(state["rank"])
        sv_cut.set_xdata([state["rank"] + 0.5, state["rank"] + 0.5])
        redraw()

    def on_mode(label):
        state["mode"] = label
        redraw()

    s_rep.on_changed(on_rep)
    s_rank.on_changed(on_rank)
    r_mode.on_clicked(on_mode)

    print(f"Interactive compare: matrix={M.shape}, auto rank={auto_rank}, "
          f"top singular values={np.round(S[:min(10, max_rank)], 4)}", file=sys.stderr)

    if save:
        fig.savefig(save, dpi=150)
        print(f"Saved figure to {save}", file=sys.stderr)
    plt.show()


def compare_static(M: np.ndarray, present_reps: List[int], sample_time_ns: float,
                   carrier_hz: float, line_broadening: float, ranks: List[int],
                   rep: int = 0, save: Path = None):
    """Static overlay of the raw spectrum and the denoised spectrum at each given rank,
    for a single time point."""
    U, S, Vh = svd_decompose(M)
    auto_rank = gavish_donoho_rank(S, M.shape)
    col = int(np.clip(rep, 0, M.shape[1] - 1))
    xaxis, xlabel = frequency_axis(M.shape[0], sample_time_ns, carrier_hz)
    ppm = xlabel.endswith("(ppm)")

    raw = fids_to_spectra(M, sample_time_ns, line_broadening)[:, col]

    fig, (ax, ax_sv) = plt.subplots(1, 2, figsize=(12, 6),
                                    gridspec_kw={"width_ratios": [3, 1]})
    ax.plot(xaxis, np.abs(raw), color="0.6", lw=1.0, label="raw")
    for r in ranks:
        den = fids_to_spectra(low_rank(U, S, Vh, r), sample_time_ns, line_broadening)[:, col]
        ax.plot(xaxis, np.abs(den), lw=1.4, label=f"rank {int(r)}")
    if ppm:
        ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Magnitude")
    ax.set_title(f"rep {present_reps[col]}  (col {col+1}/{M.shape[1]})   "
                 f"auto rank = {auto_rank}   LB = {line_broadening} Hz")
    ax.legend(loc="upper right")

    ax_sv.semilogy(np.arange(1, len(S) + 1), S, "o-", ms=3, color="0.4")
    for r in ranks:
        ax_sv.axvline(int(r) + 0.5, lw=1.0, alpha=0.7)
    ax_sv.set_title("singular values", fontsize=9)
    ax_sv.set_xlabel("index")
    fig.tight_layout()

    print(f"Static compare: matrix={M.shape}, ranks={ranks}, auto rank={auto_rank}, "
          f"top singular values={np.round(S[:min(10, len(S))], 4)}", file=sys.stderr)
    if save:
        fig.savefig(save, dpi=150)
        print(f"Saved figure to {save}", file=sys.stderr)
    plt.show()


def acquisition_reader(input: Iterable[mrd.StreamItem]) -> Iterable[mrd.Acquisition]:
    """Yield acquisitions from an mrd stream."""
    for item in input:
        if isinstance(item, mrd.StreamItem.Acquisition):
            yield item.value


def reconstruct_from_mrd(input: Union[str, BinaryIO], line_broadening: float,
                         rank: int, ranks: List[int], repetition: int, save: Path):
    with mrd.BinaryMrdReader(input) as reader:
        header = reader.read_header()
        M, present_reps, sample_time_ns = build_fid_matrix(
            header, acquisition_reader(reader.read_data()))
    carrier_hz = header.experimental_conditions.h1resonance_frequency_hz
    if ranks:
        compare_static(M, present_reps, sample_time_ns, carrier_hz,
                       line_broadening, ranks, repetition, save)
    else:
        compare_interactive(M, present_reps, sample_time_ns, carrier_hz,
                            line_broadening, rank, repetition, save)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare FID spectra with/without truncated-SVD denoising")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input mrd2 file")
    parser.add_argument("-lb", "--line-broadening", type=float, default=42, help="Line broadening factor in Hz (default: 42)")
    parser.add_argument("-r", "--rank", type=int, default=None, help="Initial rank for the interactive view (default: auto via Gavish-Donoho)")
    parser.add_argument("--ranks", type=int, nargs="+", default=None, help="Static comparison: overlay raw + denoised at each of these ranks")
    parser.add_argument("--repetition", type=int, default=0, help="Time-point (column) index to display initially (default: 0)")
    parser.add_argument("--save", type=Path, default=None, help="Save the figure to this path")
    args = parser.parse_args()

    reconstruct_from_mrd(open(args.input, "rb"), args.line_broadening,
                         args.rank, args.ranks, args.repetition, args.save)
