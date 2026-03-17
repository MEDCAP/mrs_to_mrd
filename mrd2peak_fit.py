"""
Peak fitting orchestrator.

Consumes a post-FFT mrd2 stream produced by ``mrd2recon.py`` and emits a
metabolite image stream. The Lorentzian fitting math lives in
``lorentzian_fit.py``; this module is just the pipeline:

    parse header peaks
        -> per-image phase correction (against the recon-time max spectrum)
        -> phantom single-peak fits      (per-voxel amplitude*width scaling)
        -> global multi-peak fit         (locks centers/widths/phases)
        -> per-voxel amplitude fit       (fixed peak shape, free amplitudes)
        -> write metabolites image stream

Modifier conventions on peak names (carried in the header user_parameters):
    _t  tiny peak; not eligible to be the largest peak hypothesis
    _s  source peak (e.g. injected pyruvate)
    _m  derived metabolite
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from typing import BinaryIO, Iterable, Iterator, List, Optional, Tuple

import numpy as np

import mrd

import lorentzian_fit as lf


# ---------- header / stream helpers --------------------------------------


@dataclass
class PeakSpec:
    names: List[str]
    offsets: np.ndarray
    biggest_idx: List[int]
    metabolite_idx: List[int]
    source_idx: Optional[int]
    wigglefactor: float


def _parse_peak_params(header: mrd.Header) -> PeakSpec:
    """Recover the metabolite peak list that ``mrd2recon.append_header`` stuffed
    into ``header.user_parameters``.
    """
    names: List[str] = []
    offsets: List[float] = []
    biggest: List[int] = []
    metabolites: List[int] = []
    source: Optional[int] = None
    wigglefactor = 1.0

    user = getattr(header, "user_parameters", None)
    params = getattr(user, "user_parameter_double", None) or [] if user is not None else []
    for p in params:
        if p.name == "wigglefactor":
            wigglefactor = float(p.value)
            continue
        full = p.name
        modifiers = full[full.find("_"):] if "_" in full else ""
        base = full.split("_", 1)[0]
        idx = len(names)
        names.append(base)
        offsets.append(float(p.value))
        if "t" not in modifiers:
            biggest.append(idx)
        if "m" in modifiers:
            metabolites.append(idx)
        if "s" in modifiers and source is None:
            source = idx

    return PeakSpec(
        names=names,
        offsets=np.asarray(offsets, dtype=float),
        biggest_idx=biggest,
        metabolite_idx=metabolites,
        source_idx=source,
        wigglefactor=float(wigglefactor),
    )

@dataclass
class StreamPayload:
    max_spect: Optional[np.ndarray]
    noise: Optional[float]
    phantom_arrays: List[np.ndarray]
    recon_arrays: List[np.ndarray]
    receiver_bw_hz: Optional[float]


def _read_payload(input: Iterable[mrd.StreamItem]) -> StreamPayload:
    """Split the post-FFT stream into the four buckets we need.

    See ``mrd2recon.kspace_to_ndarray`` for the producer side.
    """
    max_spect: Optional[np.ndarray] = None
    noise: Optional[float] = None
    phantom_arrays: Optional[List[np.ndarray]] = []
    recon_arrays: List[np.ndarray] = []
    receiver_bw_hz: Optional[float] = None

    for item in input:
        if isinstance(item, mrd.StreamItem.NdArrayDouble):
            noise = item.value.data
            continue
        if not isinstance(item, mrd.StreamItem.NdArrayComplexDouble):
            continue
        arr = item.value
        at = arr.head.array_type
        if at == mrd.ArrayType.USER_MAP:
            max_spect = np.asarray(arr.data, dtype=complex)
        elif at == mrd.ArrayType.PHANTOM:
            phantom_arrays.append(np.asarray(arr.data, dtype=complex))
        elif arr.head.image_type == mrd.ArrayImageType.COMPLEX:
            recon_arrays.append(np.asarray(arr.data, dtype=complex))
            if receiver_bw_hz is None:
                receiver_bw_hz = arr.meta.get("receiver bandwidth(Hz)")

    return StreamPayload(
        max_spect=max_spect,
        noise=noise,
        phantom_arrays=phantom_arrays,
        recon_arrays=recon_arrays,
        receiver_bw_hz=receiver_bw_hz,
    )


# ---------- phase correction ---------------------------------------------


def phase_correct_image(
    image: np.ndarray,
    max_spect: np.ndarray,
    *,
    noise_threshold: Optional[float] = None,
    search_range: int = 15,
) -> Tuple[np.ndarray, np.ndarray]:
    """Roll-and-rotate every (j, k) spectrum to maximize overlap with
    ``max_spect``. Returns ``(corrected_image, contribution_to_global_spect)``.

    Mirrors the per-voxel logic in ``mrd2fit_v1.epsi_recon`` (lines 213-231).
    Voxels under ``noise_threshold`` (in absolute magnitude) are left alone
    and excluded from the global accumulation.
    """
    out = np.array(image, dtype=complex, copy=True)
    contrib = np.zeros(out.shape[-1], dtype=complex)
    max_conj = np.conj(max_spect)

    for j in range(out.shape[0]):
        for k in range(out.shape[1]):
            spect = out[j, k, :]
            if noise_threshold is not None and np.max(np.abs(spect)) < noise_threshold:
                continue
            best_overlap = -np.inf
            best_r = 0
            best_th = 0.0
            for r in range(-search_range, search_range + 1):
                rolled = np.roll(spect, r)
                S0 = float(np.sum(np.real(rolled * max_conj)))
                Spi2 = float(np.sum(np.real(rolled * 1j * max_conj)))
                overlap = S0 * S0 + Spi2 * Spi2
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_r = r
                    best_th = np.pi / 2 - np.arctan2(S0, Spi2)
            corrected = np.roll(spect, best_r) * np.exp(1j * best_th)
            out[j, k, :] = corrected
            contrib += corrected

    return out, contrib


# ---------- fitting stages -----------------------------------------------


def _phantom_global_spect(image: np.ndarray) -> np.ndarray:
    """Sum of every voxel spectrum in a phantom image (no phase correction).

    The phantom is a single tube of urea, so simple summation is enough to
    estimate the line shape; matches mrd2fit_v1.epsi_recon's phantom branch.
    """
    return image.reshape(-1, image.shape[-1]).sum(axis=0)


def fit_phantom_arrays(
    phantom_arrays: List[np.ndarray],
    xscale: np.ndarray,
    wigglefactor: float,
) -> Tuple[np.ndarray, float]:
    """Single-peak full Lorentzian fit on every phantom voxel.

    Returns the per-voxel ``A * width`` map and the cross-phantom averaged
    global scaling that callers use to normalize metabolite amplitudes.
    """
    if not phantom_arrays:
        return np.zeros((0, 0, 0)), 1.0

    n_phantom = len(phantom_arrays)
    ny, nx, _ = phantom_arrays[0].shape
    phantom_recon = np.zeros((n_phantom, ny, nx))
    phantom_scaling = 0.0

    for i, image in enumerate(phantom_arrays):
        global_spect = _phantom_global_spect(image)
        scale_g = float(np.max(np.abs(global_spect)))
        if scale_g == 0.0:
            continue
        ng = global_spect / scale_g
        width = lf.estimate_width_from_fwhm(xscale, ng)
        ctx_g = lf.SpectrumCtx.from_axis(xscale, ng, wigglefactor=wigglefactor)
        ref_g = lf.single_peak_ref(xscale, ng, width)
        amp0_g = lf.guess_amplitudes_at(xscale, ng, ref_g.centers)
        gres = lf.fit_full(ctx_g, ref_g, amp0_g)
        phantom_scaling += gres.amplitudes[0] * gres.widths[0] * scale_g / n_phantom

        for j in range(ny):
            for k in range(nx):
                spect = image[j, k, :]
                scale_v = float(np.max(np.abs(spect)))
                if scale_v == 0.0:
                    continue
                nv = spect / scale_v
                ctx_v = lf.SpectrumCtx.from_axis(xscale, nv, wigglefactor=wigglefactor)
                ref_v = lf.single_peak_ref(xscale, nv, width)
                amp0_v = lf.guess_amplitudes_at(xscale, nv, ref_v.centers)
                vres = lf.fit_full(ctx_v, ref_v, amp0_v)
                phantom_recon[i, j, k] = vres.amplitudes[0] * vres.widths[0] * scale_v

    if phantom_scaling == 0.0:
        phantom_scaling = 1.0
    return phantom_recon, phantom_scaling


def fit_global_multipeak(
    global_spect: np.ndarray,
    xscale: np.ndarray,
    spec: PeakSpec,
) -> lf.PeakRef:
    """Try each candidate "biggest peak" hypothesis and return the locked
    centers/widths/phases from the best-fitting one.

    Mirrors mrd2fit_v1.py:296-313.
    """
    scale = float(np.max(np.abs(global_spect)))
    if scale == 0.0:
        raise ValueError("global spectrum is zero; cannot fit peaks")
    norm = global_spect / scale
    bw_ppm = float(np.max(xscale) - np.min(xscale) + (xscale[1] - xscale[0]))
    width = lf.estimate_width_from_fwhm(xscale, norm)
    ctx = lf.SpectrumCtx.from_axis(xscale, norm, wigglefactor=spec.wigglefactor)

    candidates = spec.biggest_idx if spec.biggest_idx else [int(np.argmax(spec.offsets))]

    best_loss = np.inf
    best_result: Optional[lf.FullFitResult] = None
    for icg in candidates:
        ref = lf.multi_peak_ref(xscale, norm, spec.offsets, icg, width, bw_ppm)
        amp0 = lf.guess_amplitudes_at(xscale, norm, ref.centers)
        result = lf.fit_full(ctx, ref, amp0)
        if result.loss < best_loss:
            best_loss = result.loss
            best_result = result

    assert best_result is not None
    return lf.PeakRef(
        centers=best_result.centers,
        widths=best_result.widths,
        phases=best_result.phases,
    )


def fit_voxel_amplitudes(
    hpimgset: np.ndarray,
    ref: lf.PeakRef,
    xscale: np.ndarray,
    *,
    wigglefactor: float = 1.0,
    noise_threshold: float = 0.0,
) -> np.ndarray:
    """Per-voxel amplitude+baseline fit with the global peak shape locked.

    Returns ``metabolites`` with shape ``(npeaks, n_meas, ny, nx)`` to match
    the layout used by ``mrd2fit_v1.generate_epsi_images``.
    """
    npeaks = ref.npeaks
    n_meas, ny, nx, _ = hpimgset.shape
    metabolites = np.zeros((npeaks, n_meas, ny, nx))

    for ide in range(n_meas):
        for j in range(ny):
            for k in range(nx):
                spect = hpimgset[ide, j, k, :]
                scale = float(np.max(np.abs(spect)))
                if scale == 0.0 or scale < noise_threshold:
                    continue
                norm = spect / scale
                ctx = lf.SpectrumCtx.from_axis(xscale, norm, wigglefactor=wigglefactor)
                amp0 = lf.guess_amplitudes_at(xscale, norm, ref.centers)
                res = lf.fit_amplitudes(ctx, ref, amp0)
                metabolites[:, ide, j, k] = res.amplitudes * scale

    return metabolites


# ---------- top-level pipeline -------------------------------------------


@dataclass
class FitOutputs:
    metabolites: np.ndarray
    phantom_recon: np.ndarray
    phantom_scaling: float
    locked_peaks: lf.PeakRef
    xscale: np.ndarray
    global_spect: np.ndarray


def run_pipeline(header: mrd.Header, payload: StreamPayload) -> FitOutputs:
    if payload.max_spect is None:
        raise ValueError("recon stream is missing the USER_MAP max-spectrum")
    if payload.receiver_bw_hz is None:
        raise ValueError("recon stream is missing 'receiver bandwidth(Hz)' meta")
    if not payload.recon_arrays:
        raise ValueError("recon stream contained no reconstructed images")

    spec = _parse_peak_params(header)

    n_freq = payload.max_spect.size
    center_freq_hz = float(header.experimental_conditions.h1resonance_frequency_hz)
    bw_ppm = payload.receiver_bw_hz[0].value / center_freq_hz * 1.0e6
    xscale = np.arange(n_freq) / n_freq * bw_ppm

    noise_threshold = (payload.noise * 3.0) if payload.noise is not None else None

    n_meas = len(payload.recon_arrays)
    ny, nx, nf = payload.recon_arrays[0].shape
    hpimgset = np.zeros((n_meas, ny, nx, nf), dtype=complex)
    global_spect = np.zeros(nf, dtype=complex)
    for i, image in enumerate(payload.recon_arrays):
        corrected, contrib = phase_correct_image(
            image, payload.max_spect, noise_threshold=noise_threshold
        )
        hpimgset[i] = corrected
        global_spect += contrib

    phantom_recon, phantom_scaling = fit_phantom_arrays(
        payload.phantom_arrays, xscale, spec.wigglefactor
    )
    print(f"phantom scaling = {phantom_scaling}", file=sys.stderr)

    locked = fit_global_multipeak(global_spect, xscale, spec)

    metabolites = fit_voxel_amplitudes(
        hpimgset,
        locked,
        xscale,
        wigglefactor=spec.wigglefactor,
        noise_threshold=(payload.noise * 3.0) if payload.noise is not None else 0.0,
    )

    return FitOutputs(
        metabolites=metabolites,
        phantom_recon=phantom_recon,
        phantom_scaling=phantom_scaling,
        locked_peaks=locked,
        xscale=xscale,
        global_spect=global_spect,
    )


# ---------- mrd2 stream output -------------------------------------------


def _generate_metabolite_images(
    header: mrd.Header,
    metabolites: np.ndarray,
    spec: PeakSpec,
) -> Iterator[mrd.StreamItem]:
    """Wrap ``metabolites`` (npeaks, n_meas, ny, nx) into an mrd2 image stream.

    Mirrors the layout of ``mrd2fit_v1.generate_epsi_images``.
    """
    nmet, nimg, _, _ = metabolites.shape
    measfreq = float(header.experimental_conditions.h1resonance_frequency_hz)
    time_between_images_ns = int(3 * 1_000_000_000)

    for ide in range(nimg):
        flags = mrd.ImageFlags(0)
        if ide == 0:
            flags = mrd.ImageFlags.FIRST_IN_SET
        elif ide == nimg - 1:
            flags = mrd.ImageFlags.LAST_IN_SET
        head = mrd.ImageHeader(
            image_type=mrd.ImageType.MAGNITUDE,
            flags=flags,
            measurement_uid=ide,
            measurement_freq=(measfreq + np.uint32(measfreq * spec.offsets / 1e6 + 0.5)),
            measurement_freq_label=np.array(spec.names, dtype=np.dtype(np.object_)),
            repetition=ide,
            acquisition_time_stamp_ns=ide * time_between_images_ns,
            image_index=ide,
            image_series_index=ide,
        )
        # mrd2 image data: (channels, slices, rows, cols, freqs)
        data = np.expand_dims(np.moveaxis(metabolites[:, ide, :, :], 0, 2), (0, 1))
        yield mrd.StreamItem.ImageDouble(mrd.Image(head=head, data=data))


def fit_spectrum(input: BinaryIO, output: Optional[BinaryIO] = None) -> FitOutputs:
    """Read a post-FFT mrd2 stream, fit Lorentzians, optionally write images."""
    with mrd.BinaryMrdReader(input) as reader:
        header = reader.read_header()
        payload = _read_payload(reader.read_data())

    # outputs = run_pipeline(header, payload)

    # if output is not None:
    #     spec = _parse_peak_params(header)
    #     with mrd.BinaryMrdWriter(output) as writer:
    #         writer.write_header(header)
    #         writer.write_data(_generate_metabolite_images(header, outputs.metabolites, spec))

    # return outputs



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit Lorentzian peaks on a reconstructed mrd2 stream.")
    parser.add_argument("-i", "--input", type=str, required=False, help="Input recon mrd2 file")
    parser.add_argument("-o", "--output", type=str, required=False, help="Output fitted mrd2 file")
    args = parser.parse_args()
    in_stream = open(args.input, "rb") if args.input else sys.stdin.buffer
    out_stream = open(args.output, "wb") if args.output else sys.stdout.buffer
    fit_spectrum(in_stream, out_stream)


# now look for specification of metabolite peaks
# BA's cirrhrat is -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8
# SZ's mouse kidney is -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -poop_tm 15.9 -hyd_tm 18.1 -lac_m 21.8
# BA's spectra is -bic_tm -0.4 -urea 2.1 -urea2_t 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8 -w 0.5
# DT's spectra is -urea 0.0 -KIC_s 8.6 -leu_tm 13.0 -hyd_tm 18.1 -?_tm 21.8 -w 1.
