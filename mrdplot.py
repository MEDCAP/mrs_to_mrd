"""
Display the contents of an mrd2 file.

Reads whatever a file happens to hold and draws the figures that apply to it, so the same
command works on a reconstruction, on a raw conversion, and on the legacy image-based output:

    python mrdplot.py -i c13mouse_epsigre_27927_recon.mrd2

mrd2recon.py writes its results as NdArrays tagged with a description, which is what picks the
figure here:

    global_spect / global_spect_fit / lorentzian_*  the fitted spectrum and its peaks
    metabolite_amplitude, metabolite_area          maps, drawn as a peak x repetition montage
    epsi_image                                     the reconstructed volumes
    phantom peak area, singular_values, noise      calibration and diagnostics

An EPSI reconstruction has one global spectrum and a metabolite map; a single voxel FID series
has one spectrum per time point and no map, so the two are told apart by what is present rather
than by a flag.

With --save the figures are written as PNGs instead of shown, for running without a display.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import BinaryIO, Dict, List

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

import mrd

ZOOM_FACTOR = 2         # image interpolation for the metabolite montage
PEAK_COLORS = ['r', 'b', 'g', 'c', 'm', 'y', 'k']


# ---------- reading ------------------------------------------------------


def meta_values(arr: mrd.NdArray, key: str) -> list:
    """Every value stored under a meta key, unwrapped from its ArrayMetaValue union."""
    entries = arr.meta.get(key)
    return [] if entries is None else [entry.value for entry in entries]


def meta_value(arr: mrd.NdArray, key: str, default=None):
    """The single value stored under a meta key."""
    values = meta_values(arr, key)
    return values[0] if values else default


def describe(arr: mrd.NdArray) -> str:
    """The description an NdArray was tagged with, which is what identifies it."""
    return meta_value(arr, "description", "")


class Contents:
    """Everything a file holds, sorted into the buckets the figures below draw from."""

    def __init__(self):
        self.header = None
        self.acquisitions: List[mrd.Acquisition] = []
        self.images: List[mrd.Image] = []
        self.arrays: Dict[str, List[mrd.NdArray]] = {}

    def first(self, description: str):
        """The first array with this description, or None."""
        found = self.arrays.get(description)
        return found[0] if found else None


def read_contents(input: BinaryIO) -> Contents:
    contents = Contents()
    with mrd.BinaryMrdReader(input) as reader:
        contents.header = reader.read_header()
        for item in reader.read_data():
            value = item.value
            if isinstance(value, mrd.Acquisition):
                contents.acquisitions.append(value)
            elif isinstance(value, mrd.Image):
                contents.images.append(value)
            elif isinstance(value, mrd.NdArray):
                contents.arrays.setdefault(describe(value), []).append(value)
    contents.acquisitions.sort(key=lambda a: a.head.acquisition_time_stamp_ns)
    return contents


# ---------- figures ------------------------------------------------------


def plot_fitted_spectrum(contents: Contents, filename: str) -> bool:
    """
    The summed spectrum, the Lorentzian model fitted to it, and where the peaks landed.

    This is the figure the fit is judged by: if the model does not sit on the data, or a peak
    center is not on a peak, the metabolite maps below are not worth reading.
    """
    spectra = contents.arrays.get("global_spect")
    if not spectra or len(spectra) > 1:
        return False    # no spectrum, or a time series, which plot_spectral_series draws

    spect = spectra[0]
    xscale = np.array(meta_values(spect, "xscale_ppm"))
    if xscale.size != spect.data.size:
        xscale = np.arange(spect.data.size)

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(f'fitted spectrum: {filename}')
    plt.plot(xscale, np.real(spect.data), 'r', label='real')
    plt.plot(xscale, np.imag(spect.data), 'g', label='imag')
    plt.plot(xscale, np.abs(spect.data), color='0.7', label='magnitude')

    fit = contents.first("global_spect_fit")
    if fit is not None:
        plt.plot(xscale, np.real(fit.data), 'k', label='fit')
        plt.plot(xscale, np.imag(fit.data), 'k', linestyle='--', label='_nolegend_')

    centers = contents.first("lorentzian_centers_ppm")
    if centers is not None:
        names = meta_values(centers, "peak_names")
        widths = contents.first("lorentzian_widths_ppm")
        top = np.max(np.abs(spect.data))
        # peaks can sit close together, so the labels run vertically and alternate height
        for ip, center in enumerate(centers.data):
            plt.axvline(center, color='0.4', linewidth=0.8)
            label = names[ip] if ip < len(names) else str(ip)
            if widths is not None:
                label += f' {center:.2f}+-{widths.data[ip]:.2f}'
            plt.text(center, top * (0.35 + 0.30 * (ip % 2)), label,
                     rotation=90, va='bottom', ha='right', fontsize=7,
                     bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.5))

    anchor = meta_value(spect, "biggest_peak_name")
    subtitle = f'anchored on {anchor}' if anchor else ''
    loss = meta_value(fit, "fit_loss") if fit is not None else None
    if loss is not None:
        subtitle += f', residual {loss:.3f}'
    plt.title(subtitle)
    plt.xlabel('frequency (ppm)')
    plt.ylabel('amplitude')
    # outside the axes, so it cannot sit on top of a peak
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    plt.tight_layout()
    return True


def plot_spectral_series(contents: Contents, filename: str) -> bool:
    """
    A single voxel FID series, as a stack of spectra and as a time-frequency image.

    The stack shows the line shape and the image shows how it evolves, which is what the
    metabolic flux experiment is actually about.
    """
    spectra = contents.arrays.get("global_spect")
    if not spectra or len(spectra) < 2:
        return False

    data = np.stack([s.data for s in spectra])
    xscale = np.array(meta_values(spectra[0], "xscale_ppm"))
    if xscale.size != data.shape[1]:
        xscale = np.arange(data.shape[1])

    fig, (stacked, image) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f'spectral time series: {filename}')

    offset = np.max(np.abs(data)) * 0.35
    for i, spect in enumerate(data):
        stacked.plot(xscale, np.abs(spect) + i * offset, color=plt.cm.viridis(i / len(data)),
                     linewidth=0.8)
    stacked.set_xlabel('frequency (ppm)')
    stacked.set_ylabel('time point')
    stacked.set_yticks([])
    stacked.set_title(f'{len(data)} spectra')

    handle = image.imshow(np.abs(data), aspect='auto', origin='lower', cmap='viridis',
                          extent=[xscale[0], xscale[-1], 0, len(data)])
    image.set_xlabel('frequency (ppm)')
    image.set_ylabel('time point')
    image.set_title('magnitude')
    fig.colorbar(handle, ax=image)

    # the peaks occupy a small part of the spectral width, so show that part rather than
    # leaving the interesting structure a few pixels wide
    occupied = np.abs(data).sum(axis=0)
    inside = np.flatnonzero(occupied > 0.01 * occupied.max())
    if inside.size:
        pad = max(1, int(0.15 * (inside[-1] - inside[0] + 1)))
        lo = xscale[max(0, inside[0] - pad)]
        hi = xscale[min(len(xscale) - 1, inside[-1] + pad)]
        if hi > lo:
            stacked.set_xlim(lo, hi)
            image.set_xlim(lo, hi)
    return True


def montage(maps: np.ndarray, zoom_factor: int = ZOOM_FACTOR) -> np.ndarray:
    """
    Lay (npeaks, nreps, ny, nx) out as one image, peaks down and repetitions across.

    Each peak is scaled by its own maximum, so a weak metabolite is still visible next to the
    substrate, and a bright line separates the rows.
    """
    npeaks, nreps, ny, nx = maps.shape
    height, width = ny * zoom_factor, nx * zoom_factor
    out = np.zeros((npeaks * height, nreps * width))
    for ipeak in range(npeaks):
        peak_max = np.max(np.abs(maps[ipeak]))
        if peak_max == 0:
            continue
        for irep in range(nreps):
            tile = zoom(maps[ipeak, irep], zoom_factor, order=2) / peak_max
            out[ipeak * height:(ipeak + 1) * height, irep * width:(irep + 1) * width] = tile
        out[ipeak * height, :] = 1
    return out


def plot_metabolite_maps(contents: Contents, filename: str) -> bool:
    """Every metabolite map in the file, one montage per map."""
    drawn = False
    for description in ("metabolite_amplitude", "metabolite_area"):
        arr = contents.first(description)
        if arr is None:
            continue
        maps = np.asarray(arr.data)
        if maps.ndim != 4:
            continue
        names = meta_values(arr, "peak_names")
        npeaks, nreps, ny, nx = maps.shape

        tiles = montage(maps)
        # voxels stay square, so the figure is sized to the montage rather than the other way
        # round, which is what keeps a wide, short montage from floating in empty space
        width_in = float(np.clip(nreps * 0.7, 6, 16))
        fig = plt.figure(figsize=(width_in, width_in * tiles.shape[0] / tiles.shape[1] + 1.2))
        fig.suptitle(f'{description}: {filename}')
        plt.imshow(tiles, cmap='gray')
        plt.xlabel('repetition')
        plt.xticks(np.arange(nreps) * nx * ZOOM_FACTOR + nx * ZOOM_FACTOR / 2,
                   [str(i) for i in range(nreps)], fontsize=6)
        plt.yticks(np.arange(npeaks) * ny * ZOOM_FACTOR + ny * ZOOM_FACTOR / 2,
                   [names[i] if i < len(names) else str(i) for i in range(npeaks)])
        plt.title('each row scaled to its own maximum', fontsize=8)
        plt.tight_layout()
        drawn = True
    return drawn


def plot_time_courses(contents: Contents, filename: str) -> bool:
    """
    Each metabolite summed over the voxels, against repetition.

    The substrate is drawn solid and the metabolites marked _m dashed, so a product rising as
    the substrate decays is visible without fitting anything.
    """
    arr = contents.first("metabolite_amplitude")
    if arr is None or np.asarray(arr.data).ndim != 4:
        return False

    maps = np.asarray(arr.data)
    names = meta_values(arr, "peak_names")
    source = meta_value(arr, "source_peak_index")
    metabolites = set(meta_values(arr, "metabolite_indices"))
    courses = maps.reshape(maps.shape[0], maps.shape[1], -1).sum(axis=2)
    scale = np.max(courses) or 1.0

    fig = plt.figure(figsize=(9, 5))
    fig.suptitle(f'metabolite time courses: {filename}')
    for ip in range(courses.shape[0]):
        name = names[ip] if ip < len(names) else str(ip)
        color = PEAK_COLORS[ip % len(PEAK_COLORS)]
        style = '-' if ip == source else ('--' if ip in metabolites else ':')
        plt.plot(courses[ip] / scale, color + style, label=name, marker='.')
    plt.xlabel('repetition')
    plt.ylabel('summed amplitude (normalized)')
    plt.legend(fontsize=8)
    return True


def plot_diagnostics(contents: Contents, filename: str) -> bool:
    """The phantom map and the SVD singular values, when the run produced them."""
    drawn = False

    phantom = contents.first("phantom peak area")
    if phantom is not None:
        maps = np.asarray(phantom.data)
        if maps.ndim == 3:
            maps = maps[np.newaxis, ...]
        fig = plt.figure(figsize=(7, 4))
        scaling = meta_value(phantom, "phantom_scaling")
        fig.suptitle(f'phantom peak area (scaling {scaling:.3f}): {filename}'
                     if scaling is not None else f'phantom peak area: {filename}')
        plt.imshow(montage(maps), cmap='gray')
        plt.xticks([])
        plt.yticks([])
        drawn = True

    singular = contents.arrays.get("singular_values")
    if singular:
        fig = plt.figure(figsize=(7, 4))
        fig.suptitle(f'SVD singular values: {filename}')
        for arr in singular:
            rank = meta_value(arr, "svd_rank")
            plt.semilogy(np.asarray(arr.data), '.-', linewidth=0.8,
                         color='0.6' if len(singular) > 1 else 'b')
            if rank is not None:
                plt.axvline(rank - 0.5, color='r', linewidth=0.8)
        plt.xlabel('index')
        plt.ylabel('singular value')
        plt.title('red line marks the retained rank')
        drawn = True

    return drawn


def plot_legacy_images(contents: Contents, filename: str) -> bool:
    """
    The mrd.Image output the legacy reconstruction wrote, kept so old files still open.

    The legacy recon streamed screenshots of its own figures as uint32 images with ARGB packed
    one pixel per word, and metabolite images as everything else. There is no image type that
    says "bitmap" - it wrote ImageType.COMPLEX for both - so the dtype is what tells them apart.
    """
    if not contents.images:
        return False

    bitmaps = [i for i in contents.images if np.asarray(i.data).dtype == np.uint32]
    metabolite_images = [i for i in contents.images if np.asarray(i.data).dtype != np.uint32]

    for image in bitmaps:
        packed = np.squeeze(image.data)
        unpacked = np.zeros(packed.shape[:2] + (4,), dtype=np.uint8)
        for channel, shift in enumerate((0, 8, 16, 24)):
            unpacked[:, :, channel] = ((packed >> shift) & 0xFF).astype(np.uint8)
        fig = plt.figure()
        fig.suptitle(f'File: {filename}')
        plt.xticks([])
        plt.yticks([])
        plt.imshow(unpacked)

    if metabolite_images:
        # each image is (rows, cols, metabolites) once the leading singleton axes are dropped
        stacked = np.stack([np.squeeze(i.data) for i in metabolite_images])
        if stacked.ndim == 4:
            fig = plt.figure()
            fig.suptitle(f'metabolic images File: {filename}')
            plt.imshow(montage(np.moveaxis(stacked, 3, 0)), cmap='gray')
            plt.xticks([])
            plt.yticks([])

    return True


def plot_acquisitions(contents: Contents, filename: str) -> bool:
    """The raw acquisition data, drawn end to end in acquisition order."""
    if not contents.acquisitions:
        return False
    fig = plt.figure(figsize=(11, 4))
    fig.suptitle(f'acquisitions: {filename}')
    for acq in contents.acquisitions:
        data = np.asarray(acq.data)
        trace = data[0] if data.ndim > 1 else data
        start = acq.head.acquisition_time_stamp_ns * 1.0e-9
        t = start + np.arange(trace.size) * acq.head.sample_time_ns * 1.0e-9
        plt.plot(t, np.real(trace), 'r', linewidth=0.4)
        plt.plot(t, np.imag(trace), 'g', linewidth=0.4)
    plt.xlabel('time (s)')
    plt.ylabel('signal')
    return True


# ---------- driver -------------------------------------------------------


def plot_mrd(input: BinaryIO, filename: str, *, raw: bool = False, save: Path = None) -> None:
    contents = read_contents(input)

    counts = ', '.join(f'{len(v)} {k or "untagged"}' for k, v in contents.arrays.items())
    print(f'found {len(contents.acquisitions)} acquisitions, {len(contents.images)} images'
          + (f', {counts}' if counts else ''), file=sys.stderr)

    drawn = False
    for figure in (plot_fitted_spectrum, plot_spectral_series, plot_metabolite_maps,
                   plot_time_courses, plot_diagnostics, plot_legacy_images):
        drawn |= figure(contents, filename)
    if raw:
        drawn |= plot_acquisitions(contents, filename)

    if not drawn:
        print('nothing in this file has a figure to draw; pass --raw to see the acquisitions',
              file=sys.stderr)
        return

    if save is not None:
        save.mkdir(parents=True, exist_ok=True)
        stem = filename or 'mrd'
        for num in plt.get_fignums():
            fig = plt.figure(num)
            title = (fig._suptitle.get_text().split(':')[0] if fig._suptitle else f'figure{num}')
            slug = re.sub(r'_+', '_', re.sub(r'[^0-9A-Za-z]+', '_', title)).strip('_')
            path = save / f'{stem}_{slug or f"figure{num}"}.png'
            fig.savefig(path, dpi=150, bbox_inches='tight')
            print(f'wrote {path}', file=sys.stderr)
        plt.close('all')
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot MRD file contents')
    parser.add_argument('-i', '--input', type=str, required=False, help='Input file, defaults to stdin')
    parser.add_argument('--raw', action='store_true', help='Also plot the raw acquisition data')
    parser.add_argument('-s', '--save', type=Path, default=None, help='Write the figures to this directory instead of showing them')
    args = parser.parse_args()

    if args.save is not None:
        matplotlib.use('Agg')

    if args.input is None:
        input = sys.stdin.buffer
        filename = ''
        plot_mrd(input, filename, raw=args.raw, save=args.save)
    else:
        filename = Path(args.input).stem
        with open(args.input, 'rb') as input:
            plot_mrd(input, filename, raw=args.raw, save=args.save)
