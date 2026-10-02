"""
Display the contents of an mrd2 file.

Three figures, each drawn when the file carries what it needs, so the same command works on a
converted stream, on a shift-corrected one and on a reconstruction:

    python mrdplot.py -i ischemia_179_epsigre_26575_recon.mrd2

    acquisitions                       the k-space itself, folded on the gradient switch with
                                       the window the fft keeps drawn on it, so a converter or a
                                       shift output plots with nothing else in the file. A recon
                                       carries its acquisitions through unchanged, so it shows
                                       the same samples the fit was run on. The averaged prescan
                                       is left out, and --switches folds a file whose header
                                       records no switch count
    metabolite_global_spect / _fit     the summed spectrum and the Lorentzian model fitted to it,
                                       sample by sample, with each peak's offset and the distance
                                       the fit had to move it
    metabolite_amplitude               one map per metabolite across the repetitions

With --save the figures are written as PNGs instead of shown, for running without a display, and
the transformed readouts are written beside them as a .mat.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import BinaryIO, Dict, List

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.io import savemat

import mrd

# the recon fits the series and the averaged prescan separately into one stream and names every
# array for its encoding. The series is what these figures are about
LABEL = "metabolite"


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
            elif isinstance(value, mrd.NdArray):
                contents.arrays.setdefault(describe(value), []).append(value)
    contents.acquisitions.sort(key=lambda a: a.head.acquisition_time_stamp_ns)
    return contents


def header_nswitches(header: mrd.Header) -> int:
    """The switch count the converter recorded, or 0 for a file with no EPSI readout."""
    if header is None or header.user_parameters is None:
        return 0
    for item in header.user_parameters.user_parameter_long:
        if item.name == "nswitches":
            return int(item.value)
    return 0


# ---------- figures ------------------------------------------------------


def plot_kspace(contents: Contents, filename: str, *, switches: int = 0) -> bool:
    """
    The k-space the acquisitions carry, as an image, for a file that holds any.

    Folded on the gradient switch when the switch count is known: every sample of an EPSI readout
    falls at one of `total` positions inside a switch, and the reconstruction keeps the same
    stretch of positions out of each of them, so the picture the window is judged by is position
    within the switch against switch. The echo should sit on one column in every row, inside the
    window and on the expected echo position, and an echo that walks across the columns is the
    drift mrd2shift takes out. That makes this the figure to run on a converter or a shift
    output, where there is nothing else in the file yet.

    A file whose header records no switch count cannot be folded, so its readouts are drawn end
    to end instead, one row per acquisition. Most conversions older than the current converter
    are in that state; --switches folds them anyway.

    The averaged prescan is left out. It calibrates a reconstruction rather than being
    reconstructed, and its own drift says nothing about the series beside it.
    Args:
        - switches: the switch count to fold on, overriding whatever the header records
    """
    if not contents.acquisitions:
        return False
    nswitch = switches or header_nswitches(contents.header)

    groups: Dict[int, List[mrd.Acquisition]] = {}
    for acq in contents.acquisitions:
        groups.setdefault(acq.head.encoding_space_ref or 0, []).append(acq)

    drawn = False
    for ref, acqs in sorted(groups.items()):
        if acqs[0].head.flags & mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT:
            print(f'leaving the {len(acqs)} prescan acquisitions of encoding {ref} out',
                  file=sys.stderr)
            continue
        samples = acqs[0].samples()
        # acquisitions at another geometry belong in no row of this image
        readouts = [np.abs(np.asarray(acq.data)[0]) for acq in acqs if acq.samples() == samples]
        total = samples // nswitch if nswitch > 1 else 0

        figure, axes = plt.subplots(figsize=(7, 9))
        if total >= 2:
            used = nswitch * total
            # summed over views and repetitions, which is the aggregate a single readout is too
            # noisy to show
            signal = np.sum([row[:used].reshape(nswitch, total) for row in readouts], axis=0)
            discard_pre = acqs[0].head.discard_pre or 0
            # what make_buffer keeps: whatever the two discards leave of the switch, which is the
            # readout axis the reconstruction transforms
            kept = total - discard_pre - (acqs[0].head.discard_post or 0)
            figure.suptitle(f'k-space per switch of encoding {ref}: {filename}')
            image = axes.imshow(signal, vmin=0, origin='lower', aspect='auto',
                                interpolation='nearest')
            axes.plot(np.argmax(signal, axis=1), np.arange(nswitch), 'xr', ms=4,
                      label='brightest sample')
            if kept > 0:
                # edges as well as a wash, so the window is legible without dimming the samples
                # inside it
                axes.axvspan(discard_pre - 0.5, discard_pre + kept - 0.5, color='w', alpha=0.12,
                             label=f'the {kept} points the fft reads')
                for edge in (discard_pre - 0.5, discard_pre + kept - 0.5):
                    axes.axvline(edge, color='w', lw=1.2, ls='--', alpha=0.8)
                # the echo belongs at the middle of the kept window, where k-space crosses zero
                echo = discard_pre + kept // 2
                axes.axvline(echo, color='r', lw=1.5, alpha=0.6,
                             label=f'expected echo at {echo}')
            axes.set_xlabel(f'position within the {total} point switch')
            axes.set_ylabel('switch')
            axes.set_title(f'{len(readouts)} of {len(acqs)} acquisitions summed, '
                           f'{nswitch} switches of {total} points, '
                           f'discard_pre {discard_pre}', fontsize=8)
            axes.legend(fontsize=8, loc='upper right')
            figure.colorbar(image, ax=axes, label='signal summed over views and repetitions')
        else:
            figure.suptitle(f'readouts of encoding {ref}: {filename}')
            rows = np.stack(readouts)
            # a handful of samples at the head of the first repetitions are hundreds of times the
            # rest, and scaling to them leaves every other row flat
            image = axes.imshow(rows, vmin=0, vmax=np.percentile(rows, 99.5), origin='lower',
                                aspect='auto', interpolation='nearest')
            axes.set_xlabel(f'sample of the {samples} point readout')
            axes.set_ylabel('acquisition, in the order the file holds them')
            axes.set_title(f'{len(readouts)} of {len(acqs)} acquisitions, no switch count '
                           f'recorded: pass --switches to fold them on the switch', fontsize=8)
            figure.colorbar(image, ax=axes, label='signal')
        figure.tight_layout()
        drawn = True

    return drawn


def placed_centers(spect: np.ndarray, xscale: np.ndarray, offsets: np.ndarray,
                   biggest: int) -> np.ndarray:
    """
    Where the fit put each peak before it was allowed to move, in ppm.

    The same arithmetic candidate_centers ran: the offsets are a rigid pattern of known
    chemistry, so the fit starts by anchoring the peak it believes is the tallest on the tallest
    point of the spectrum and hanging the others off it. Wrapped modulo the spectral width,
    because the axis is relative and a peak past the end folds round rather than falling off.
    Args:
        - spect: the summed spectrum the fit was run on
        - xscale: its ppm axis, evenly spaced
        - offsets: the peak offsets in ppm, in peak order
        - biggest: which peak the winning hypothesis anchored on
    Returns:
        - one placed center per peak, in ppm
    """
    bw_ppm = float(xscale[-1] - xscale[0] + (xscale[1] - xscale[0]))
    anchor = float(xscale[int(np.argmax(np.abs(spect)))])
    return (anchor - (np.asarray(offsets, dtype=float) - offsets[biggest])) % bw_ppm


def plot_lorentzian_fit(contents: Contents, filename: str) -> bool:
    """
    The summed spectrum and the model fitted to it, sample by sample.

    Drawn as the samples that were fitted rather than as a curve through them: the spectrum is
    one point per switch, so a smooth line would be drawing resolution the data does not have,
    and a peak two samples wide is exactly the case worth seeing honestly.

    Each peak is labelled with the offset it was named by and the distance the fit moved it from
    where that offset placed it. That delta is what says whether the pattern landed: a peak that
    had to walk a long way to fit is a peak fitted to somebody else's line, whatever the residual
    says about the model as a whole.
    """
    spectra = contents.arrays.get(f"{LABEL}_global_spect")
    if not spectra:
        return False
    if len(spectra) > 1:
        print(f'{len(spectra)} {LABEL}_global_spect arrays; drawing the first', file=sys.stderr)

    spect = spectra[0]
    data = np.asarray(spect.data)
    xscale = np.array(meta_values(spect, "xscale_ppm"))
    if xscale.size != data.size:
        xscale = np.arange(data.size)

    figure = plt.figure(figsize=(11, 5))
    figure.suptitle(f'fitted spectrum: {filename}')
    plt.plot(xscale, np.real(data), 'r.-', ms=4, lw=0.8, label='real')
    plt.plot(xscale, np.imag(data), 'g.-', ms=4, lw=0.8, label='imag')

    fit = contents.first(f"{LABEL}_global_spect_fit")
    if fit is not None:
        model = np.asarray(fit.data)
        plt.plot(xscale, np.real(model), 'k.--', ms=3, lw=0.8, label='fit real')
        plt.plot(xscale, np.imag(model), '.--', ms=3, lw=0.8, color='0.5', label='fit imag')

    centers = contents.first(f"{LABEL}_lorentzian_centers_ppm")
    if centers is not None:
        fitted = np.asarray(centers.data).real
        names = meta_values(centers, "peak_names")
        # the offsets the peaks were named by, which only the maps carry
        maps = contents.first(f"{LABEL}_amplitude") or contents.first(f"{LABEL}_area")
        offsets = np.array(meta_values(maps, "peak_offsets_ppm")) if maps is not None else None
        biggest = meta_value(spect, "biggest_peak_index")
        deltas = None
        if offsets is not None and offsets.size == fitted.size and biggest is not None:
            placed = placed_centers(data, xscale, offsets, int(biggest))
            bw_ppm = float(xscale[-1] - xscale[0] + (xscale[1] - xscale[0]))
            # the short way round the axis, since a center that folded past the end moved a
            # sample rather than a whole spectral width
            deltas = (fitted - placed + bw_ppm / 2) % bw_ppm - bw_ppm / 2

        blended = plt.gca().get_xaxis_transform()
        # peaks can sit close together, so the labels run vertically and alternate height
        for ip, center in enumerate(fitted):
            plt.axvline(center, color='0.4', linewidth=0.8)
            label = names[ip] if ip < len(names) else str(ip)
            if deltas is not None:
                label += f' {deltas[ip]:+.3f}'
            label += f' @{center:.2f}'
            # placed in axis fractions and clipped, so a long label stays inside the axes
            # rather than climbing over the title
            plt.text(center, 0.22 + 0.34 * (ip % 2), label, transform=blended,
                     rotation=90, va='bottom', ha='right', fontsize=7, clip_on=True,
                     bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.5))

    anchor = meta_value(spect, "biggest_peak_name")
    subtitle = f'anchored on {anchor}' if anchor else ''
    loss = meta_value(fit, "fit_loss") if fit is not None else None
    if loss is not None:
        subtitle += f', residual {loss:.3f}'
    plt.title(f'{subtitle}. one marker per sample, d is how far the fit moved each peak',
              fontsize=8)
    plt.xlabel('frequency (ppm)')
    plt.ylabel('amplitude')
    # outside the axes, so it cannot sit on top of a peak
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    plt.tight_layout()
    return True


def montage(maps: np.ndarray) -> np.ndarray:
    """
    Lay (npeaks, nreps, ny, nx) out as one image, peaks down and repetitions across.

    Each peak is scaled by its own maximum, so a weak metabolite is still visible next to the
    substrate, and a bright line separates the rows. The voxels are laid out at the size they
    were fitted at: a map is a dozen voxels across, and interpolating it up would draw structure
    the reconstruction never produced.
    """
    npeaks, nreps, ny, nx = maps.shape
    out = np.zeros((npeaks * ny, nreps * nx))
    for ipeak in range(npeaks):
        peak_max = np.max(np.abs(maps[ipeak]))
        if peak_max == 0:
            continue
        for irep in range(nreps):
            out[ipeak * ny:(ipeak + 1) * ny, irep * nx:(irep + 1) * nx] = (
                maps[ipeak, irep] / peak_max)
        out[ipeak * ny, :] = 1
    return out


def plot_metabolite_maps(contents: Contents, filename: str) -> bool:
    """Each metabolite's map across the repetitions, one row per peak."""
    arr = contents.first(f"{LABEL}_amplitude") or contents.first(f"{LABEL}_area")
    if arr is None:
        return False
    maps = np.asarray(arr.data).real
    if maps.ndim != 4:
        return False

    names = meta_values(arr, "peak_names")
    npeaks, nreps, ny, nx = maps.shape
    tiles = montage(maps)

    # voxels stay square, so the figure is sized to the montage rather than the other way round,
    # which is what keeps a wide, short montage from floating in empty space
    width_in = float(np.clip(nreps * 0.7, 6, 16))
    figure = plt.figure(figsize=(width_in, width_in * tiles.shape[0] / tiles.shape[1] + 1.2))
    figure.suptitle(f'{describe(arr)}: {filename}')
    plt.imshow(tiles, cmap='gray', interpolation='nearest')
    plt.xlabel('repetition')
    plt.xticks(np.arange(nreps) * nx + nx / 2, [str(i) for i in range(nreps)], fontsize=6)
    plt.yticks(np.arange(npeaks) * ny + ny / 2,
               [names[i] if i < len(names) else str(i) for i in range(npeaks)])
    plt.title('each row scaled to its own maximum', fontsize=8)
    plt.tight_layout()
    return True


def kspace_cube(acqs: List[mrd.Acquisition], nswitch: int) -> np.ndarray:
    """
    One encoding's acquisitions as the complex cube the reconstruction fills.

    Indexed (repetition, view, readout point, switch), which is the arrangement every stage
    downstream reads: the readout axis is the points the fft keeps out of each switch, and the
    switch axis is the spectral one, since one spectral point is acquired per switch. Views and
    repetitions come off the indices the converter set rather than off their counts, so a group
    missing a view still lands in the right row.

    The samples are the ones the file holds, with no line broadening: that is a reconstruction
    parameter rather than a property of the data.
    Args:
        - acqs: the acquisitions of one encoding, all at the same geometry
        - nswitch: the switch count the readout was acquired at
    Returns:
        - the cube, or None when the discards leave no readout to keep
    """
    samples = acqs[0].samples()
    total = samples // nswitch
    discard_pre = acqs[0].head.discard_pre or 0
    kept = total - discard_pre - (acqs[0].head.discard_post or 0)
    if total < 1 or kept <= 0:
        return None

    views = sorted({acq.head.idx.kspace_encode_step_1 or 0 for acq in acqs})
    reps = sorted({acq.head.idx.repetition or 0 for acq in acqs})
    view_at = {view: i for i, view in enumerate(views)}
    rep_at = {rep: i for i, rep in enumerate(reps)}

    used = nswitch * total
    cube = np.zeros((len(reps), len(views), kept, nswitch), dtype=complex)
    for acq in acqs:
        if acq.samples() != samples:
            continue
        # reshaped through the sample axis, which comes first, so this is the switch-major order
        # the samples were acquired in, then cut to the window the fft reads
        body = np.asarray(acq.data)[0, :used].reshape(nswitch, total)
        cube[rep_at[acq.head.idx.repetition or 0],
             view_at[acq.head.idx.kspace_encode_step_1 or 0]] = (
                 body[:, discard_pre:discard_pre + kept].T)
    return cube


def save_mat(contents: Contents, stem: str, save: Path, *, switches: int = 0) -> bool:
    """
    Write the transformed data beside the figures as a .mat, for reading it somewhere else.

    One array per encoding and nothing else: `rawdata` for the series and `phantom` for the
    averaged prescan, each complex and indexed (repetition, view, readout point, switch). The
    prescan is written here though the figures skip it, since it is data somebody may want; what
    it is not is evidence about where the series' echo sits.

    Transformed the way mrd2recon transforms it, one repetition at a time over all three of its
    axes at once, because a DFT along the view axis is a sum over every view and no single
    acquisition holds more than one of them. The shape is unchanged by that, so the axes stay as
    named above and read as (repetition, y, x, frequency) afterwards.

    No line broadening is applied. That is a reconstruction parameter rather than a property of
    the samples, and applying one here would bake a choice the recon still has to make into the
    export.
    Returns:
        - True when a file was written
    """
    nswitch = switches or header_nswitches(contents.header)
    if nswitch <= 1:
        return False

    groups: Dict[int, List[mrd.Acquisition]] = {}
    for acq in contents.acquisitions:
        groups.setdefault(acq.head.encoding_space_ref or 0, []).append(acq)

    payload: Dict[str, np.ndarray] = {}
    for ref, acqs in sorted(groups.items()):
        cube = kspace_cube(acqs, nswitch)
        if cube is None:
            continue
        for rep in range(cube.shape[0]):
            axes = (0, 1, 2)
            cube[rep] = np.fft.fftshift(np.fft.fftn(cube[rep], axes=axes), axes=axes)
        key = 'phantom' if acqs[0].head.flags & mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT \
              else 'rawdata'
        # a file with two encodings of the same kind would otherwise write one over the other
        payload[key if key not in payload else f'{key}_{ref}'] = cube

    if not payload:
        return False
    path = save / f'{stem}.mat'
    savemat(str(path), payload)
    print(f'wrote {path}: {", ".join(payload)}', file=sys.stderr)
    return True


# ---------- driver -------------------------------------------------------


def plot_mrd(input: BinaryIO, filename: str, *, switches: int = 0,
             save: Path = None) -> None:
    contents = read_contents(input)

    drawn = plot_kspace(contents, filename, switches=switches)
    drawn |= plot_lorentzian_fit(contents, filename)
    drawn |= plot_metabolite_maps(contents, filename)

    if not drawn:
        print('nothing in this file has a figure to draw', file=sys.stderr)
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
        save_mat(contents, stem, save, switches=switches)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot MRD file contents')
    parser.add_argument('-i', '--input', type=str, required=False, help='Input file, defaults to stdin')
    parser.add_argument('--switches', type=int, default=0, help='Fold the readouts on this many gradient switches, for a file whose header records no switch count')
    parser.add_argument('-s', '--save', type=Path, default=None, help='Write the figures to this directory instead of showing them')
    args = parser.parse_args()

    if args.save is not None:
        matplotlib.use('Agg')

    if args.input is None:
        input = sys.stdin.buffer
        filename = ''
        plot_mrd(input, filename, switches=args.switches, save=args.save)
    else:
        filename = Path(args.input).stem
        with open(args.input, 'rb') as input:
            plot_mrd(input, filename, switches=args.switches, save=args.save)
