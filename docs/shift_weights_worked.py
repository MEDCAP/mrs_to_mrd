"""
Build docs/shift-methods-weights.pdf: how roll, alloc and contiguous fill one output sample.

Every method reads output sample j of switch i from the fractional source position
`i*total + j - offset_i` on the contiguous (zero-padded) readout, and differs only in the weights
it puts on the real samples around that position.  This draws those weights over real ischemia data,
for the real and the imaginary part, and asserts that the hand-computed sums equal what
mrd2shift.apply_roll / apply_alloc / apply_contiguous return.

    MPLBACKEND=Agg python docs/shift_weights_worked.py [-i file.mrd2] [-o out.pdf]
"""
import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mrd  # noqa: E402
import mrd2shift as ms  # noqa: E402

DEFAULT_INPUT = Path.home() / "dev/data/epsi_kidney_data/ischemia/ischemia_187_1_raw.mrd2"
HALF = 8
COLOURS = {"roll": "#d1495b", "alloc": "#edae49", "contiguous": "#00798c", "none": "#888888"}


def kernels(source: float) -> dict:
    """{method: (indices, weights)} for one fractional source position, restated from mrd2shift."""
    base = int(np.floor(source))
    frac = source - base
    taps = np.arange(base - HALF + 1, base + HALF + 1)
    x = source - taps
    return {"roll": (np.array([int(np.rint(source))]), np.array([1.0])),
            "alloc": (np.array([base, base + 1]), np.array([1.0 - frac, frac])),
            "contiguous": (taps, np.sinc(x) * np.sinc(x / HALF))}


def kernel_curve(method: str, x: np.ndarray, source: float) -> np.ndarray:
    """The continuous weight a sample at x would get, whose value at the integers is the stem."""
    d = x - source
    if method == "roll":
        return (np.abs(d) <= 0.5).astype(float)
    if method == "alloc":
        return np.clip(1.0 - np.abs(d), 0.0, None)
    return np.where(np.abs(d) < HALF, np.sinc(d) * np.sinc(d / HALF), 0.0)


def text_page(pdf, title, paragraphs, table=None):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.08, 0.94, title, fontsize=17, weight="bold", va="top")
    y = 0.89
    for para in paragraphs:
        mono = para.startswith("    ")
        lines = para.splitlines() if mono else textwrap.wrap(para, 92)
        for line in lines:
            fig.text(0.08, y, line, fontsize=10, va="top",
                     family="monospace" if mono else None)
            y -= 0.02
        y -= 0.012
    if table is not None:
        ax = fig.add_axes([0.08, 0.04, 0.84, max(0.05, y - 0.07)])
        ax.axis("off")
        t = ax.table(cellText=table[1], colLabels=table[0], loc="upper center", cellLoc="right")
        t.auto_set_font_size(False)
        t.set_fontsize(8)
        t.scale(1, 1.25)
    pdf.savefig(fig)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output", type=Path,
                        default=Path(__file__).resolve().parent / "shift-methods-weights.pdf")
    args = parser.parse_args()

    header, items = ms.read_stream(str(args.input))
    acqs = [it.value for it in items if isinstance(it, mrd.StreamItem.Acquisition)]
    groups = {}
    for acq in acqs:
        groups.setdefault(ms.switch_layout(header, acq), []).append(acq)
    layout = max(groups, key=lambda k: len(groups[k]))
    nswitch, total, kept, start = layout
    pool = [a for a in groups[layout] if not ms.is_prescan(a)]

    cube, _ = ms.acquisition_cube(pool, nswitch, total)
    raw_signal = np.abs(cube).sum(axis=(2, 3))
    drift = ms.measure_drift(cube, nswitch, total)
    slope = drift["slope"]
    pad = ms.measure_pad(ms.readout_anchor(drift["lines"], raw_signal, total), total, start, kept)
    offsets = ms.switch_offsets(nswitch, slope)
    used = nswitch * total
    target = start + kept // 2
    print(f"{args.input.name}: {nswitch}x{total}, window {start}..{start + kept - 1}, "
          f"slope {slope:+.4f}, pad {pad:+d}")

    # the brightest readout (centre of k-space, early repetition) is the one worth drawing
    lines = ms.pooled_lines(pool)
    brightest = int(np.argmax([np.abs(line).sum() for line in lines]))
    padded = ms.prepend_zeros(lines[brightest], pad).astype(np.complex128)
    samples = padded[0, :used]
    outputs = {"roll": ms.apply_roll(padded, offsets, nswitch, total)[0],
               "alloc": ms.apply_alloc(padded, offsets, nswitch, total)[0],
               "contiguous": ms.apply_contiguous(padded, offsets, nswitch, total)[0]}

    def source_of(i):
        return i * total + target - offsets[i]

    fracs = np.array([source_of(i) % 1 for i in range(nswitch)])
    # a middle switch whose fraction is near one half, where the kernels disagree most, and one
    # near a quarter; skip the first few, whose drift is still under a sample
    candidates = np.arange(8, nswitch - 8)
    half_switch = int(candidates[np.argmin(np.abs(fracs[candidates] - 0.5))])
    quarter_switch = int(candidates[np.argmin(np.abs(fracs[candidates] - 0.25))])
    chosen = [half_switch, quarter_switch]

    # the hand-computed sums are the methods themselves, or this document explains something else
    worked = {}
    for i in chosen:
        s = source_of(i)
        k = kernels(s)
        sums = {m: complex(np.sum(w * samples[idx])) for m, (idx, w) in k.items()}
        for m, v in sums.items():
            real = outputs[m][i * total + target]
            assert abs(v - real) < 1e-9 * max(1.0, abs(real)), (m, i, v, real)
        worked[i] = (s, k, sums)
    print(f"switches {chosen}: hand sums match apply_roll/apply_alloc/apply_contiguous")

    with PdfPages(args.output) as pdf:
        # ---- page 1: setup
        text_page(pdf, "How roll, alloc and contiguous fill one sample", [
            f"Data: {args.input.name}, brightest of {len(pool)} series readouts (the averaged "
            f"prescan is never shifted). {nswitch} switches of {total} samples; the reconstruction "
            f"reads positions {start}..{start + kept - 1} of each switch, with the echo expected at "
            f"{target} (ramp {start} + {kept}/2).",
            f"Measured drift: {slope:+.4f} samples per switch, so the gradient period is really "
            f"{total + slope:.4f} where the stream is cut every {total}. The echo walks "
            f"{slope * (nswitch - 1):+.2f} samples over the train. A {pad:+d} zero pad first puts "
            f"switch 0's echo on position {target}; the method then pulls every later switch back.",
            "Switch i is displaced by offset_i = -slope * i (anchored on switch 0). Output sample j "
            "of switch i is therefore read from the fractional position on the contiguous readout",
            "    source = i * total + j - offset_i = i * total + j + slope * i",
            "which in general falls between two measured samples. The three methods differ only in "
            "which real samples they read around `source` and with what weight:",
            "    roll        w = 1 at round(source)                       1 tap\n"
            "    alloc       w = 1-f at floor(source), f at floor+1       2 taps (linear)\n"
            "    contiguous  w(x) = sinc(x) * sinc(x/8),  x = source-n    16 taps (Lanczos-8)",
            "f is the fractional part of source. The weights are real, so the same weights are "
            "applied to the real and to the imaginary part of the complex signal; the plots show "
            "both. All three read along the contiguous readout, so a tap past the end of a switch "
            "reads the real neighbouring switch rather than wrapping into its own ramp.",
            f"Worked switches: {half_switch} (f = {fracs[half_switch]:.3f}, near one half, where the "
            f"kernels disagree most) and {quarter_switch} (f = {fracs[quarter_switch]:.3f}), both at "
            f"output position j = {target}. Every value in this document is checked against "
            "mrd2shift.apply_roll / apply_alloc / apply_contiguous when it is built.",
        ])

        # ---- page 2: weight curves on the discrete readout
        fig, axes = plt.subplots(3, 2, figsize=(11, 8.5), sharey=True)
        for col, i in enumerate(chosen):
            s, k, _ = worked[i]
            lo, hi = int(np.floor(s)) - HALF - 1, int(np.floor(s)) + HALF + 2
            xs = np.linspace(lo, hi, 2000)
            for row, m in enumerate(("roll", "alloc", "contiguous")):
                ax = axes[row, col]
                ax.axhline(0, color="k", lw=0.6)
                # switch boundaries and the kept window of the source switch
                for b in range(lo - lo % total, hi + 1, total):
                    ax.axvline(b - 0.5, color="k", ls=":", lw=0.8)
                w0 = (int(s) // total) * total
                ax.axvspan(w0 + start - 0.5, w0 + start + kept - 0.5, color="#dddddd", zorder=0)
                ax.plot(xs, kernel_curve(m, xs, s), color=COLOURS[m], lw=1.2, alpha=0.6)
                ints = np.arange(lo, hi + 1)
                ax.plot(ints, np.zeros_like(ints), "o", ms=3, color="#bbbbbb")
                idx, w = k[m]
                ml, sl, _ = ax.stem(idx, w, basefmt=" ")
                plt.setp(ml, color=COLOURS[m], ms=6)
                plt.setp(sl, color=COLOURS[m])
                for n, wn in zip(idx, w):
                    if abs(wn) > 0.02:
                        ax.annotate(f"{wn:+.3f}", (n, wn), textcoords="offset points",
                                    xytext=(0, 6 if wn >= 0 else -12), ha="center", fontsize=6.5)
                ax.axvline(s, color="k", lw=1.2)
                ax.set_xlim(lo, hi)
                ax.set_ylim(-0.3, 1.2)
                ax.set_ylabel(m)
                if row == 0:
                    ax.set_title(f"switch {i}: source = {s:.4f}  (f = {s % 1:.3f})", fontsize=10)
                if row == 2:
                    ax.set_xlabel("sample index on the contiguous readout")
        fig.suptitle("Weight each method puts on the real samples around the fractional source "
                     "(black line).\nCurve = continuous kernel, stems = the weights actually used; "
                     "grey band = the source switch's kept window, dotted = switch boundaries",
                     fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        pdf.savefig(fig)
        plt.close(fig)

        # ---- page 3: real and imaginary parts, each method's result at the fractional position
        for i in chosen:
            s, k, sums = worked[i]
            lo, hi = int(np.floor(s)) - HALF - 1, int(np.floor(s)) + HALF + 2
            ints = np.arange(lo, hi + 1)
            fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
            for ax, part, name in ((axes[0], np.real, "real"), (axes[1], np.imag, "imaginary")):
                ax.axhline(0, color="k", lw=0.5)
                ax.plot(ints, part(samples[ints]), "-", color="#bbbbbb", lw=1, zorder=1)
                # each measured sample, ringed by method, ring size by |weight|
                for m, marker, shift in (("contiguous", "o", 0), ("alloc", "s", 0), ("roll", "D", 0)):
                    idx, w = k[m]
                    ax.scatter(idx, part(samples[idx]), s=30 + 400 * np.abs(w), facecolors="none",
                               edgecolors=COLOURS[m], lw=1.5, marker=marker, label=f"{m} taps",
                               zorder=2)
                ax.plot(ints, part(samples[ints]), "o", color="k", ms=3, zorder=3)
                for m in ("roll", "alloc", "contiguous"):
                    ax.plot([s], [part(sums[m])], "*", ms=16, color=COLOURS[m], mec="k",
                            zorder=4, label=f"{m} result {part(sums[m]):+.4g}")
                ax.axvline(s, color="k", lw=1)
                for b in range(lo - lo % total, hi + 1, total):
                    ax.axvline(b - 0.5, color="k", ls=":", lw=0.8)
                ax.set_ylabel(f"{name} part")
                ax.legend(fontsize=7, ncol=2, loc="lower left", markerscale=0.7)
            axes[1].set_xlabel("sample index on the contiguous readout")
            fig.suptitle(f"Switch {i}, output position {target}: measured samples (black dots), the "
                         f"taps each method reads (rings, size = |weight|),\nand what each returns "
                         f"at source = {s:.4f} (stars). Same real weights on both parts.",
                         fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            pdf.savefig(fig)
            plt.close(fig)

        # ---- page: the table of taps for the half switch
        s, k, sums = worked[half_switch]
        cidx, cw = k["contiguous"]
        rows = []
        for n, wc in zip(cidx, cw):
            wr = 1.0 if n == k["roll"][0][0] else 0.0
            wa = dict(zip(k["alloc"][0], k["alloc"][1])).get(n, 0.0)
            v = samples[n]
            rows.append([str(n), f"{n % total}", f"{wr:.3f}", f"{wa:.3f}", f"{wc:+.4f}",
                         f"{v.real:+.4g}", f"{v.imag:+.4g}"])
        for m in ("roll", "alloc", "contiguous"):
            v = sums[m]
            rows.append([f"{m} sum", "", "", "", "", f"{v.real:+.4g}", f"{v.imag:+.4g}"])
        ang = {m: np.degrees(np.angle(sums[m])) for m in sums}
        text_page(pdf, f"Switch {half_switch}: every tap and the weighted sums", [
            f"source = {s:.4f}, f = {s % 1:.4f}. Columns: sample index on the contiguous readout, its "
            "position within its switch, the weight roll, alloc and contiguous put on it, and its "
            "real and imaginary value. Each sum is Σ w·(re + i·im) and equals the method's output.",
            "Magnitude / phase of each result:    "
            + "   ".join(f"{m} {abs(sums[m]):.4g} / {ang[m]:+.1f}°" for m in sums),
        ], table=(["index", "pos", "roll w", "alloc w", "contiguous w", "real", "imag"], rows))

        # ---- page: across the train
        frac_err = offsets - np.rint(offsets)
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        ax = axes[0, 0]
        ax.plot(-offsets, ".-", color="k", ms=3, label="true displacement  slope·i")
        ax.plot(-np.rint(offsets), drawstyle="steps-mid", color=COLOURS["roll"], label="roll: round")
        ax.set_xlabel("switch i")
        ax.set_ylabel("samples")
        ax.legend(fontsize=8)
        ax.set_title("displacement per switch")
        ax = axes[0, 1]
        ax.plot(frac_err, ".-", color=COLOURS["roll"], ms=3)
        for i in chosen:
            ax.axvline(i, color="k", ls=":", lw=0.8)
        ax.set_xlabel("switch i")
        ax.set_ylabel("samples")
        ax.set_title("roll error: true - rounded, a ±0.5 sawtooth")
        ax = axes[1, 0]
        fs = np.linspace(0, 1, 201)
        for f in (0.0, 0.25, 0.5):
            xx = np.linspace(-0.5, 0.5, 400)
            # response of each kernel to frequency x (cycles/sample) at this fraction
            alloc_r = np.abs((1 - f) + f * np.exp(-2j * np.pi * xx))
            ax.plot(xx, alloc_r, color=COLOURS["alloc"], alpha=0.4 + f, label=f"alloc f={f}")
        taps = np.arange(-HALF + 1, HALF + 1)
        for f in (0.25, 0.5):
            w = np.sinc(f - taps) * np.sinc((f - taps) / HALF)
            xx = np.linspace(-0.5, 0.5, 400)
            ax.plot(xx, np.abs(np.exp(-2j * np.pi * np.outer(xx, taps)) @ w),
                    color=COLOURS["contiguous"], alpha=0.4 + f, label=f"contiguous f={f}")
        ax.axhline(1, color=COLOURS["roll"], label="roll: gain 1, but position off by up to ±0.5")
        ax.set_xlabel("frequency along the readout (cycles/sample)")
        ax.set_ylabel("|gain|")
        ax.set_ylim(0, 1.15)
        ax.legend(fontsize=7)
        ax.set_title("what each kernel does to the readout's content")
        axes[1, 1].axis("off")
        maps = fig.add_gridspec(2, 8)[1, 4:].subgridspec(1, 4, wspace=0.08)
        for n, m in enumerate(("none", "roll", "alloc", "contiguous")):
            sig = ms.pooled_signal(lines, nswitch, total, lambda line, m=m: ms.apply_base_and_method(
                line, pad, offsets, slope, nswitch, total, m))
            ax = fig.add_subplot(maps[0, n])
            ax.imshow(sig, aspect="auto", origin="lower", cmap="magma",
                      extent=(-0.5, total - 0.5, -0.5, nswitch - 0.5))
            ax.axvspan(start - 0.5, start + kept - 0.5, fc="none", ec="w", ls="--", lw=0.8)
            ax.set_title(m if m != "none" else "pad only", fontsize=9, color=COLOURS[m])
            ax.set_xlabel("position", fontsize=8)
            ax.tick_params(labelsize=7)
            if n:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("switch i", fontsize=8)
        fig.suptitle(f"Across the train: {nswitch} switches, slope {slope:+.4f}", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)

        # ---- last page: takeaways
        text_page(pdf, "What the weights mean", [
            "roll reads one measured sample, the nearest one. It is exact only when f is 0; elsewhere "
            "it reads a point up to half a sample from where the echo really is. Its gain is 1 at "
            "every frequency, so it does not blur, but it leaves a position error that saws between "
            "-0.5 and +0.5 along the train. The switch index is the spectral axis, so that "
            "per-switch error lands in the spectrum as a periodic phase error.",
            "alloc is a two-tap triangle: 1-f on the sample below, f on the sample above. It lands "
            "on the right position on average, but a straight line between two samples is a "
            "low-pass filter. At f = 0.5 it simply averages two neighbours and passes nothing at "
            "the Nyquist frequency. The loss changes with f, so it too varies from switch to switch.",
            "contiguous is the band-limited answer: a sinc windowed to 16 taps (Lanczos-8). The weights "
            "alternate in sign and fall off with distance, the negative lobes being what restores "
            "the high-frequency content the triangle loses. Its gain is flat to near Nyquist at any f, "
            "so every switch is moved by the right fraction with the same response. Its taps near "
            "the ends of a switch reach into the real neighbouring switch; only at the two physical "
            "ends of the readout are they zero.",
            "All three are the same operation, Σ w_n · x_n at the same fractional source, with a "
            "different set of weights. That is the whole difference between them.",
        ])
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
