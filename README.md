# todo
- verify recon and lorn fit from main is inherited into add-ndarray branch
- verify other methods of 

# EPSI echo drift correction

`mrd2shift.py` is the middle stage of the pipeline:

    tar -> MRStomrd2.py -> mrd2shift.py -> mrd2recon.py -> mrdplot.py

It takes the echo drift out of a converted MRD v2 stream. This note explains what the drift is, what
gets measured, and what each `--method` does to the readout. Every number below is measured on
`cirrhrat_43_1`, 64 switches of 28 samples, reading 4..15.

## The defect

The scanner records `no_samples` and `no_switches` separately, so everything downstream divides:
1792 / 64 = 28 samples per switch. That division is wrong. The real gradient period is 28.1937
samples, so switch *i* of the stream starts 0.1937 × *i* samples out of step with the gradient cycle
it is meant to describe. Across a 64 switch train the echo walks +12.2 samples.

The reconstruction reads the same window out of every switch, so early switches are read correctly
and late ones collect gradient ramp instead of signal.

Two axes, and which is which matters. Position within a switch is **kx**, the readout direction. The
switch index is **time**, and becomes the spectral axis after the transform.

## What is measured

Nothing here is supplied by hand. Two numbers come off the data.

**The slope.** `fit_peak_lines` fits the per-switch peaks as two parallel lines sharing one slope.
Two lines because a switch crosses k-space centre twice — the readout echo on the gradient plateau,
and the rephasing echo after it — so the brightest position in a switch is whichever of them won
there, and what looks like scatter is two orderly families. Sharing the slope means every switch
constrains it, not only the ones that happened to peak on the readout echo. On cirrhrat this reads
**+0.1937 per switch**, a period of 28.1937, rms 0.62.

The separation between the two lines is pinned at `total/2` rather than fitted. The two crossings sit
`ppsw/2 + 2·ramp` apart and `period = ppsw + 4·ramp`, so that spacing is identically `period/2` and
carries no information. Worse, `period/2` is exactly where a two-line model cannot resolve: a line
plus `period/2` and a line minus `period/2` are the same line relabelled, so least squares is
*repelled* from the truth. Measured on this scan the loss reads 22.2 at separation 12, 81.3 at the
true 14, and 22.2 again at 16.

**The pad.** `readout_anchor` reads where switch 0's readout echo sits, and `measure_pad` returns the
gap between that and where the sequence puts it, `discard_pre + npoints_per_switch/2` = 4 + 6 = 10.
On cirrhrat the anchor is 4.05, so the pad is **+6**.

Picking the right echo matters. The two crossings are exactly `total/2` apart, so both candidate pads
tie in magnitude with opposite signs and proximity cannot separate them. Brightness can: the readout
echo is on the plateau where the signal is, so the family whose peaks average brighter is the readout
family. Reading it off the fitted line rather than off switch 0 itself matters just as much — on this
scan switch 0's own argmax lands on the *rephasing* echo at 17, and anchoring there would have put
the readout echo at position 1 of 28.

The measured +6 sits one sample off the `ramp + 3` = +7 an older converter hard-coded. That is a real
disagreement between the hardcode and the data, not a bug.

## The base: zero-padding

Before any drift method runs, `prepend_zeros` pushes `pad` zeros into the front of the whole readout
and drops as many off the tail. This is a prepend over the contiguous readout, not a roll: what falls
off the end of the train is gone rather than wrapped back to its start, and the positions switch 0
gains at its front were never measured — zero says that, where a wrapped sample would claim the tail
of the train was acquired before its head.

Switch 0 carries the least drift and the most signal, so once it is on target the first few switches
can be trusted, and everything left to correct is the accumulation along the train. That is what the
methods below deal with.

**The pad and the drift correction must share an anchor, and the anchor is switch 0.** `switch_offsets`
returns `delta_i = -slope · i`, so switch 0 does not move and every later switch is pulled back onto
it. An anchor on the middle of the train would hold the *mean* echo position instead, and the two
would compose into a readout sitting half the drift — six samples of a twelve point window here —
past where the pad aimed it.

The cost is that the largest displacement now falls at the end of the train, `slope · (M-1)` = 12.2
samples rather than half that. The methods that wrap inside a switch pay it; the ones that read from
the neighbouring switch do not.

## What is not corrected

**The prescan is never shifted.** It is calibration rather than data — a separate acquisition that
happens to travel in the same file — and the scaling the reconstruction takes from it has to mean the
same thing before and after a correction. Left alone, cirrhrat's phantom scaling reads 625.678
whatever `--method` ran. Corrected along with the series it read 599.589 under `phase`, 494.786 under
`regrid` and 284.540 under `alloc`: a calibration drifting with the thing it is supposed to calibrate,
and under `alloc` less than half its true value.

The prescan is also the case that shows why the series' numbers cannot simply be scaled onto another
layout. It divides into 64 switches of 34 samples where the series has 28, and its drift cannot be
measured: fitted on its own data it gives +0.0937 per switch at rms 4.81, well past the 1.5 that
separates a readout with an echo from one where the argmax is landing on noise. Scaling the series'
slope by switch length would give +0.2352, which fits the prescan's own peaks *worse* than its own
estimate does (two-line cost 1171 against 1026, and both are bad). Neither number is trustworthy, so
neither is applied.

It is also what turned up a bug in `switch_layout`. `npoints_per_switch` is written once for a whole
file, so reading `kept` off the header gave the prescan the series' 12, and `4·4 + 12 = 28` does not
account for a 34 point switch — which put the prescan's window and its expected echo position six
samples from where they are. Taken from each acquisition's own discards instead, `34 − 4 − 12 = 18`,
and `4·4 + 18 = 34` closes exactly. Both of this scan's layouts are well formed; the header was the
only thing that was not. `geometry_closes` stays as the guard for layouts that genuinely do not close
— the kidney data, where `4×2.5 + 12` is 22 against a 20 point switch, is one — and any layout failing
it is left alone.

Any non-prescan layout is measured on its own data rather than by extrapolation, and left alone when
that measurement is unusable. Being left alone is the better failure.

## The methods

### `regrid` — resample onto the period actually acquired

The other methods displace each switch. This one re-grids the whole readout once.

The defect is that the stream is divided by `total` where the gradient period is `total + slope`, so
switch *i* of the real sequence begins at sample `i·(total + slope)`, not `i·total`. `regrid` reads
there and writes to `i·total`:

    source[i, p] = i · (total + slope) + p        for p in 0 .. total-1
    output[i, p] = resample(readout, source[i, p])

Every switch lands back on the grid the reconstruction assumes. Because the source positions are
fractional, the read is done by `resample`: a 16-tap Lanczos-windowed sinc (`half_width = 8`), which
is band-limited interpolation rather than rounding. A readout is one continuous time series, so a sample
between two of its points is genuinely recoverable that way; the window keeps truncating the kernel
from ringing, and anything outside the readout reads as zero, since at the two physical ends there is
no neighbour and zero says "not measured" where clamping would repeat an edge.

Anchored on switch 0 like everything else: reading straight from `i·(total + slope)` leaves every echo
at the within-period offset switch 0 already had, which the pad has already put on target. `lead`
would move that grid and is left at zero.

What it buys: no per-switch discontinuity at the boundaries, and nothing wrapped. What it costs: the
source positions run 0.0 to 1803.2 while the readout holds 1792 samples, so **13 samples at the tail
of the last switch read as zero**. That is the honest price of reading a longer period out of a fixed
number of samples.

### `contiguous` — displace along the contiguous readout

`shift` is cyclic within a switch, so the samples it pulls in at the ends of the train — where the
displacement is largest — are that switch's own ramp and rephasing points. Those are junk, and
burying the weak metabolites under them is what costs the wrapping methods their hydrate. The readout
is one continuous time series, so the samples that really sit beside the window are the neighbouring
switch's. `contiguous` takes those, through the same `resample`:

    source[i, p] = (i · total + p) - offsets[i]     where offsets[i] = -slope · i

**With the switch-0 anchor this is the same operation as `regrid`.** Substituting the offsets gives
`i·total + p + slope·i = i·(total + slope) + p`, which is `regrid`'s source exactly; measured, the two
agree to 2×10⁻¹³ and produce identical reconstructions down to the phantom scaling (494.786 both).
They were distinct under the old mid-train anchor, which offset `regrid` by a constant `middle`. Both
are kept because they are written from different starting points — one from the period, one from the
per-switch displacement — and either could be the one to generalise later.

### `shift` — exact, cyclic within the switch

The displacement applied over the whole switch by the Fourier shift theorem, so it is exact for a
fractional number of samples. Because it is cyclic, what leaves one end of a switch arrives at the
other. That is fine while the displacement is small and increasingly wrong towards the end of the
train, where 12.2 samples of a 28 point switch wrap the switch's own ramp into the window.

Moving the whole switch is what re-centres the k-space line on the echo and fixes the asymmetric
truncation that broadens the spatial point spread — it sharpens the metabolite map, which a phase
confined to the window cannot.

### `roll` — `shift` with the fractional part rounded away (control)

`np.roll` by `round(offset)`. This is the hardcoded stair-step correction (`shift=0` for `i<3`, `1`
for `i<6`, …) generated from the measured slope instead of written out by hand. It is a control: the
difference between it and `shift` is what the fractional part is worth.

### `alloc` — `shift` with the fractional part split between neighbours (control)

Linear interpolation between the two neighbouring integer rolls. The two-tap split is a convolution
along kx and therefore a multiplication across the image, so it does not blur — it shades the readout
edges of the field of view to |1−2f| of the centre.

### `phase` — dropped

There was a sixth method: a phase ramp applied to just the samples the reconstruction reads, rotating
the k-space line within itself so no sample enters from outside the window. It was meant to repair the
phase relationship along the switch train while leaving the metabolite map alone.

It is gone, because leaving the readout alone is the whole problem. It cannot straighten the echo by
construction — measured, it left the walk at sd 3.97 against the raw 3.92 — and through a full
reconstruction it made the global fit *worse* than no correction at all, residual 1.648 against 1.284.
A correction that cannot move the echo has nothing to contribute here.

## Results on cirrhrat_43_1

Readout echo position per switch, windowed to 10 ± 6 so the rephasing echo cannot win the argmax:

| | sw0 | sw1 | sw15 | sw30 | sw45 | sw63 | sd |
|---|---|---|---|---|---|---|---|
| raw | 16 | 4 | 7 | 10 | 13 | 16 | 3.92 |
| pad only | 9 | 9 | 13 | 16 | 4 | 5 | 4.80 |
| `shift` | 9 | 9 | 10 | 11 | 10 | 16 | 1.57 |
| `regrid` | 9 | 9 | 10 | 10 | 10 | 10 | **1.29** |
| `contiguous` | 9 | 9 | 10 | 10 | 10 | 10 | **1.29** |
| `roll` | 9 | 9 | 10 | 10 | 10 | 16 | 1.66 |
| `alloc` | 9 | 9 | 10 | 10 | 10 | 16 | 1.61 |

The pad puts switch 0 on target and leaves the walk untouched, which is what it is for. `regrid` and
`contiguous` then hold every switch on 9–10; the wrapping methods hold it until the far end and lose
switch 63.

Through a full reconstruction — global fit residual, and each metabolite's area summed over the
signal-bearing repetitions 1–8:

| method | residual | bic | urea | pyr | ala | hyd | lac |
|---|---|---|---|---|---|---|---|
| before | 1.284 | 1871 | 8233 | 8500 | 728 | 1485 | 746 |
| `shift` | 1.025 | 844 | 5595 | 7950 | 722 | 2031 | 792 |
| `regrid` | **1.024** | 1237 | 5713 | 7880 | 844 | 2465 | 783 |
| `contiguous` | **1.024** | 1237 | 5713 | 7880 | 844 | 2465 | 783 |
| `roll` | 2.660 | 930 | 5635 | 7953 | 1403 | 743 | 1165 |
| `alloc` | 1.338 | 1021 | 5184 | 7300 | 818 | 2255 | 648 |

`regrid` and `contiguous` give the lowest residual and the largest hydrate. `roll` is the worst of
all, worse than no correction: rounding the displacement to whole samples costs more than the drift
it removes, which is the argument for keeping it only as a control.

**One caveat on reading that table.** A lower global-fit residual does not by itself prove a better
reconstruction — it can also mean a small peak has slid onto a strong neighbour. Here the hydrate
centre moves from 3.40 ppm before correction to 2.63 ppm after `regrid`, 0.77 ppm towards urea, while
urea's own area falls to 0.69× and hydrate's rises to 1.66×. Some of that hydrate gain is real
recovery and some may be borrowed from urea. Check the peak placement in the fitted-spectrum figures
before quoting the hydrate number.

## Figures

`docs/` holds every figure referenced here, regenerated by the commands in the next section.

- `echo_position_all_methods.png` — position within a switch across, switch number up, for the raw
  readout and each of the five methods over the zero-pad base. A drifting echo is a slanted stripe and
  a corrected one is vertical. The reconstruction window is shaded and the expected echo position is a
  dotted line. The base is not drawn on its own: it translates the whole picture by `pad` and leaves
  the walk exactly as it was, which the raw panel already shows.
- `<method>_recon_fitted_spectrum.png` — the summed spectrum, the Lorentzian model fitted to it, and
  where the peaks landed. The figure the fit is judged by: if the model does not sit on the data, or
  a peak centre is not on a peak, the maps are not worth reading.
- `<method>_recon_metabolite_area.png`, `<method>_recon_metabolite_amplitude.png` — the maps as a
  peak × repetition montage, each row scaled to its own maximum so a weak metabolite stays visible
  next to the substrate.
- `<method>_recon_metabolite_time_courses.png` — each metabolite summed over voxels against
  repetition, substrate solid and metabolites dashed.
- `<method>_recon_phantom_peak_area_scaling_625_678.png` — the prescan map. The scaling is in the
  filename and reads 625.678 for every method, because the prescan is not shifted.

`<method>` is `before`, `shift`, `regrid`, `contiguous`, `roll` or `alloc`.

## Usage

    # correct a stream, measuring both numbers off the data (--method regrid is the default)
    python mrd2shift.py -i raw.mrd2 -o straight.mrd2

    # as a pipeline stage; both default to $INPUT_PIPE / $OUTPUT_PIPE
    cat raw.mrd2 | python mrd2shift.py -i - -o - | ...

    # dry run: draw every method and write no stream
    python mrd2shift.py -i raw.mrd2 -p docs/echo_position_all_methods.png

    # override the measured pad, or turn the base off with --pad 0
    python mrd2shift.py -i raw.mrd2 -o straight.mrd2 --pad 7

Options beyond `--method`, `--pad` and `-p/--plot`: `--drop-first-switch` zeroes the first switch,
whose sample is the integral of the spectrum and fits a sum of Lorentzians badly; `--force` corrects a
stream that already records a drift taken out of it. The slope is always measured and has no flag.

What was applied is recorded on the header as `echo_drift_slope_applied`, `echo_drift_pad_applied` and
`echo_drift_method`, so a second pass refuses rather than correcting twice.

Every diagnostic goes to stderr — anything on stdout lands in the middle of the stream and the next
stage dies on its magic bytes.

To regenerate everything in `docs/` from a converted stream:

    for m in shift regrid contiguous roll alloc; do
        python mrd2shift.py -i before.mrd2 -o $m.mrd2 --method $m
    done
    for m in before shift regrid contiguous roll alloc; do
        python mrd2recon.py -i $m.mrd2 -o ${m}_recon.mrd2 \
            -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8
        python mrdplot.py -i ${m}_recon.mrd2 -s docs/
    done
    python mrd2shift.py -i before.mrd2 -p docs/echo_position_all_methods.png

---

## Notes on reading MRS data

For reading mrs data everywhere
- copy MRSreader.py as a class definition (generic)
- MRStomrd2.py but does not have to be generic since all sequences will be in .MRD2 format
- If necessary, make a slim version of mrs converter that does not depend on organize
- MRStomrd2.py can be hyper locally specific to this case
