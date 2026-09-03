# Tests and sample data

## Running the tests

```
python3 -m unittest discover -s test -t .
```

numpy and scipy are the only requirements. stdlib `unittest` rather than pytest, because
this repo has no test dependencies and adding one to run a single suite is a poor trade;
pytest collects these fine if it is ever added.

`test_epsi_parity.py` compares the current EPSI reconstruction against the legacy
`mrd2_recon_to_incorporate.py`, stage by stage. It splits into two kinds of test:

- **identical** — the twelve stages the two are meant to agree on, asserted with no
  tolerance against `legacy_reference.py`. A failure means the current code has drifted
  from the reconstruction it inherited.
- **known delta** — the six stages where they genuinely differ. The current behaviour is
  the reference of record in every case, so these pin it and name the difference, keeping
  it a decision rather than a surprise.

`legacy_reference.py` is a transcription of the legacy arithmetic, not an import of it. The
legacy file cannot be run against anything this repo produces today: it imports `acqtypes`
and `lorn`, which are not on this branch; it reads `acq.data` as (samples, coils) where the
converter has written (coils, samples) since 68a5773; it reads the centre frequency from
`acquisition_center_frequency`, which nothing writes any more, so its ppm axis would divide
by zero; and it calls `plt.show()` throughout.

Nothing in the suite reads a scan from disk, and `conftest.py` supplies a stand-in for the
`mrd` package when it is not installed, so the suite runs anywhere. It prefers the real
package whenever it is importable, so it cannot mask a change in the mrd schema. To run
against the real one, use an interpreter that has it:

```
~/.local/share/mamba/envs/mrd/bin/python -m unittest discover -s test -t .
```

## Sample data

A sample folder tar of an epsi experiment, and a single file of fid spectral data.

- `ischemia_121.tar`
    a folder with typical epsi experiment, where single .MRD represent one repetition.
    The subfolders of `scan_id` represent a whole measurement of repetitions, which needs to be combined into a single file.
    There could be subfolders with different dimension of data, such as naverages>1 for phantom data. Convert this .MRD file as a single file.
- `cirrhrat_39_4.tar`
    A folder where single .MRD file represent the whole repetitions. The first few folders may have nrepetition=1 as they are test scan.
    It also carries a `cirrhrat_39_4_recon.mrd2`, a reference recon output.
