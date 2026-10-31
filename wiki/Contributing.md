# Contributing

Bug reports, questions, ideas and pull requests are all welcome. The full guide is in
[CONTRIBUTING.md](https://github.com/duceppemo/BrucellaID/blob/main/CONTRIBUTING.md); this page is the short
version, plus what the code looks like from the inside.

## Reporting something
[Open an issue](https://github.com/duceppemo/BrucellaID/issues/new/choose) with the version
(`brucellaid --version`), the command you ran and the output of a run with `--verbose`. For a species or a
sequence type that looks wrong, add the row of `brucellaid.tsv` and the rows of `brucellaid_mlst.tsv` and
`brucellaid_mash.tsv` for that sample: they show what the call rests on.
[Troubleshooting](Troubleshooting) covers the usual causes.

## Development setup
```bash
git clone https://github.com/duceppemo/BrucellaID
cd BrucellaID
conda env create -f environment.yml
conda activate brucellaid
pip install --no-deps -e .
pytest --cov
```
The suite covers 100% of the package; the runs that need MASH or Seal are skipped when the tool is missing.
Lint with `ruff check .`, or install the hook with `pre-commit install`. `ruff format` is not used.

## The code
| Module | What it does |
| --- | --- |
| `cli.py` | Arguments, logging, and the order of the run |
| `samples.py` | Finds the fastq files of a folder and pairs R1 with R2 |
| `mash.py` | Pipes the reads into `mash dist` and reads the distances back |
| `references.py` | The reference table, and what a MASH sketch ID stands for |
| `seal.py` | Runs `seal.sh`, including the paths it cannot be given as they are |
| `mlst.py` | Turns Seal's per-allele counts into allele calls and a sequence type |
| `pipeline.py` | One sample, then all of them, with the errors kept per sample |
| `report.py` | The four reports, and the legacy xlsx one |

Two rules the code follows throughout: anything the user can fix is a `BrucellaIDError`, reported without a
traceback, and a sample that fails keeps its error in its own result instead of ending the run.

## The reference data
`brucellaid/data/` holds the MASH sketch, the reference table, the MLST alleles and the MLST profiles.
`scripts/make_mlst_data.py` rebuilds the MLST files from PubMLST and `scripts/make_reference_data.py`
rebuilds the sketch from the accessions of the reference table. Both print what changed; say in the pull
request which day the data was downloaded, and add a line to the [Changelog](Changelog). A sketch must keep
MASH's defaults (k = 21, 1000 hashes), because the threshold of 150 shared hashes is read against that
sketch size.

## Documentation
This wiki lives in the [`wiki/`](https://github.com/duceppemo/BrucellaID/tree/main/wiki) folder of the
repository and is published from there when changes reach `main`. **Do not edit it on GitHub directly**:
those changes are overwritten by the next sync.
