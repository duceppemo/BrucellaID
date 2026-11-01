# Contributing to BrucellaID

Thanks for your interest! Bug reports, questions, ideas and pull requests are all welcome.
Please follow the [code of conduct](CODE_OF_CONDUCT.md).

## Reporting a bug or asking a question
[Open an issue](https://github.com/duceppemo/BrucellaID/issues/new/choose) and pick the matching template.
For bugs, please include the version (`brucellaid --version`), the command you ran and the messages it printed
(ideally from a run with `--verbose`). For a species or a sequence type that looks wrong, include the row of
`brucellaid.tsv` and the rows of `brucellaid_mlst.tsv` and `brucellaid_mash.tsv` for that sample: they show what
the call rests on. Check the
[troubleshooting page](https://github.com/duceppemo/BrucellaID/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/BrucellaID/wiki/FAQ) first.

## Development setup
```
git clone https://github.com/duceppemo/BrucellaID
cd BrucellaID
conda env create -f environment.yml
conda activate brucellaid
pip install --no-deps -e .
pytest --cov
```
BrucellaID itself needs nothing but the Python standard library; the environment is for MASH, BBTools and the test
tools. The end-to-end tests in `tests/test_integration.py` run MASH and Seal on a synthetic genome and on reads
tiled over the real MLST alleles, and are skipped if `mash` or `seal.sh` is not installed. The suite covers 100% of
the package, so please keep it there.

Lint with [ruff](https://docs.astral.sh/ruff/): `ruff check .`, or install the git hook with `pre-commit install`.
`ruff format` is **not** used: the code is wrapped by hand to read as prose, and reformatting it would rewrite
most of the files for nothing.
The tests run on GitHub Actions for Python 3.10, 3.12 and 3.14.

To see it work on data you can look at, run `bash examples/run_example.sh`: it builds a read set from the MLST
alleles of a known sequence type, runs BrucellaID on it and prints the reports.

## Reference data
The five files in `brucellaid/data/` are the reference data:

| File | What it is |
| --- | --- |
| `references.msh` | The MASH sketch of the 16 *Brucella* reference genomes, as vSNP assembled it |
| `references.tsv` | What each genome of the sketch is: species, biovar, strain, NCBI accessions |
| `mlst_alleles.fasta` | The 309 alleles of the nine MLST-9 loci |
| `mlst_profiles.tsv` | The 122 sequence types, as allele numbers |
| `mlst_st_species.tsv` | The species of each sequence type, counted over PubMLST's isolate database |

`scripts/make_mlst_data.py` rebuilds the two MLST files from PubMLST, `scripts/make_st_species.py` rebuilds the ST
to species table from PubMLST's isolate database, and `scripts/make_reference_data.py` rebuilds the sketch from
the accessions of `references.tsv`. Run them, check what changed, and say in the pull request which
day the data was downloaded; `wiki/Changelog.md` records it too. A sketch must keep MASH's defaults (k = 21,
1000 hashes), because the threshold of 150 shared hashes is read against that sketch size.

## Pull requests
* Branch from `main` and keep each pull request focused on one change.
* Add or update tests for any change in behaviour.
* Match the style of the surrounding code.
* Update the documentation in `wiki/` if the change affects users (see below), and add a line to
  `wiki/Changelog.md` under "Unreleased" (the section is created by the first change after a release).

## Documentation
The [wiki](https://github.com/duceppemo/BrucellaID/wiki) is maintained in the [`wiki/`](wiki/) folder of this
repository and published automatically when changes reach `main`. **Do not edit the wiki on GitHub directly**:
those changes would be overwritten. Edit the files in `wiki/` in your pull request instead.

## Releases (maintainers)
1. Update the version in `brucellaid/__init__.py` and `CITATION.cff` (`version` and `date-released`).
2. Rename "Unreleased" in `wiki/Changelog.md` to the version and date.
3. Refresh the example report of the wiki if the report changed: run BrucellaID over the four read sets of the
   [tutorial](wiki/Tutorial.md), then
   `python3 scripts/make_example_report.py <output folder>`, which writes `assets/example_report.html` with the
   machine it ran on removed. Never commit a report straight from a run: it names the host, the account and
   every path.
4. Commit, tag (`git tag -a vX.Y.Z`), push the commit and the tag, then create the GitHub release once the tests pass.
5. Publishing the GitHub release automatically uploads the package to
   [PyPI](https://pypi.org/project/brucellaid/) (`.github/workflows/publish.yml`, trusted publishing: no token
   needed; the release tag must match the package version), and Zenodo archives it and mints a DOI. Add that DOI to
   `CITATION.cff` and to the badge in `README.md`.
6. Bioconda: update `recipe/meta.yaml` (version and sha256 of the new GitHub tarball). Bioconda's bot usually opens
   the pull request to `bioconda-recipes` on its own for new releases; otherwise copy the recipe and open it by hand.
