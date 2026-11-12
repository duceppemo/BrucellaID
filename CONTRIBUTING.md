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
Neither is in `environment.yml`, so `pip install ruff pre-commit` first.
`ruff format` is **not** used: the code is wrapped by hand to read as prose, and reformatting it would rewrite
most of the files for nothing.
The tests run on GitHub Actions for Python 3.10, 3.12 and 3.14.

To see it work on data you can look at, run `bash examples/run_example.sh`: it builds a read set from the MLST
alleles of a known sequence type, runs BrucellaID on it and prints the reports.

## Reference data
The reference data in `brucellaid/data/`, plus the logo the reports carry:

| File | What it is |
| --- | --- |
| `references_extended.msh` | The MASH sketch of the 431 curated genomes a run identifies against by default from 2.0: all 30 named species of the clade, plus the lineages with no name yet |
| `references_extended.tsv` | What each of those genomes is. Its `code` is the species epithet (`abortus5`, `anthropi`), not a vSNP code |
| `references.msh` | The MASH sketch of the 16 *Brucella* reference genomes, as vSNP assembled it. `--reference-set vsnp` selects it; it was the default before 2.0 |
| `references.tsv` | What each genome of that sketch is: vSNP species code, species, biovar, strain, NCBI accessions |
| `mlst_alleles.fasta` | The 309 alleles of the nine MLST-9 loci |
| `mlst_profiles.tsv` | The 122 sequence types, as allele numbers |
| `mlst_st_species.tsv` | The species of 118 of the 122 sequence types, counted over the 3,403 isolates of PubMLST's isolate database that carry one |
| `logo.svg`, `logo.png` | The logo of the reports, written by `assets/make_logo.py` |

`scripts/make_mlst_data.py` rebuilds the two MLST files from PubMLST and `scripts/make_st_species.py` rebuilds
the ST to species table from PubMLST's isolate database. The two sketches have a script each, and they are not
interchangeable: `scripts/make_extended_references.py` rebuilds the 431-genome default from the accessions of
`references_extended.tsv`, and `scripts/make_reference_data.py` rebuilds the 16-genome vSNP sketch from the
accessions of `references.tsv`. Run them, check what changed, and say in the pull request which
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
3. Refresh the example report in `docs/` if the report changed: fetch the four read sets of the
   [tutorial](wiki/Tutorial.md), then run
   `python3 scripts/make_example_report.py <folder of the tutorial reads>`, which analyses them itself and writes
   `docs/example_report.html` and `docs/example_report.pdf` with the computer they ran on left out, refusing
   to write either if a value it had to remove survived. Never commit a report straight from a run: it names
   the host, the account and every path. `docs/` is published at https://duceppemo.github.io/BrucellaID/ by
   GitHub Pages.
4. Push the commit to `main` and wait for `tests.yml` to go green **before** tagging: that workflow does not
   run on tags, so a tag is only as tested as the branch commit it points at. Then
   `git tag -a vX.Y.Z`, push the tag, and create the GitHub release.
5. Publishing the GitHub release automatically uploads the package to
   [PyPI](https://pypi.org/project/brucellaid/) (`.github/workflows/publish.yml`, trusted publishing: no token
   needed; the release tag must match the package version), and Zenodo archives it and mints a DOI. Add that DOI to
   `CITATION.cff` and to the badge in `README.md` in a follow-up commit, since neither exists until the release
is published.
6. Bioconda: update `recipe/meta.yaml` (version and sha256 of the new GitHub tarball). Bioconda's bot usually opens
   the pull request to `bioconda-recipes` on its own for new releases, but only once the recipe is already in
that repository; until the first one is merged, copy the recipe and open the pull request by hand.
