# Changelog

## Unreleased
* The Zenodo DOI of the archived releases is recorded in `CITATION.cff`, in the README badge and in the
  citation on this wiki: [10.5281/zenodo.23845956](https://doi.org/10.5281/zenodo.23845956) resolves to
  the latest version.

## 1.0.0 (2026-10-31)
BrucellaID becomes an installable, tested package. It answers the same two questions as the first version —
which *Brucella* species a read set is, and which sequence type — and it answers the first one exactly as
before. What changed is everything around that, and how the MLST is called.

**Packaging and interface**
* Installable with `pip install brucellaid`, or from bioconda: a `brucellaid` command, and
  `python -m brucellaid`.
* **The two OLC dependencies are gone**: `accessoryFunctions` (OLCTools), which every function imported, and
  `MLSTsippr` (sipprverse), which ran the MLST. Neither installs from PyPI, so neither could be packaged.
  BrucellaID now needs nothing but the Python standard library, plus MASH and BBTools.
* **The MLST database is in the package**, so `-d` is gone. The input arguments are `-r1`/`-r2` for one
  sample or `-i` for a folder, instead of `-s`.
* The main table is **printed to standard output**, with the messages on standard error, so
  `brucellaid ... > results.tsv` works.
* A folder is searched recursively, R1/R2 files are paired, and samples run `-p` at a time.
* A sample that fails no longer stops the run: its error is a column of the report, and the exit code says
  that something failed.
* fasta input is **refused with a reason** instead of being typed badly: both steps count reads.

**Reports**
* `brucellaid.tsv` (one row per sample), `brucellaid_mlst.tsv` (one row per locus, with the reads behind
  each call and the runner-up), `brucellaid_mash.tsv` (every distance to every reference) and
  `brucellaid.json` (everything, including the versions of MASH, BBTools and Java).
* `--xlsx` still writes `stat_alignment_summary_<date>.xlsx`, with the same four columns and the vSNP
  species code in the `species` column, for the pipelines that read it.

**Identification**
* Unchanged in substance: the same 16 *Brucella* reference genomes from vSNP, the same MASH sketch, the same
  threshold of 150 shared hashes out of 1000. The reads of a sample are still sketched as one read set.
* What is new: the species, biovar, strain and NCBI accessions of each reference are in a table that ships
  with the package, the distance to **every** reference is reported, and a reference that is nearly as close
  as the closest one is called out in a note. The sketch and the table are checked against each other at
  startup.

**MLST**
* The alleles are called with **Seal** from BBTools instead of MLSTsippr's bowtie2 and samtools. Seal counts
  a read for every allele of a locus it shares the highest number of 25-mers with, and for no other
  (`ambiguous=all`, `clearzone=0`, no mismatch), so the allele the sample carries is counted for every read
  of the locus and the others only where they agree with it. The allele with the highest count wins, and a
  tie is reported rather than broken. See [How it works](How-it-works).
* Every call now comes with the reads behind it and the runner-up, so the margin is visible.
* **The database is PubMLST's**, downloaded by `scripts/make_mlst_data.py` on 2026-10-31 (the MLST (9 locus)
  scheme, data submitted up to 2024-12-31): 309 alleles and 122 sequence types, where the copy shipped with
  the old dependencies folder had 263 alleles and 81. That copy had also drifted: it still held a `gyrB_23`
  that PubMLST has deleted, and that sequence matched the *B. abortus* 9-941 genome better than the real
  gyrB_1, which turned ST 1 into a new profile. The 24 other gyrB alleles of the old copy are identical to
  PubMLST's.

**Around the code**
* 194 tests, 100% statement and branch coverage, run on Python 3.10, 3.12 and 3.14 by GitHub Actions.
* A [Validation](Validation) page, and `validation/run_validation.sh`, which reproduces it from public
  genomes.
* A [wiki](Home), a [tutorial](Tutorial) on public read sets, a bioconda recipe, `CITATION.cff`, and
  `examples/run_example.sh`, which runs in a few seconds on data it makes itself.
* The 193 MB of per-reference vSNP `script_dependents` folders are no longer in the repository: no code read
  them. They remain in the history and in the
  [v0.1 release](https://github.com/duceppemo/BrucellaID/releases/tag/v0.1), which is the last state of the
  clone-and-run version.
