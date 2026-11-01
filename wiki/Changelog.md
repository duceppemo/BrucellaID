# Changelog

## Unreleased
**A quality assurance report.** Every run now writes `brucellaid_report.html`: one self-contained file —
no script, no stylesheet, no font, nothing fetched from anywhere — holding what a result rests on rather
than only the result.

* **`brucellaid_report.html`**, written on every run (`--no-html` skips it): a summary with the samples
  grouped by what they were called and a box for a reviewer's name, date and signature; the notes,
  warnings and errors one row each; a block per sample with **every** distance MASH measured (the closest
  highlighted, the ones near enough to it shaded), how the species was decided and what the sequence type
  contributed, the nine allele calls with the margin each had over the allele it beat, and the input files
  with their size, date and MD5; what the numbers mean; and the run itself — operator, computer, command,
  versions, and **the MD5 of every reference data file**. It prints to clean pages as it stands.
* **`brucellaid_mqc.json`**, also on every run: a MultiQC custom content table, so `multiqc results/` shows
  a BrucellaID section with the species, what decided it, the sequence type and the status of each sample.
* **`--pdf`** writes the same report as a PDF, with "Page 3 of 9" in every footer. It needs reportlab,
  which is an optional extra: `pip install 'brucellaid[pdf]'`. BrucellaID itself still installs with
  nothing but the standard library, and the HTML report holds the same information.
* **`--no-md5`** skips the checksums of the input files, which cost a second or two per gigabyte, and
  **`--operator NAME`** puts a name other than the user account on the report.
* The report carries the BrucellaID logo, inlined as svg in the HTML (one drawing, themed with CSS
  variables, so it follows the light and the dark page) and drawn from a png in the PDF. Both ship in
  the package, written by `assets/make_logo.py`.
* An example report, [HTML](https://duceppemo.github.io/BrucellaID/example_report.html) and
  [PDF](https://duceppemo.github.io/BrucellaID/example_report.pdf): the four read sets of the
  [Tutorial](Tutorial), rendered with the computer they ran on left out by
  `scripts/make_example_report.py`, and published with GitHub Pages from `docs/` so that a browser
  renders them instead of showing their source.
* Each sample now has a status — `ok`, `warning` (not identified, no ST, a `new` profile, or the two kinds
  of evidence disagreeing) or `failed` — in the reports and in the MultiQC table.
* `brucellaid.json` gains, per sample, the `status`, the input files with their checksums and, in
  `species_call`, `st_not_used`: why the sequence type was left out of the species call when it was. Its
  run block gains the operator, the user, the computer, the operating system, the working directory and
  the checksummed `reference_data`.

**The sequence type now decides the species MASH cannot.** *B. canis* ATCC 23365 and *B. suis* bv. 4 strain
40 are 0.00024 apart in the sketch — the closest pair of the sixteen, closer than two *B. suis* references
are to each other — so which one a *canis* or a *suis* bv. 4 read set lands on first is noise. The two
species share no sequence type at all, so the ST settles it.

* New data file `brucellaid/data/mlst_st_species.tsv`: the species of each ST of the scheme, counted over
  the 3,426 isolates of the PubMLST *Brucella* isolate database, with the share of them that agree and the
  biovars recorded. 118 STs, 105 of them seen in one species only. Built by the new
  `scripts/make_st_species.py`.
* Where a reference of another species is within 0.0005 of the closest one, the species is no longer the
  closest reference's by default: a sequence type whose isolates agree chooses between the species MASH left
  within reach, and the biovar and the vSNP code follow the species chosen. An ST that contradicts a species
  MASH picked cleanly is reported but never obeyed. See [How it works](How-it-works).
* Two new columns in `brucellaid.tsv`: **`SpeciesSource`** (`MASH`, `MASH+ST (agree)`,
  `MASH+ST (arbitrated)`, `ST (arbitrated)`, `MASH (unresolved)`, `MASH (conflict)`) and **`STSpecies`**.
  `brucellaid.json` gains a `species_call` object. Anything reading the main table by column *number* has to
  be updated; by column name, nothing changes.
* New options: `--no-st-species` (decide the species as 1.0.0 did), `--min-st-purity` (0.95),
  `--min-st-isolates` (3) and `--st-species FILE`.
* An ST is only used when at least `--min-st-purity` of at least `--min-st-isolates` isolates agree on a
  species, and never when that species is "*Brucella* sp.". The labels are the submitters' and some are
  wrong, so a split ST is left out of the decision rather than settled by a majority.
* `CONTRIBUTING.md` listed the allele and sequence type counts of vSNP's old copy of the database, 263 and
  81; they are PubMLST's 309 and 122.
* The wiki front page carries the logo, in the light and the dark variant, as the README does.
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
