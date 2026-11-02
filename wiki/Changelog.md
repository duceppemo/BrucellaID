# Changelog

## Unreleased
* **Concatenated gzip files are sized honestly.** The coverage estimate read the gzip trailer, which names
  only the *last* member of a file, so lane files joined with `cat` and the 64 KB blocks `bgzip` writes
  were estimated at a fraction of their size -- a deep bgzipped library could read as under 1x, keep the
  k-mer copy threshold at 2, and let the error k-mers back into the sketch. The head of the file is now
  decompressed member by member: several members mean the trailer is ignored outright, and a trailer that
  contradicts the measured compression ratio is ignored too.
* The reports are written in UTF-8 whatever the locale: a sample name with an accent no longer loses the
  reports of a finished run under `LC_ALL=C`.
* A gzip file with a corrupt stream inside makes the sample fail with the real error instead of crashing
  the run while its content was being sniffed for Seal.
* A coverage below 10x is written with one decimal, and a trace of reads reads "under 0.1x" -- never
  "roughly 0x coverage" beside a real base count.
* **Report tables are only as wide as their content**, in the HTML and the PDF alike, every column is
  left-aligned, and a rule opens and closes each table. A table stretched to the full page put unrelated
  amounts of air between its columns, and a right-aligned number sat away from its label. In the PDF, a
  column now gets the width its longest value asks for, and when a table is wider than the page only the
  columns with the most to say wrap.
* The `About ... coverage: sketching the k-mers seen at least N time(s)` log line now begins with the name
  of the sample it describes. With `--parallel`, two samples identify at once and their lines interleave
  with no way to tell whose coverage was whose.
* 335 tests. 1.2.0 claimed 100% coverage while four lines of the depth change had none; three tests close
  that, and the number is true again.

## 1.2.0 (2026-11-02)
**The k-mer copy threshold rises with the sequencing depth.** `mash dist -m 2` keeps every k-mer the read
set holds at least twice, and an error k-mer turns up about depth × error rate times, so past roughly 50x
the errors reach two copies on their own and crowd the true k-mers out of the bottom-1000 sketch. The
distance then grows with the depth instead of with the genetic distance: on 244 real libraries, the median
distance to the correct reference was 0.0219 at 100–150x where the assembly says 0.0001, and above about
240x the shared hashes fell under the 150 that `--min-hashes` asks for and the sample got **no call at
all**. The threshold is now derived from each sample — one more required copy per 20x of coverage, never
below 2, never above 20 — with the coverage estimated from the read files in a tenth of a second, without
decompressing them (the gzip trailer gives the size, the first thousand records give the bases per byte).

* Below about 50x the derived threshold is the same 2 as before and **nothing changes**, which is what the
  low end needs: a flat 10 would destroy a 5x sketch.
* At depth, everything does: on real libraries the medians moved from 666/1000 shared hashes at distance
  0.010637 to 966/1000 at 0.000843, and the 0.0005 close-match margin — which the sequence-type
  arbitration sits inside — means something again. The verdict changed on 14 of the 244 libraries, 6 of
  which had been reporting close matches off a margin of exactly zero.
* The reports say which threshold was used and the coverage it came from, and `brucellaid.json` carries
  `min_kmer_copies` and `estimated_coverage` under `closest_reference`.
* `--min-kmer-copies N` pins any fixed value, including the old 2.

**The validation covers depth now, and real reads.** Everything it simulated was 20x, where the derived
threshold returns the same 2 the old fixed value gave: the bug above would have sailed through it. Two
genomes are now simulated deep — *B. melitensis* 16M at 150x and *B. suis* 1330 at 250x — where the
derived threshold keeps 1000/1000 sketch hashes and the fixed 2 kept 860–902 and 728–748 in a six-sample
sweep, so `check_results.py` holds a floor of 950 and also checks the estimated coverage against the depth
the reads were simulated at and the threshold against the coverage. And the four public ENA libraries of
the [Tutorial](Tutorial) run next to the simulated samples: each must give the species its submitter gave
it, the reference and sequence type it is known to give, and the same species and reference as the
simulated sample of its own reference genome. All 22 samples pass, and the fix shows on real data:
SRR8550484 (~69x) moved from 422/1000 hashes at distance 0.0248 to 775/1000 at 0.0065, while ERR2136546
(47x) stayed exactly where it was. See [Validation](Validation).

**Upgrading:** no report or table changed shape. Results change only on read sets above about 50x, where
they were the ones being distorted: distances shrink toward what an assembly would give, shared hashes
rise, and a sample deep enough to have been `ND` can get its call back.

* 331 tests, 100% statement and branch coverage.

## 1.1.0 (2026-10-31)
**The biovar column is gone.** It never was a biovar call — it held the biovar of the closest
reference genome — and it was read as one. Three things rule out a real biovar from this data, and
none of them is a matter of trying harder:

* **Coverage.** The sketch holds *B. abortus* biovars 1 and 3 and no others, so **43%** of the
  biovar-labelled *B. abortus* isolates of the PubMLST isolate database are of a biovar it has no
  reference for. *B. suis* bv 5 is missing too.
* **Resolution.** The two *B. melitensis* **biovar 1** references of the sketch are **0.000756** apart,
  *further* than the biovar 2 and biovar 3 references are from each other (0.000706). Two genomes of one
  biovar are further apart than genomes of two biovars, so the ranking carries no biovar signal.
* **Biology.** A biovar is a phenotype and is not monophyletic: the best any rule from sequence type to
  biovar could manage is **66%** for *B. abortus* and **63%** for *B. melitensis*. ST 8 alone is biovar 3
  in 47 isolates, biovar 1 in 34 and biovar 2 in 22. Only *B. suis* is tidy, at 96%.

So `Biovar` is out of `brucellaid.tsv`, out of the MultiQC table, out of `brucellaid.json` and out of the
reports. What stays, because it is a fact about a reference rather than a claim about a sample: the biovar
of each reference genome in `brucellaid_mash.tsv`, and the biovars recorded for the sample's sequence type,
with their counts, in the quality assurance report. For a biovar, use classical biotyping or
[brucetyper](https://github.com/duceppemo/brucetyper). See [How it works](How-it-works#no-biovar).

**The report is easier to read.**

* A sample block lists the **five closest** reference genomes, plus any other near enough to the closest
  that the distance cannot separate them, instead of all sixteen. `brucellaid_mash.tsv` still holds every
  distance, and the report says how many it left out.
* **One sample to a printed page.** The fact lists are laid out in two columns and the tables across the
  page, which is what makes a sample fit.
* The summary is one table with a labelled row per group rather than one table per group, statuses are
  shown as coloured pills, sections are numbered, and the type, spacing and tables were reworked in both
  the light and the dark theme.
* The allele table merges the runner-up and its read count into one column.
* **The PDF was rebuilt to match**, rather than merely holding the same words: the same palette,
  numbered sections, one grouped summary table with coloured statuses, two columns of facts above each
  table, quieter tables with a rule under the header instead of a grid, a tinted box to sign, a rule
  above the footer, and one sample to a page there too.

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
* Installable with `pip install brucellaid` (and from bioconda once the recipe is merged): a
  `brucellaid` command, and `python -m brucellaid`.
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
