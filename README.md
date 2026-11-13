<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo.svg" width="520" alt="BrucellaID: Brucella species identification and MLST from reads">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/duceppemo/BrucellaID/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/BrucellaID/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <!-- The branch is in the badge URL because Codecov still has "master" as this repository's default
       branch, and master is gone: without it the badge reads "unknown". -->
  <a href="https://app.codecov.io/gh/duceppemo/BrucellaID/tree/main"><img src="https://codecov.io/gh/duceppemo/BrucellaID/branch/main/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/BrucellaID/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/BrucellaID?cacheSeconds=3600" alt="Release"></a>
  <a href="https://pypi.org/project/brucellaid/"><img src="https://img.shields.io/pypi/v/brucellaid?cacheSeconds=3600" alt="PyPI"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="https://github.com/duceppemo/BrucellaID/blob/main/LICENSE"><img src="https://img.shields.io/github/license/duceppemo/BrucellaID" alt="License"></a>
  <a href="https://github.com/duceppemo/BrucellaID/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23845956"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23845956.svg" alt="DOI"></a>
</p>

Identify the *Brucella* species of a fastq sample by its [MASH](https://github.com/marbl/mash) distance to
429 curated reference genomes, covering all 30 named species of the clade, and type it with the nine-locus
[Whatmore MLST scheme](https://doi.org/10.1186/1471-2180-7-34), whose alleles are called from the reads with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools. Point it at a sample, or a folder of samples.

```
$ brucellaid -i reads/ -o results/
Sample      Species           SpeciesSource         ClosestReference                           SharedHashes  SequenceType  STSpecies
ERR2136546  Brucella canis    MASH+ST (arbitrated)  GCA_000530495.1_Oliveri_genomic.fna.gz     995/1000      20            Brucella canis
ERR3454415  Brucella abortus  MASH+ST (agree)       GCA_000742195.1_ASM74219v1_genomic.fna.gz  994/1000      1             Brucella abortus
```
(Two of the four [tutorial](https://github.com/duceppemo/BrucellaID/wiki/Tutorial) read sets, against the
429-genome default set. Columns elided to fit; `brucellaid.tsv` also holds `MashDistance`,
`AlleleProfile`, `Reads` and `Notes`.)

Both the reference sketch and the MLST database come with the package: nothing to download, index or
configure.

## Features
* **Two independent answers**: the species, from the closest of 429 curated *Brucella* genomes covering
  all 30 named species of the clade, and the
  sequence type, from the nine MLST alleles. Both come with the numbers they rest on, and neither is
  quietly rewritten by the other.
* **The ST decides the one pair a distance cannot.** *B. canis* and *B. suis* bv. 4 are 0.00024 apart —
  closer than two *B. suis* references are to each other, and in the 429-genome default set the nearest
  *canis*–*suis* pair, at 0.000264, is still the closest pair of different species there is — so which one a
  *canis* read set lands on first is noise. The two species share no sequence type at all, so where MASH
  leaves both within reach, the ST chooses, from a table of the species of every ST counted over the 3,403
  PubMLST isolates that carry one. A `SpeciesSource` column says what each call rested on, and an ST that contradicts a
  species MASH picked cleanly is reported, never obeyed.
* **Every answer comes with its evidence**: the distance to every reference, the reads behind every allele
  call and the runner-up that call had to beat, a note when a second reference was nearly as close, a note
  when the closest one is not close at all, and a tie reported as a tie rather than broken by a coin toss.
  Where a note would have to list a hundred near-identical genomes it gives a count instead, and
  `brucellaid_mash.tsv` still holds every distance.
* **A quality assurance report you can sign**: every run writes `brucellaid_report.html`, one
  self-contained file — no script, no stylesheet, nothing fetched from anywhere — holding the closest
  reference genomes and not only the winning one, how the species was decided, every allele call with the
  margin it had, the input files with their MD5, the exact command, and **the checksum of every reference
  data file the run measured against**, so two runs a year apart can be shown to have used the same data.
  One sample to a printed page, with a box for a reviewer's name, date and signature. `--pdf` writes the
  same thing as a PDF.
  **[See an example report](https://duceppemo.github.io/BrucellaID/example_report.html)**, or [the same thing as a PDF](https://duceppemo.github.io/BrucellaID/example_report.pdf).
* **Fits your pipeline**: a JSON report with everything in it, and a
  [MultiQC](https://docs.seqera.io/multiqc/custom_content) file, so `multiqc results/` folds the run into
  the report your pipeline already produces.
* **Batch mode**: point it at a folder; fastq files are found recursively, R1/R2 paired, and several samples
  run at a time. A sample that fails leaves the others readable.
* **Current reference data**: the MLST database is PubMLST's, 309 alleles and 122 sequence types, plus the
  species of 118 of them counted over PubMLST's isolate database — all rebuilt by scripts in the
  repository. Two sketches ship: the 429-genome clade set that 2.0 uses by default, and vSNP's 16 genomes
  under `--reference-set vsnp`, each with the species, biovar, strain and accession of every genome.
* **The whole clade, by default**: from 2.0 a run measures against 429 curated genomes
  covering **all 30 named species** — the former *Ochrobactrum* included — plus the unnamed lineages, so
  an atypical isolate is reported as *Brucella sp.* rather than forced onto the nearest classical species.
  Nothing in INSDC is *labelled* *B. papionis*, the thirtieth species — its genome sits there under
  *Brucella* sp., which is [how we found it](https://github.com/duceppemo/BrucellaID/blob/main/curation/papionis/README.md).
  [Validated](https://github.com/duceppemo/BrucellaID/wiki/Validation-at-scale) twice. On simulated reads
  of 2,022 NCBI assemblies: scored against what the sequences demonstrably are, the extended set got
  **every one of the 1,543 held-out named samples right** — five of them are mislabelled in NCBI, which the
  curation itself caught, and scored against the labels as written it is 1,538 of 1,543. And on **every
  usable public *Brucella* read set in ENA**, 4,024 real libraries, where it matches the submitter on 98.7%
  against the vSNP set's 96.2%, and wins all 95 of the disagreements that have a right answer without
  losing one.
* **Backwards compatible where it counts**: the MASH identification is vSNP's, down to the threshold of
  150 shared hashes (`--no-st-species` reports it alone), and `--xlsx` still writes the
  `stat_alignment_summary_<date>.xlsx` of the first version.
* **A short dependency list**: the Python standard library, plus MASH and BBTools. The two optional output
  formats, xlsx and PDF, are the only things that want a package, and neither is needed.
* **Tested**: 392 tests, 100% coverage, end to end on a synthetic genome and on reads tiled over the real
  MLST alleles.
* **Validated genome for genome against the vSNP set** (`--reference-set vsnp`), by a suite you can
  re-run yourself (`validation/run_validation.sh`): reads
  simulated from 13 of its 16 genomes, 7 of its 8 species, with every species and closest reference right
  and all 14 sequence types the one the genome carries — read off the assembly by exact match, which
  shares no code with the read-based caller. The *B. canis* genome is the one call the distance did not
  make on its own, and the ST made it; not one came back as a disagreement. The nearest neighbouring
  species shares 4 of 1000 hashes and comes back as not identified
  ([Validation suite](https://github.com/duceppemo/BrucellaID/wiki/Validation-suite)).

## What it does not do
Assembly, SNP calling, phylogeny, AMR, biovar typing and MLVA. It names the species and the sequence type,
and writes down what both answers rest on. **It does not handle biovars**, because a biovar is a
phenotype that does not follow the tree — no distance and no sequence type can settle one, so BrucellaID
does not guess ([the detail](https://github.com/duceppemo/BrucellaID/wiki/How-it-works#no-biovar)). Use
classical biotyping or [brucetyper](https://github.com/duceppemo/brucetyper) for that.

Two limits of the validation are worth knowing before you rely on it. **No mixed or metagenomic sample is
tested anywhere**: every measurement is a single isolate, so a sample holding two organisms is outside what
has been measured, and the loose-match note is the only warning you get. And **every validated run is
paired Illumina WGS** — nothing here measures a long-read library. The
[Validation](https://github.com/duceppemo/BrucellaID/wiki/Validation) pages state the limits of each
measurement under their own heading.

## Installation
```
conda create -n brucellaid -c conda-forge -c bioconda mash "gsl=2.7.0" bbmap "python>=3.10"
conda activate brucellaid
pip install brucellaid
```
MASH and BBTools (which provides `seal.sh`) are not on PyPI, so they come from bioconda; BrucellaID itself
needs nothing else. A bioconda package that installs all three in one command is
[awaiting review](https://github.com/bioconda/bioconda-recipes/pull/71040), so `conda install brucellaid`
does not find it yet. See [Installation](https://github.com/duceppemo/BrucellaID/wiki/Installation) for
the details and for running it from a clone.

## Quick start
```
brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/  # Paired-end reads
brucellaid -r1 sample.fastq.gz -o results/                            # Single-end reads
brucellaid -i reads/ -o results/ -p 4                                 # A folder, 4 samples at a time
brucellaid -i reads/ -o results/ --no-mlst                            # Species only
```
The results are printed and saved in `results/`: one row per sample, one row per locus, every MASH distance,
a JSON report, a MultiQC file, and `brucellaid_report.html` — a self-contained quality assurance report with
the evidence behind every call and a box to sign. To see it work on data with a known answer, run
`bash examples/run_example.sh`, which builds its own reads from the alleles of a sequence type.

**What the two answers mean**: the species is the species of the closest reference genome, except where a
reference of another species is just as close and the sequence type chooses between the two — the
`SpeciesSource` column says which happened. The sequence type is `new` when the nine alleles are called but
the combination is not in the table, and `ND` when a locus could not be called. There is no biovar, for
the reason given above. See
[How it works](https://github.com/duceppemo/BrucellaID/wiki/How-it-works).

## Documentation
The **[wiki](https://github.com/duceppemo/BrucellaID/wiki)** covers
[installation](https://github.com/duceppemo/BrucellaID/wiki/Installation), a
[tutorial](https://github.com/duceppemo/BrucellaID/wiki/Tutorial) on public read sets,
[usage and options](https://github.com/duceppemo/BrucellaID/wiki/Usage),
[output files](https://github.com/duceppemo/BrucellaID/wiki/Output-files),
[how it works](https://github.com/duceppemo/BrucellaID/wiki/How-it-works),
[validation](https://github.com/duceppemo/BrucellaID/wiki/Validation),
[troubleshooting](https://github.com/duceppemo/BrucellaID/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/BrucellaID/wiki/FAQ).

## Citing
If you use BrucellaID, please cite it together with the tools and the data it rests on:

> Duceppe M-O. BrucellaID: *Brucella* species identification and MLST from sequencing reads.
> Zenodo. https://doi.org/10.5281/zenodo.23845956

This DOI always points to the latest version; each release also has its own DOI, listed on
[Zenodo](https://doi.org/10.5281/zenodo.23845956).

> Ondov BD, Treangen TJ, Melsted P, et al. Mash: fast genome and metagenome distance estimation using
> MinHash. *Genome Biology* 2016;17:132. https://doi.org/10.1186/s13059-016-0997-x

> Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/

> Whatmore AM, Perrett LL, MacMillan AP. Characterisation of the genetic diversity of *Brucella* by
> multilocus sequencing. *BMC Microbiology* 2007;7:34. https://doi.org/10.1186/1471-2180-7-34

> Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the
> PubMLST.org website and their applications. *Wellcome Open Research* 2018;3:124.
> https://doi.org/10.12688/wellcomeopenres.14826.1

The 429 genomes of the default reference set are curated representatives of the assembly collection behind
the mashID *Brucella* pack, which has its own DOI:

> Duceppe M-O. *Brucella* database for mashID. figshare 2026. https://doi.org/10.6084/m9.figshare.33970168

GitHub's **"Cite this repository"** button gives the citation in APA and BibTeX formats, from
[`CITATION.cff`](https://github.com/duceppemo/BrucellaID/blob/main/CITATION.cff).

## Related
* [brucetyper](https://github.com/duceppemo/brucetyper) — *in silico* Bruce-ladder, Suis-ladder, AMOS and
  RB51 junction PCR for *Brucella*, from assemblies. Where BrucellaID gives the species from a distance and
  the sequence type from the alleles, brucetyper reads the assays a laboratory would run on a gel.
* [seqcounter](https://github.com/duceppemo/seqcounter) — count the occurrences of any reference sequences
  in reads or assemblies.
* [spoligotyper](https://github.com/duceppemo/spoligotyper) — *in silico* spoligotyping of the
  *M. tuberculosis* complex.

## Contributing
Bug reports, questions and pull requests are welcome: see [CONTRIBUTING.md](https://github.com/duceppemo/BrucellaID/blob/main/CONTRIBUTING.md).

## Author
Marc-Olivier Duceppe, Canadian Food Inspection Agency (CFIA): marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](https://github.com/duceppemo/BrucellaID/blob/main/LICENSE)
