<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo.svg" width="520" alt="BrucellaID: Brucella species identification and MLST from reads">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/duceppemo/BrucellaID/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/BrucellaID/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://codecov.io/gh/duceppemo/BrucellaID"><img src="https://codecov.io/gh/duceppemo/BrucellaID/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/BrucellaID/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/BrucellaID?cacheSeconds=3600" alt="Release"></a>
  <a href="https://pypi.org/project/brucellaid/"><img src="https://img.shields.io/pypi/v/brucellaid?cacheSeconds=3600" alt="PyPI"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/BrucellaID" alt="License"></a>
  <a href="https://github.com/duceppemo/BrucellaID/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23845956"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23845956.svg" alt="DOI"></a>
</p>

Identify the *Brucella* species of a fastq sample by its [MASH](https://github.com/marbl/mash) distance to 16
reference genomes, and type it with the nine-locus
[Whatmore MLST scheme](https://doi.org/10.1186/1471-2180-7-34), whose alleles are called from the reads with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools. Point it at a sample, or a folder of samples.

```
$ brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
Sample  Species           Biovar  ClosestReference           MashDistance  SharedHashes  SequenceType  AlleleProfile
sample  Brucella abortus  1       NC_006932-NC_006933.fasta  0.000024      999/1000      1             2-1-1-2-1-3-1-1-1
```

Both the reference sketch and the MLST database come with the package: nothing to download, index or
configure.

## Features
* **Two independent answers**: the species, from the closest of 16 *Brucella* reference genomes, and the
  sequence type, from the nine MLST alleles. Neither call corrects the other, and both come with the
  numbers they rest on.
* **Nothing is hidden**: the distance to every reference, the reads behind every allele call and the
  runner-up that call had to beat, a note when a second reference was nearly as close, and a tie reported as
  a tie rather than broken by a coin toss.
* **Batch mode**: point it at a folder; fastq files are found recursively, R1/R2 paired, and several samples
  run at a time. A sample that fails leaves the others readable.
* **Current reference data**: the MLST database is PubMLST's, 309 alleles and 122 sequence types, rebuilt by
  a script in the repository. The sketch is vSNP's 16 *Brucella* genomes, with the species, biovar, strain
  and accession of each.
* **Backwards compatible where it counts**: the species identification is vSNP's, down to the threshold of
  150 shared hashes, and `--xlsx` still writes the `stat_alignment_summary_<date>.xlsx` of the first
  version.
* **Nothing to configure**: only the Python standard library, plus MASH and BBTools.
* **Tested**: 194 tests, 100% coverage, end to end on a synthetic genome and on reads tiled over the real
  MLST alleles.
* **Validated** on 16 genomes with a known answer: every species, biovar and closest reference right, and
  all 14 sequence types the one the genome carries — read off the assembly by exact match, which shares no
  code with the read-based caller. The nearest neighbouring species shares 4 of 1000 hashes and comes back
  as not identified ([Validation](https://github.com/duceppemo/BrucellaID/wiki/Validation)).

## Installation
```
conda create -n brucellaid -c conda-forge -c bioconda brucellaid
conda activate brucellaid
```
MASH and BBTools (which provides `seal.sh`) are not on PyPI, so they come from bioconda; BrucellaID itself
needs nothing else. See [Installation](https://github.com/duceppemo/BrucellaID/wiki/Installation) for pip
and for running it from a clone.

## Quick start
```
brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/  # Paired-end reads
brucellaid -r1 sample.fastq.gz -o results/                            # Single-end reads
brucellaid -i reads/ -o results/ -p 4                                 # A folder, 4 samples at a time
brucellaid -i reads/ -o results/ --no-mlst                            # Species only
```
The results are printed and saved in `results/`: one row per sample, one row per locus, every MASH distance,
and a JSON report. To see it work on data with a known answer, run `bash examples/run_example.sh`, which
builds its own reads from the alleles of a sequence type.

**What the two answers mean**: the species is the species of the closest reference genome, and the biovar is
that reference's biovar, which is not the same as a biovar call — the *Brucella* biovar references differ by
very little. The sequence type is `new` when the nine alleles are called but the combination is not in the
table, and `ND` when a locus could not be called. See
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

GitHub's **"Cite this repository"** button gives the citation in APA and BibTeX formats, from
[`CITATION.cff`](CITATION.cff).

## Related
* [brucetyper](https://github.com/duceppemo/brucetyper) — *in silico* Bruce-ladder, Suis-ladder, AMOS and
  RB51 junction PCR for *Brucella*, from assemblies. Where BrucellaID gives the species from a distance and
  the sequence type from the alleles, brucetyper reads the assays a laboratory would run on a gel.
* [seqcounter](https://github.com/duceppemo/seqcounter) — count the occurrences of any reference sequences
  in reads or assemblies.
* [spoligotyper](https://github.com/duceppemo/spoligotyper) — *in silico* spoligotyping of the
  *M. tuberculosis* complex.

## Contributing
Bug reports, questions and pull requests are welcome: see [CONTRIBUTING.md](CONTRIBUTING.md).

## Author
Marc-Olivier Duceppe, Canadian Food Inspection Agency (CFIA): marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](LICENSE)
