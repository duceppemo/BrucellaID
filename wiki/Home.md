<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo.svg" width="520" alt="BrucellaID: Brucella species identification and MLST from reads">
  </picture>
</p>

# BrucellaID

BrucellaID reads a *Brucella* fastq sample and answers two questions:

* **Which species is it?** The reads are compared with [MASH](https://github.com/marbl/mash) to 432 curated
  *Brucella* genomes covering all 30 named species of the clade, *Ochrobactrum* included, and the closest
  one names the species and the reference itself. `--reference-set vsnp` measures against vSNP's 16 genomes
  instead, which was the default before 2.0 ([How it works](How-it-works#the-reference-sets)).
* **Which sequence type is it?** The nine alleles of the
  [Whatmore MLST-9 scheme](https://doi.org/10.1186/1471-2180-7-34) are called from the reads with
  [Seal](https://sourceforge.net/projects/bbmap/) from BBTools, and looked up in the
  [PubMLST](https://pubmlst.org/organisms/brucella-spp) profile table.

Both the reference sketch and the MLST database come with the package: there is nothing to download, index or
configure. Where two species are too close for a sketch to choose between them — *B. canis* and
*B. suis* bv. 4 are 0.00024 apart — the sequence type chooses, and the `SpeciesSource` column says so.

```
$ brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
Sample  Species             SpeciesSource         ClosestReference           SharedHashes  SequenceType  STSpecies
sample  Brucella abortus    MASH+ST (agree)       NC_006932-NC_006933.fasta  999/1000      1             Brucella abortus
canis   Brucella canis      MASH+ST (arbitrated)  NC_010103-NC_010104.fasta  1000/1000     20            Brucella canis
```
(Columns elided to fit; `brucellaid.tsv` also holds `MashDistance`, `AlleleProfile`, `Reads` and `Notes`.)

## Start here
* [Installation](Installation) — conda, pip, or a clone
* [Tutorial](Tutorial) — one real sample, start to finish
* [Usage](Usage) — every option
* [Output files](Output-files) — what each report holds

## Then
* [How it works](How-it-works) — what the two calls rest on, and where each of them stops
* [Validation](Validation) — 4,042 public read sets, and 16 genomes read for read
* [Troubleshooting](Troubleshooting) and the [FAQ](FAQ)
* [Changelog](Changelog)

## Citing
> Duceppe M-O. BrucellaID: *Brucella* species identification and MLST from sequencing reads.
> Zenodo. https://doi.org/10.5281/zenodo.23845956

That DOI resolves to the latest version. Cite it together with the tools and data BrucellaID rests on: MASH,
BBTools, the [MLST scheme](https://doi.org/10.1186/1471-2180-7-34) and
[PubMLST](https://doi.org/10.12688/wellcomeopenres.14826.1). The full list, ready to paste, is in the
[README](https://github.com/duceppemo/BrucellaID#citing), and GitHub's "Cite this repository" button gives
it in APA and BibTeX, from
[CITATION.cff](https://github.com/duceppemo/BrucellaID/blob/main/CITATION.cff).
