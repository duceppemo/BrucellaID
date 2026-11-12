<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/BrucellaID/main/assets/logo.svg" width="520" alt="BrucellaID: Brucella species identification and MLST from reads">
  </picture>
</p>

# BrucellaID

BrucellaID reads a *Brucella* fastq sample and answers two questions:

* **Which species is it?** The reads are compared with [MASH](https://github.com/marbl/mash) to 431 curated
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
$ brucellaid -i reads/ -o results/
Sample      Species           SpeciesSource         ClosestReference                           SharedHashes  SequenceType  STSpecies
ERR2136546  Brucella canis    MASH+ST (arbitrated)  GCA_000530495.1_Oliveri_genomic.fna.gz     995/1000      20            Brucella canis
ERR3454415  Brucella abortus  MASH+ST (agree)       GCA_000742195.1_ASM74219v1_genomic.fna.gz  994/1000      1             Brucella abortus
```
(Two of the four [tutorial](Tutorial) read sets, against the 431-genome default set. Columns elided to
fit; `brucellaid.tsv` also holds `MashDistance`, `AlleleProfile`, `Reads` and `Notes`.)

## Start here
* [Installation](Installation) — conda, pip, or a clone
* [Tutorial](Tutorial) — one real sample, start to finish
* [Usage](Usage) — every option
* [Output files](Output-files) — what each report holds

## Then
* [How it works](How-it-works) — what the two calls rest on, and where each of them stops
* [Validation](Validation) — 1,543 held-out genomes, 4,024 public read sets, and a 22-sample suite you can
  re-run, plus [the k-mer copy threshold](K-mer-threshold) and how it was measured
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
