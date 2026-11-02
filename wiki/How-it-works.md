# How it works

A run does two independent things to the same reads: it measures the distance to 16 reference genomes, and
it calls the nine MLST alleles. The two answers are kept apart in the reports, with the numbers each one
rests on, and a third step puts them together only where it has to — for the one pair of species a sketch
of any size cannot separate. When they disagree anywhere else, both numbers are in the reports and the
disagreement is the finding.

## 1. Which species: the MASH distance to 16 reference genomes

The reads of a sample — R1 and R2 together, as one read set — are piped into
`mash dist -r -m N references.msh -`. MASH reduces the read set to the 1000 smallest hashes of its 21-mers
and compares that sketch to the sketch of each reference genome, which gives a distance and the number of
hashes the two share.

* `-m N` drops the k-mers the read set holds fewer than N times. Most of those are sequencing errors, and
  they would otherwise fill the sketch with noise. **N rises with the coverage** and BrucellaID derives it
  from the size of the input: one more required copy per 20x, from 2 to 20
  (`--min-kmer-copies` overrides it). An error k-mer turns up about `depth x error rate` times, so at 30x a
  single copy is almost certainly an error and two copies almost certainly real, while at 400x erroneous
  k-mers reach two copies routinely. A fixed `-m 2`, which versions up to 1.1.0 used, starts losing the
  true reference above roughly 240x: the shared hashes fall under `--min-hashes` and nothing is reported,
  on a sample that has more data than it needs. A fixed high value is no better, because it empties the
  sketch of a shallow sample. The report says which value was used.
* `-r` tells MASH the input is reads, so the genome size is estimated from k-mer multiplicity.
* The k-mer size and the sketch size are inherited from the reference sketch, which is why a custom sketch
  has to keep MASH's defaults.

The closest reference names the species, and the reference table also gives that genome's own biovar,
strain and NCBI accessions — facts about the reference, not calls about the sample. A reference is only reported when the reads share at least **150 of its 1000 hashes**
(`--min-hashes`). Below that the sample is reported as `ND`: there is some *Brucella* in it, or there is
not, but the closest of 16 genomes is not an answer. The first version of BrucellaID used the same
threshold.

**What this tells you, and what it does not.** The *Brucella* species are extraordinarily close — they are
over 99% identical — so a MASH distance separates them by a small margin. Read the species as the answer.
`brucellaid_mash.tsv` lists every distance, and a note on the sample says when another reference was
nearly as close.

<a name="no-biovar"></a>

### There is no biovar

BrucellaID reports a species and a sequence type. It does not report a biovar, and version 1.1 removed the
column that used to hold the biovar of the closest reference genome, because readers took it for a call.
Three things rule a biovar out, and none of them is fixable by trying harder:

| | |
| --- | --- |
| **Coverage** | The sketch holds *B. abortus* biovars 1 and 3 and no others. Of the biovar-labelled *B. abortus* isolates in the PubMLST isolate database, **43%** are of a biovar there is no reference for (bv 6, 2, 5, 4, 9, 7); for those, the nearest reference can only ever name the wrong one. *B. suis* bv 5 is missing too |
| **Resolution** | The two *B. melitensis* **biovar 1** references of the sketch are **0.000756** apart — *further* than the biovar 2 and biovar 3 references are from each other (0.000706). Two genomes of the same biovar are further apart than genomes of different biovars, so ranking by distance carries no biovar signal |
| **Biology** | A biovar is a phenotype (CO<sub>2</sub>, H<sub>2</sub>S, dye sensitivity, A/M sera, phage lysis) and is not monophyletic. Over the isolate database, the best any ST → biovar rule could manage is **66%** for *B. abortus* and **63%** for *B. melitensis*: ST 8 alone is biovar 3 in 47 isolates, biovar 1 in 34 and biovar 2 in 22. Only *B. suis* is tidy, at 96% |

What is still reported, because it is a fact rather than a call: `brucellaid_mash.tsv` gives the biovar of
each *reference genome*, and the quality assurance report shows the biovars recorded for the sample's ST,
with their counts, so you can see for yourself how little they settle.

**For a biovar**, use classical biotyping, or
[brucetyper](https://github.com/duceppemo/brucetyper), which reads the Bruce-ladder and Suis-ladder assays
*in silico*.

Of the 16 references, two (*B. ceti* Cudo and *B. suis* bv. 4 strain 40) come from vSNP under their own
names rather than from an accession; the other 14 are NCBI genomes, listed in `references.tsv`.

**One pair this cannot decide.** *B. canis* ATCC 23365 and *B. suis* bv. 4 strain 40 are **0.00024** apart
— the closest pair of the sixteen, closer than two *B. suis* references are to each other, and inside the
variation within either species. Which of the two a *canis* or a *suis* bv. 4 read set lands on first is
noise, and no sketch, at any size, with any number of references, will change that: they are that close. The
next closest cross-species pair in the sketch is 0.0007 apart, which MASH decides on its own. Step 3 is
about that one pair.

## 2. Which sequence type: the nine MLST alleles, called with Seal

The same reads go to `seal.sh` with the 309 alleles of the nine loci as references, and these settings:

```
k=25 hdist=0 rcomp=t maskmiddle=f ambiguous=all clearzone=0
```

Seal counts, for each read, which reference sequences share the **highest number of 25-mers** with it, and
counts the read for each of those and for no other. The consequence is what the call rests on:

* the allele the sample really carries shares every k-mer of every read of that locus, so it is counted for
  **all** of them;
* another allele of the same locus is only counted for the reads that fall where the two agree, and loses
  every read that covers a position where they differ.

So for each locus the allele with the highest count wins, and `brucellaid_mlst.tsv` prints the runner-up and
its count next to it, which is the margin. Two alleles with the same count are reported as a tie
(`ND` for that locus) rather than broken by a coin toss: at a low depth, two alleles can agree everywhere
the reads happen to cover. A locus with fewer than `--min-reads` reads (5) is not called either.

`hdist=0` is deliberate. The alleles of a locus differ by a few bases — `trpE_3` and `trpE_5` differ at one
position of 486 — and a k-mer allowed a mismatch would match both.

The nine allele numbers are then looked up in the profile table:

| Result | What it means |
| --- | --- |
| A number | That profile is in the table: the sample is that ST |
| `new` | All nine alleles were called, but the combination is not in the table |
| `ND` | At least one locus has no call. The `Notes` column says which, and why |

**Novel alleles.** An allele that is not in the database is called as its nearest neighbour, which usually
shows up as a lower count and a close runner-up, and often as a `new` profile. BrucellaID does not assemble,
so it cannot tell you the new sequence: assemble the sample and submit the locus to
[PubMLST](https://pubmlst.org/organisms/brucella-spp) if you suspect one.

**Depth.** The alleles are 422 to 589 bp long. At 20x, a locus is covered by roughly 80 reads and the margin
over the runner-up is comfortable; below about 5x the ties start, and they are reported as ties rather than
guessed. Coverage at the very ends of a locus matters as much as in the middle: two alleles differing at a
single position near one end can only be told apart by a read that reaches it.

## 3. Which species, when the distance cannot say: the sequence type

The two answers are different kinds of evidence, not sharper and blunter versions of the same one. A sketch
measures the whole genome and smears the two species together; the nine MLST loci measure 4,396 bp and keep
them apart. In the PubMLST *Brucella* isolate database, *B. canis* is ST 20, 21 and 106 and nothing else,
*B. suis* is fourteen other STs, and **not one ST is shared between them**.

So the species is decided like this:

1. MASH names the species of the closest reference genome.
2. If a reference of **another** species is within 0.0005 of the closest one, MASH has not chosen, and the
   species it named is held open.
3. The ST then chooses among the species MASH left within reach — and only among those.

The result is in the `SpeciesSource` column:

| `SpeciesSource` | What happened |
| --- | --- |
| `MASH` | One species was within reach, and no usable ST. The species is the closest reference's |
| `MASH+ST (agree)` | One species was within reach, and the ST says the same thing |
| `MASH+ST (arbitrated)` | Two species were within reach, and the ST confirms the closest one |
| `ST (arbitrated)` | Two species were within reach, and the ST names the other one. **The species reported is the ST's**, and so is the vSNP code, which is taken from the closest reference of that species |
| `MASH (unresolved)` | Two species were within reach and the ST could not choose. The species is the closest reference's, and the note says why the ST was no help |
| `MASH (conflict)` | The ST names a species that was not within reach at all. **The MASH call stands** — a thousand min-hashes against nine loci — and the note says to check the sample |

An ST only speaks when the isolates carrying it agree on a species: at least **95%** of them
(`--min-st-purity`), over at least **3** isolates (`--min-st-isolates`). The species labels in the isolate
database are the submitters', and some of them are wrong — the same thing that turns up in NCBI assemblies
— so an ST whose isolates disagree is left out of the decision rather than settled by a majority of three.
Of the 118 STs with a species, 105 are seen in one species only and 113 clear the 95% bar; the exceptions
are things like ST 1, which is *B. abortus* in 463 of 466 isolates and *B. melitensis* in 3. STs recorded
only as "*Brucella* sp." say nothing and are ignored.

`--no-st-species` turns the whole step off and reports the closest reference's species, whatever else was
just as close.

**This does not make the ST a species call.** It breaks one tie, in one direction, inside the set MASH
already narrowed to. An ST that disagrees with a species MASH picked cleanly never wins; it is reported, and
it is worth looking into.

## What the reference data is

| File | Source |
| --- | --- |
| `references.msh`, `references.tsv` | The *Brucella* part of [vSNP](https://github.com/USDA-VS/vsnp)'s reference set, 16 genomes |
| `mlst_alleles.fasta`, `mlst_profiles.tsv` | The MLST (9 locus) scheme of [PubMLST](https://pubmlst.org/organisms/brucella-spp), downloaded with `scripts/make_mlst_data.py` |
| `mlst_st_species.tsv` | The species of each ST, counted over the 3,426 isolates of the PubMLST *Brucella* isolate database, with `scripts/make_st_species.py` |

The MLST scheme itself is [Whatmore, Perrett and MacMillan (2007)](https://doi.org/10.1186/1471-2180-7-34):
nine loci, 4,396 bp in all.

## The extended reference set

The 16 genomes of the default set cover 8 species — the classical six plus *B. microti* and *B. neotomae*.
The 2020 taxonomy moved the genus *Ochrobactrum* into *Brucella*, and sequencing keeps finding isolates of
the atypical lineages, so the package also ships a second set: 429 curated genomes covering all 29 named
species of the clade, selected with

```bash
brucellaid -i reads/ -o results/ --reference-set extended
```

With it, a *B. anthropi* or *B. intermedia* sample is named instead of coming back `ND`, a marine-mammal
or vole-adjacent species lands on its own references instead of the nearest classical one, and an isolate
of a lineage that has no name yet is reported as *Brucella sp.* — an honest "inside the clade, outside the
named species" — rather than forced onto the closest named species. What does **not** change: the species
is still the species of the closest genome, the ST still arbitrates where two species are within reach
(*B. canis* sits inside *B. suis* whichever set measures it), and there is still
[no biovar](How-it-works#no-biovar).

The set is derived from the 2,022 *Brucella* (NCBI taxon 234) assemblies behind the
[mashID *Brucella* pack](https://doi.org/10.6084/m9.figshare.33970168), dereplicated at a MASH distance of
0.0003 and then curated: genomes whose NCBI label their own sequence contradicts were dropped, unlabelled
twins of named species were dropped, and one cluster got the label its neighbours agree on. Every decision,
with its reason, is
[`scripts/extended_references_curation.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/extended_references_curation.tsv),
and [Validation](Validation#at-scale-the-extended-reference-set) shows what the set does on simulated
reads of all 2,022 assemblies. The default set remains the default: it is vSNP's, every answer in the
field so far was measured against it, and the extended set is one release old.

## Why not an assembly?

Both steps count reads. MASH sketches a read set and filters the k-mers seen fewer times than the copy
threshold derived from the coverage, and the allele calling
works because a locus is covered by many reads that each vote for an allele; a contig votes once, and a
contig holding several loci would be assigned to one of them. A fasta file is refused rather than quietly
mis-typed. For assemblies, [brucetyper](https://github.com/duceppemo/brucetyper) does *in silico* PCR
instead.
