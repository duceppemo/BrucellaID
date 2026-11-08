# How it works

A run does two independent things to the same reads: it measures the distance to 432 reference genomes, and
it calls the nine MLST alleles. The two answers are kept apart in the reports, with the numbers each one
rests on, and a third step puts them together only where it has to — for the one pair of species a sketch
of any size cannot separate. When they disagree anywhere else, both numbers are in the reports and the
disagreement is the finding.

<a name="the-mash-distance"></a>

## 1. Which species: the MASH distance to 432 reference genomes

The reads of a sample — R1 and R2 together, as one read set — are piped into
`mash dist -r -m N references.msh -`. MASH reduces the read set to the 1000 smallest hashes of its 21-mers
and compares that sketch to the sketch of each reference genome, which gives a distance and the number of
hashes the two share.

* `-m N` drops the k-mers the read set holds fewer than N times. Most of those are sequencing errors, and
  they would otherwise fill the sketch with noise. **N rises with the coverage** and BrucellaID derives it
  from the size of the input: one more required copy per 5x, from 2 to 48
  (`--min-kmer-copies` overrides it). An error k-mer turns up about `depth x error rate` times, so at 30x a
  single copy is almost certainly an error and two copies almost certainly real, while at 400x erroneous
  k-mers reach two copies routinely. Measured over 63 real libraries, a fixed `-m 2` — which versions up to
  1.1.0 used — keeps a median 896 of the 1000 hashes below 100x, 477 from 100x to 250x, and 174 above that;
  on four libraries past 1400x it falls under the 150 that `--min-hashes` asks for, so nothing is reported
  at all on a sample that has more data than it needs. A fixed high value is no better, because it empties
  the sketch of a shallow sample. The report says which value was used, and
  [The k-mer copy threshold](K-mer-threshold) is how the rule was arrived at.
* `-r` tells MASH the input is reads, so the genome size is estimated from k-mer multiplicity.
* The k-mer size and the sketch size are inherited from the reference sketch, which is why a custom sketch
  has to keep MASH's defaults.

The closest reference names the species, and the reference table also gives that genome's own biovar,
strain and NCBI accessions — facts about the reference, not calls about the sample. A reference is only reported when the reads share at least **150 of its 1000 hashes**
(`--min-hashes`). Below that the sample is reported as `ND`: there is some *Brucella* in it, or there is
not, but the closest genome of the set is not an answer. The first version of BrucellaID used the same
threshold.

**What this tells you, and what it does not.** The *Brucella* species are extraordinarily close — they are
over 99% identical — so a MASH distance separates them by a small margin. Read the species as the answer.
`brucellaid_mash.tsv` lists every distance, and a note on the sample says when another reference was
nearly as close.

What the distance does **not** tell you is how much to trust the ranking. It orders the references; it
carries no confidence, and against the 432-genome default set the second line is usually close behind the
first, because a sample's nearest neighbours are mostly other genomes of its own species. A note that
another reference was nearly as close is therefore the normal case, not a warning — the field to alert on
is `Status`, not `Notes`.

A second note says when the closest reference is **not close at all**. Being the nearest of several
hundred genomes is not the same as being a good match, and the distance says which: past
`--loose-distance`, 0.001 by default, the call carries a note with the distance and the shared hashes.
It does not change the species, because a loose match does not predict a wrong one: over every public
*Brucella* read set in ENA, agreement with the submitter is flat in distance — 99.8% within 0.0005 and
between 95% and 100% in every band above it, out to past 0.02. The threshold itself is the 95th percentile
of the distance of agreeing calls, which is 0.00066, and
[the derivation](Validation-at-scale#the-loose-match-distance) has both tables. A loose match has three
causes, and
the note names all three: the reads hold more than *Brucella*; the sample belongs to a lineage with no near
relative in the sketch; or a minority *Brucella* population is being filtered out along with the
sequencing errors, which `--min-kmer-copies` either side of the derived value will show — see
[The k-mer copy threshold](K-mer-threshold). Any of the three is worth knowing, and none is visible from
the species name alone.

<a name="no-biovar"></a>

### There is no biovar

BrucellaID reports a species and a sequence type. It does not report a biovar, and version 1.1 removed the
column that used to hold the biovar of the closest reference genome, because readers took it for a call.
Three things rule a biovar out, and none of them is fixable by trying harder:

| | |
| --- | --- |
| **Coverage** | A biovar is recorded for only **7 of the 432** genomes of the default set, because the NCBI records of the rest carry none: *B. abortus* bv. 5, *B. melitensis* bv. 1 and 2, and *B. suis* bv. 1, 2 and 3. Of the 240 *B. abortus* isolates with a numbered biovar in the PubMLST isolate database, **93%** are of a biovar no reference carries; for those the nearest reference could only ever name the wrong one. vSNP's 16 genomes are better covered here and still miss 42% of them |
| **Resolution** | The two *B. suis* **biovar 2** genomes of the default set are **0.00163** apart, while the *B. suis* biovar 3 genome sits **0.00150** from one of them — *closer* than the two biovar 2 genomes are to each other. Two genomes of the same biovar are further apart than genomes of different biovars, so ranking by distance carries no biovar signal. The same held in vSNP's 16: its two *B. melitensis* biovar 1 references are 0.000756 apart, its biovar 2 and 3 references 0.000706 |
| **Biology** | A biovar is a phenotype (CO<sub>2</sub>, H<sub>2</sub>S, dye sensitivity, A/M sera, phage lysis) and is not monophyletic. Over the isolate database, the best any ST → biovar rule could manage is **66%** for *B. abortus* and **63%** for *B. melitensis*: ST 8 alone is biovar 3 in 47 isolates, biovar 1 in 34 and biovar 2 in 22. Only *B. suis* is tidy, at 96% |

What is still reported, because it is a fact rather than a call: `brucellaid_mash.tsv` gives the biovar of
each *reference genome*, and the quality assurance report shows the biovars recorded for the sample's ST,
with their counts, so you can see for yourself how little they settle.

**For a biovar**, use classical biotyping, or
[brucetyper](https://github.com/duceppemo/brucetyper), which reads the Bruce-ladder and Suis-ladder assays
*in silico*.

The 432 genomes of the default set are all NCBI assemblies, named in `references_extended.tsv` by the GCA
accession they were downloaded under. Of vSNP's 16, two (*B. ceti* Cudo and *B. suis* bv. 4 strain 40) come
from vSNP under its own names rather than from an accession; the other 14 are NCBI genomes, listed in
`references.tsv`.

**One pair this cannot decide.** *B. canis* and *B. suis* are closer to each other than either is to
anything else, in both sets: *B. canis* ATCC 23365 and *B. suis* bv. 4 strain 40 are **0.00024** apart in
vSNP's 16, and the closest *canis*–*suis* pair of the 432 is **0.000264** — the closest pair of different
species in either sketch, closer than two *B. suis* references are to each other, and inside the variation
within either species. Which of the two a *canis* or a *suis* read set lands on first is noise, and no
sketch, at any size, with any number of references, will change that: they are that close. The three
closest pairs of different species in the default set are all *canis*–*suis*; the closest that is not,
*B. ceti* against *B. pinnipedialis*, is 0.00058. Step 3 is about the *canis*–*suis* pair.

<a name="the-mlst-alleles"></a>

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

**Depth, and what a comfortable margin looks like.** The alleles are 422 to 589 bp long, so at 20x a locus
is covered by roughly 80 reads. Across the 20x simulated samples of the
[validation suite](Validation-suite), every called allele had 64 to 150 reads behind it (median 108) and
beat its runner-up by 14 reads at the narrowest, median 34 — that is the scale "close behind" should be read
against. Below about 5x the ties start, and they are reported as ties rather than
guessed. Coverage at the very ends of a locus matters as much as in the middle: two alleles differing at a
single position near one end can only be told apart by a read that reaches it.

<a name="the-species-call"></a>

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
| `references_extended.msh`, `references_extended.tsv` | 432 curated genomes of the whole clade, the set a run uses by default from 2.0: dereplicated representatives of the assemblies behind the [mashID *Brucella* pack](https://doi.org/10.6084/m9.figshare.33970168), rebuilt with `scripts/make_extended_references.py` |
| `references.msh`, `references.tsv` | The *Brucella* part of [vSNP](https://github.com/USDA-VS/vsnp)'s reference set, 16 genomes, reached with `--reference-set vsnp` and rebuilt with `scripts/make_reference_data.py` |
| `mlst_alleles.fasta`, `mlst_profiles.tsv` | The MLST (9 locus) scheme of [PubMLST](https://pubmlst.org/organisms/brucella-spp), downloaded with `scripts/make_mlst_data.py` |
| `mlst_st_species.tsv` | The species of each ST, counted over the 3,403 isolates of the PubMLST *Brucella* isolate database that carry an ST, with `scripts/make_st_species.py` |

The MLST scheme itself is [Whatmore, Perrett and MacMillan (2007)](https://doi.org/10.1186/1471-2180-7-34):
nine loci, 4,396 bp in all.

## The reference sets

The package ships two, and **the extended set is the default from 2.0**:

```bash
brucellaid -i reads/ -o results/                      # 432 genomes, all 30 named species
brucellaid -i reads/ -o results/ --reference-set vsnp # vSNP's 16 genomes, 8 species
```

vSNP's 16 genomes cover 8 species — the classical six plus *B. microti* and *B. neotomae*. The 2020
taxonomy moved the genus *Ochrobactrum* into *Brucella*, and sequencing keeps finding isolates of the
atypical lineages, so a sample of any of the other 22 species could only land on a wrong neighbour or on
`ND`. The extended set is 432 curated genomes covering all 30 named species of the clade.

With it, a *B. anthropi* or *B. intermedia* sample is named instead of coming back `ND`, a marine-mammal
or vole-adjacent species lands on its own references instead of the nearest classical one, and an isolate
of a lineage that has no name yet is reported as *Brucella sp.* — "inside the clade, outside the
named species" — rather than forced onto the closest named species. What does **not** change: the species
is still the species of the closest genome, the ST still arbitrates where two species are within reach
(*B. canis* sits inside *B. suis* whichever set measures it), and there is still
[no biovar](How-it-works#no-biovar).

**Why the default moved.** Measured over every usable public *Brucella* read set in ENA — 4,042 runs — the
two sets disagree on 146, and the extended set is right on 91 of them and vSNP's 16 on none. There is no
run in that archive where the smaller set gave a better answer, and across the 2,018 *B. melitensis* and
1,118 *B. abortus* runs that carry most real traffic the two return the same call, so the extra genomes add
coverage without disturbing what the 16 already did well.
[Validation](Validation-at-scale#every-usable-public-read-set-in-ena) has the numbers.

**Why `vsnp` is still there.** It is what every answer from this tool before 2.0 was measured against, it
is the set vSNP itself uses, and its `code` column carries vSNP's own species codes (`ab1`, `mel1`), which
the extended table does not — see [Usage](Usage) if a pipeline reads those.

The set is derived from the 2,022 *Brucella* (NCBI taxon 234) assemblies behind the
[mashID *Brucella* pack](https://doi.org/10.6084/m9.figshare.33970168), dereplicated at a MASH distance of
0.0003 and then curated: genomes whose NCBI label their own sequence contradicts were dropped, unlabelled
twins of named species were dropped, and one cluster got the label its neighbours agree on. Every decision,
with its reason, is
[`scripts/extended_references_curation.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/extended_references_curation.tsv),
and [Validation](Validation-at-scale#simulated-reads-of-2022-assemblies) shows what the set does on simulated
reads of all 2,022 assemblies and on every usable public *Brucella* read set in ENA.

One of the thirty was nearly missed. Nothing in INSDC is labelled *B. papionis*, the species described
from baboons in 2014: NCBI's assembly database answers nothing for it and so does ENA's. Its genome is
nevertheless there, filed as "*Brucella* sp. NVSL 07-0026", and the first pass of the curation dropped it
as a redundant unnamed twin of *B. microti* 0.0020 away — a species-level gap in the classical core, not
a within-species one. It is now the species' reference, relabelled on its sequence type and its
distances, which
[`curation/papionis/README.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/papionis/README.md)
lays out. Only one *B. papionis* strain has ever been sequenced, so the species is covered rather than
tested.

## Why not an assembly?

Both steps count reads, so a fasta file is refused rather than quietly mis-typed: see
[the FAQ](FAQ#can-it-read-assemblies), which also gives the two things to use on an assembly instead.
