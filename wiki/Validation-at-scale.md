# Validation at scale

The [validation suite](Validation-suite) is reproducible and small. These two measurements are
neither: they are reported here because they are the only way to say how a reference set behaves
across the whole clade and across other people's libraries.

## Simulated reads of 2,022 assemblies

The [extended reference set](How-it-works#the-reference-sets) was validated separately, at a scale
the repository's script is not meant for: 150 bp paired-end reads simulated at 20x (0.2% substitutions)
from **all 2,022 *Brucella* assemblies** of the collection it was curated from, identified against three
reference sets on the same machine. Of the 2,022, 95 are labelled *Brucella sp.* in NCBI and have no right
answer to score, and 432 are genomes of the extended set itself, where a self-match tells us nothing. Those
two exclusions overlap by 48, the unnamed genomes that are themselves references, so what is left is
**1,543 held-out samples carrying a named species** (2,022 − 95 unnamed − 384 named references), and the
numbers below are theirs.

| Reference set | Correct species | Wrong | Not identified |
| --- | --- | --- | --- |
| vSNP set (16 genomes, 8 species) | 1,451 (94.0%) | 25 | 67 |
| mashID pack as published (389 genomes) | 1,532 (99.3%) | 11 | 0 |
| Extended set (432 genomes, curated) — **the default from 2.0** | 1,538 (99.68%) | 5 | 0 |

The vSNP set's misses are the species it does not carry — *B. pinnipedialis* and the porpoise lineage
of *B. ceti* land on *B. ceti* or *B. microti*, the whole former *Ochrobactrum* clade comes back `ND` —
and all but two were flagged `unresolved` or `conflict`. The two that were not are the reason the
extended set exists: two "*B. abortus*" samples reported as *B. melitensis* with no warning at all.

Both turned out to be **mislabelled in NCBI**, and so are the extended set's five "wrong" calls, every
one a genome the curation had already thrown out of the reference list for contradicting its own label:
the "*B. abortus*" that is 0.0007 from *B. melitensis* and carries the melitensis ST 8, the two atypical
"*B. abortus* bv. 1" BCB genomes, and a submission batch whose *B. suis* and *B. inopinata* sequences are
swapped (its "*B. suis* bv. 4 str. 40" is 0.022 from every *B. suis* and 0.0004 from *B. inopinata* BO1).
Each of those records is tracked with its evidence in
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md),
together with the public read sets whose submitted species the distance and the sequence type both
contradict, so they can be reported rather than only excluded.

Scored against what those five genomes demonstrably are, the extended set got **every one of the 1,543
held-out named samples right**, including all 69 *B. canis*, all 24 *B. ceti* and all 8
*B. pinnipedialis* — the boundaries the distance alone cannot draw. Of the 95 *Brucella sp.* samples, 53
were reported as *Brucella sp.*, 40 as the named species they sit next to (36 of them with a note saying
how close the call was), and the two genomes that are nowhere near the clade came back `ND`.

The set was revised after this was measured, and the figures above are recomputed against the set as it
now stands: three genomes were added (the *B. papionis* reference and two unnamed lineages the curation
had wrongly dropped as twins), and six genomes of the former *Ochrobactrum* group were relabelled by ANI
to each species' type strain. Every one of the 2,022 assemblies was re-measured. Only four samples are
nearer to an added genome than to the reference they had matched, and **all four are labelled *Brucella*
sp. in NCBI**, so none was ever in the named universe. The six relabels move no held-out named sample
either, because no such sample has one of them as its closest reference. The two neighbours with the
least margin, a *B. microti* and a *B. pinnipedialis* read set, return the same species, closest
reference and distance before and after.

What this does not establish is everything the [validation suite](Validation-suite#what-it-does-not-establish) does not, plus one thing: the reads were simulated from the very
collection the set was curated from, so a lineage absent from GenBank in October 2026 is absent from both
sides of the measurement. That is what the next section is for.

## Every usable public read set in ENA

Simulated reads cannot show how a reference set behaves on real libraries, so both sets were run over
**every usable public *Brucella* read set in ENA**: each paired Illumina WGS run under NCBI taxon 234,
filtered to a proper `_1`/`_2` pair between 10 MB and 6 GB — 4,044 runs, 1.6 TB of reads. These are other
people's samples with other people's labels, so a disagreement is a question rather than automatically an
error.

**4,024 of the 4,044 were analysed.** Twenty could not be completed: 18 never finished downloading a
second mate after two passes, and two of the manifest's runs were already unavailable. A run with only one
mate is left out rather than analysed half-depth, because the k-mer copy threshold is derived from
coverage and a half pair would read as a shallow library.

| Reference set | Matches the submitter | Does not |
| --- | --- | --- |
| vSNP set (16 genomes, 8 species) | 3,871 (96.2%) | 153 |
| Extended set (432 genomes) — **the default from 2.0** | 3,972 (98.7%) | 52 |

The two sets disagree on 168 runs. **The extended set is right on 95 of them and the vSNP set on none** —
there is no run in the sweep where the sixteen genomes gave a better answer. (Counted as the runs where
the extended call matches the submitter and the vSNP call does not; it is 101 before setting aside six
whose submitter label is only *Brucella* sp. and so has no species-level answer to be right about.)

The 168 split two ways. **89 are runs the vSNP set could not name at all**, coming back `ND` at 0.15–0.30
while the extended set names them: 35 *B. anthropi*, 21 *B. intermedia*, 7 *B. pituitosa*,
7 *B. pseudogrignonensis* and the rest of the former *Ochrobactrum* clade. **The other 79 are runs both
sets named but named differently**, and the extended set corrects 20 of them, most of them marine: the
vSNP set scored **0 of 11** on *B. pinnipedialis*, sending seven of them to *B. ceti*, one to
*B. microti* and leaving three unidentified, where the extended set gets 8 of 11 with the other three not
identified.

All 168 are listed in
[`scripts/public_sweep_disagreements.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/public_sweep_disagreements.tsv),
and `scripts/score_public_sweep.py` derives every figure in these two paragraphs from it.

What the extended set's 52 non-matches are:

* **23 share too little with the clade to be named.** They hold 0 to 8 of 1000 hashes against every genome
  of the set, four of them exactly zero, where a correctly identified library shares a median 994 and all
  but 2% of them more than 900. BrucellaID named none of them, which is what the threshold of 150 shared
  hashes is for. Twenty-two
  of the twenty-three come from a handful of submissions and are tracked with their evidence in
  [`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md);
  only `SRR35815234` is not yet listed there.
* **29 are named, but not as the submitter labelled them.** For most the MASH distance and the sequence
  type — which share no code and no reference data — both contradict the submitted species and agree with
  each other. Several sit at 0.0001 or closer to the species the sequence says they are, two of them at
  zero. Eleven come from a single submission. The rest are the honest limits of the set: former
  *Ochrobactrum* species whose nearest reference is far enough away to be a thin reference side rather
  than a wrong answer, and an *inopinata*-adjacent lineage that lands on an unnamed representative.

**Not one of the 52 is a case where the set held the right reference and named the wrong species.**

What this does establish is that the extended set costs nothing on the species the vSNP set already
covers: across 2,038 *B. melitensis* and 1,122 *B. abortus* runs — the bulk of real traffic — the two sets
return the same answer on all but two, so the extra 416 genomes add coverage without disturbing the two
most common species in the archive.

### The loose-match distance

`--loose-distance` came out of this sweep rather than out of a guess, and the rule is one line: **the
threshold is the 95th percentile of the distance of calls that agree with the submitter**, taken to a
round number. Past that, the closest reference is further away than a correctly identified library
normally sits, which is worth a note.

It was **0.005 in 2.0 and is 0.001 from 2.1.0**, because the re-derived
[k-mer copy threshold](K-mer-threshold) shrank the read-derived distances — the typical run by about
2.4x and the tail that this threshold is drawn from by nearly 7x — and the old value had been calibrated
on the inflated ones. Measured over all 4,001 identified runs:

| | 2.0 | 2.1.0 onwards |
| --- | --- | --- |
| 95th percentile of agreeing calls | 0.00464 | **0.00068** |
| threshold taken from it | 0.005 | **0.001** |
| share of identified runs it flags | 4.9% | 3.4% |
| of those, still agreeing with the submitter | 96% | 92% |

So the note stays rare, and the old value left where it was would have missed the thin-reference cases it
exists to catch — a *B. pituitosa* run at 0.0029 among them. The choice is not on a knife edge: every
threshold from 0.0005 to 0.002 flags 2.0% to 8.2% of runs, and ordinary libraries are nowhere near it.

**What says the note should not change the species call** is a separate measurement, and the only one here
that is not a restatement of the threshold: how often a call at a given distance is still the submitter's
species.

| distance | runs | agree with submitter |
| --- | --- | --- |
| 0 – 0.0005 | 3,673 | 99.5% |
| 0.0005 – 0.001 | 191 | 99.5% |
| 0.001 – 0.005 | 87 | 95.4% |
| 0.005 – 0.02 | 32 | 81.2% |
| past 0.02 | 18 | 94.4% |

**Agreement falls with distance rather than staying flat.** A loose match therefore does carry some signal
that the species is wrong — a run between 0.005 and 0.02 matches its submitter 81% of the time against
99.5% for a close one. The signal is weak, and the worst band is still mostly right, which is why the note
is reported rather than acted on; but it is not nothing, and an earlier version of this page said
agreement was flat across every band. That was measured on a 207-run sample and the full sweep does not
support it.

> Two cautions on this table. The denominator is the submitter's label, not the truth, so some of the
> disagreement in the loose bands is the submitter being wrong rather than BrucellaID — the former
> *Ochrobactrum* runs cluster there, and they are thin-reference cases, not errors. And the bands beyond
> 0.005 hold 32 and 18 runs, so a few samples move them by several points.

The 4,001 identified runs are in
[`scripts/loose_distance_calibration.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/loose_distance_calibration.tsv),
and `scripts/recalibrate_loose_distance.py` prints both tables from it. Until the sweep was re-run under
the current threshold that file held a weighted 207-run sample instead, which put the 95th percentile at
0.00066 against the 0.00068 measured here — close enough that the threshold chosen from the sample needed
no change.

### The three-file runs, and a staging mistake worth recording

Eighteen of these runs publish **three** fastq files rather than two. That is how ENA presents a run
holding both paired and unpaired reads: `_1` is the unpaired ones and `_2`/`_3` are the mates, which carry
identical read counts. On `ERR10445439` the three hold 16,700, 1,163,751 and 1,163,751 reads.

The first pass of this sweep staged every run as `_1` plus `_2`, so for these eighteen it fed the orphan
file in as R1 and the real R1 as R2. The two then held different numbers of reads, Seal refused the pair,
and all eighteen came back with an error and no sequence type. An earlier version of this page reported
that as a property of the archive — "unequal numbers of reads in their two mates, as published" — which
was wrong. All three files match the md5 ENA publishes because there is nothing wrong with them.

Re-run on `_2`/`_3`, **all eighteen type cleanly**: seventeen as *B. melitensis* ST 8 and one,
`SRR9164843`, as *B. microti* ST 53, each agreeing with the species MASH called. The species figures above
are unaffected, because MASH had named all eighteen correctly either way and their distances moved by at
most one shared hash.

The lesson is about the denominator rather than the tool: a sweep that globs `_1` and assumes `_2` is its
mate will mis-pair exactly the runs that carry unpaired reads, and it will do so silently. BrucellaID now
warns when a `_3` file sits beside a pair it has formed, which is the only point at which the layout is
visible.

## What this does not establish

The sweep scores agreement with submitters, not truth, and the 20 mislabelled read sets it found are a
reminder that the denominator has errors in it too. Every one is listed with its evidence in
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md).

It also establishes nothing about read sets that are not paired Illumina WGS of a single isolate: that is
what the run selection filtered to, so no long-read, metagenomic or deliberately mixed library is measured
anywhere here.
