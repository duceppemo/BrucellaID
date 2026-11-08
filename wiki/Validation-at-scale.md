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
filtered to a proper `_1`/`_2` pair between 10 MB and 6 GB — 4,044 runs, 1.6 TB of reads, of which 4,042
were analysed. The two that are missing could not be downloaded after two passes; none failed to analyse.
These are other people's samples with other people's labels, so a disagreement is a question rather than
automatically an error.

| Reference set | Matches the submitter | Does not |
| --- | --- | --- |
| vSNP set (16 genomes, 8 species) | 3,890 (96.2%) | 152 |
| Extended set (432 genomes) — **the default from 2.0** | 3,990 (98.7%) | 52 |

The two sets disagree on 146 runs. **The extended set is right on 94 of them and the vSNP set on none** —
there is no run in the sweep where the sixteen genomes gave a better answer. (Counted as the runs where the
extended call matches the submitter and the vSNP call does not; it is 100 before setting aside six whose
submitter label is only *Brucella* sp. and so has no species-level answer to be right about.)

The 146 split two ways. **89 are runs the vSNP set could not name at all**, coming back `ND` at 0.15–0.30
while the extended set names them: 35 *B. anthropi*, 21 *B. intermedia*, 7 *B. pituitosa*,
7 *B. pseudogrignonensis* and the rest of the former *Ochrobactrum* clade. Four of those runs re-measured
against the vSNP set under the 2.1.0 copy threshold sit at 0.168–0.244 and share 3 to 15 of 1000 hashes, so
this band is one the threshold change leaves where it was. **The other 57 are runs both sets named but
named differently**, and the extended set corrects 19 of them, most of them marine: the vSNP set scored
**0 of 11** on *B. pinnipedialis*, sending six of them to *B. ceti*, one to *B. melitensis*, one to
*B. microti* and leaving three unidentified, where the extended set gets 8 of 11 with the other three not
identified.

All 146 are listed in
[`scripts/public_sweep_disagreements.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/public_sweep_disagreements.tsv),
and `scripts/score_public_sweep.py` derives every figure in these two paragraphs from it. They are shipped
because three of those figures had drifted: they were computed before the sweep's retry pass filled in its
last runs and were never recomputed, which put the head-to-head win at 91 rather than 94, the unnameable
runs at 72 rather than 89, and *B. intermedia* at 15 rather than 21.

What the extended set's 52 non-matches are, one by one:

* **20 are mislabelled read sets.** For 19 of them the MASH distance and the sequence type — which share
  no code and no reference data — both contradict the submitted species and agree with each other; for the
  twentieth, both reference sets agree independently. Re-measured under the 2.1.0 threshold, several sit at
  0.0001 or closer to the species the sequence says they are, and two sit at zero. Eleven come from a
  single submission.
* **22 contain little or no *Brucella*.** Re-measured under the 2.1.0 copy threshold they share 0 to 6 of
  1000 hashes with every genome of the clade, four of them exactly zero, against the 980 or more a real
  *Brucella* library shares with its own species. BrucellaID named none of them, which is what the
  threshold of 150 shared hashes is for. Under the threshold the sweep ran with, the range was 0 to 146 and
  one run sat 4 hashes under the cut; the re-derivation took that borderline case to 2 of 1000, because
  asking one copy per 5x of a deep mixed library removes a trace population along with the errors.
* **10 are the honest limits of the set.** Eight are former *Ochrobactrum* species whose nearest reference
  is far enough away to be a thin reference side rather than a wrong answer; one is the *B. papionis* read
  set, which had no reference to match at the time; one is an *inopinata*-adjacent lineage that landed on
  an unnamed representative. The sweep measured that band at 0.005–0.02; across the four of them
  re-measured under the 2.1.0 threshold it is 0.0029–0.0088. All four are still flagged, because
  `--loose-distance` came down with them — at the old 0.005 the *B. pituitosa* run at 0.0029 would have
  been missed, which is what the recalibration below exists to prevent.

The sweep ran before the *B. papionis* genome was added. Re-running that one read set with it in the set
reports *B. papionis* at 0.00010, which would make the extended figure 3,991 of 4,042; most of the other
runs were not re-measured, so the table above is what was actually measured. The 48 read sets behind
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md)
have since been re-measured under the 2.1.0 threshold: every species call is unchanged, and every distance
is smaller.

**Not one of the 52 is a case where the set held the right reference and named the wrong species.** Scored
against what the sequences demonstrably are, the extended set is right on 4,032 of 4,042.

What this does establish is that the extended set costs nothing on the species the vSNP set already
covers: across 2,018 *B. melitensis* and 1,118 *B. abortus* runs — the bulk of real traffic — the two sets
return the same answer, so the extra 416 genomes add coverage without disturbing the two most common
species in the archive.

This sweep also found the gap it needed to: for the diverse part of the clade there was no distance-based
caution at all. A sample 0.0089 from its nearest *B. anthropi* was reported as a plain `MASH` call with the
same confidence as one at 0.0003, and only one of those ten limit cases carried any flag. That is now
fixed. Past `--loose-distance`, 0.001 by default, the call carries a note giving the distance and the
shared hashes.

### The loose-match distance

`--loose-distance` came out of this sweep rather than out of a guess, and the rule is one line: **the
threshold is the 95th percentile of the distance of calls that agree with the submitter**, taken to a
round number. Past that, the closest reference is further away than a correctly identified library normally
sits, which is worth a note.

It was **0.005 in 2.0 and is 0.001 from 2.1.0**, because the re-derived
[k-mer copy threshold](K-mer-threshold) shrank every read-derived distance by about sevenfold and the old
value was calibrated on the inflated ones. Re-measured against the current threshold:

| | 2.0 | 2.1.0 |
| --- | --- | --- |
| 95th percentile of agreeing calls | 0.00464 | **0.00066** |
| threshold taken from it | 0.005 | **0.001** |
| share of the 4,020 identified runs it flags | 4.9% | 4.0% |
| of those, still agreeing with the submitter | 96% | 98% |

So the note stays as rare as it was, and the old value left where it was would have missed the
thin-reference cases it exists to catch — the *B. pituitosa* run two bullets up, at 0.0029, among them.
The choice is not on a knife edge: anything from 0.0005 to 0.002 flags 2.8% to 6.6% of runs with 96% to 98%
of them still agreeing, and ordinary libraries are nowhere near it, the 3,226 runs already inside 0.001
before now sitting inside 0.00046.

**What says the note should not change the species call** is a separate measurement, and the only one here
that is not a restatement of the threshold: agreement with the submitter is flat in distance.

| distance | runs | agree with submitter |
| --- | --- | --- |
| 0 – 0.0005 | 3,756 | 99.8% |
| 0.0005 – 0.001 | 103 | 100% |
| 0.001 – 0.005 | 85 | 95.3% |
| 0.005 – 0.02 | 53 | 100% |
| past 0.02 | 23 | 100% |

A loose match is a diluted library or a thin reference side; it does not predict a wrong species, which is
why it is reported and not acted on.

> Earlier versions of this page offered a second fact as corroboration — that the agreeing calls past the
> threshold all shared fewer than 900 of their 1000 hashes. **That was not evidence.** MASH derives the
> distance from the shared-hash count, so the two are one measurement in different units: over these 4,042
> runs, 317 distinct shared-hash values map to 317 distinct distances and not one maps to two. "Past 0.005"
> and "sharing fewer than about 820" are the same statement, so the shared-hash version is entailed by the
> threshold rather than supporting it. The flat agreement rate above is what does the work.

The 207 re-measured runs are in
[`scripts/loose_distance_calibration.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/scripts/loose_distance_calibration.tsv),
and `scripts/recalibrate_loose_distance.py` prints both tables from it. Measuring all 4,020 again was not
practical, so the population was binned by each run's old distance, a sample drawn from each bin, and every
row carries its bin's weight; the percentiles are computed on the weighted sample. The method was checked
by re-deriving the 2.0 figures from the original sweep with it, which gave back the published 0.00464 and
4.9%. The script also applies the same weighting to the pre-2.1 distances it ships, which gives 0.00496 and
4.8% — the gap from those published figures is the sampling error of 207 runs standing for 4,020.

## What this does not establish

The sweep scores agreement with submitters, not truth, and the 20 mislabelled read sets it found are a
reminder that the denominator has errors in it too. Every one is listed with its evidence in
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md).

It also establishes nothing about read sets that are not paired Illumina WGS of a single isolate: that is
what the run selection filtered to, so no long-read, metagenomic or deliberately mixed library is measured
anywhere here.
