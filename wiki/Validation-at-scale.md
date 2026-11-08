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

The two sets disagree on 146 runs. **The extended set is right on 91 of them and the vSNP set on none** —
there is no run in the sweep where the sixteen genomes gave a better answer. Seventy-two of those are
species the vSNP set cannot name, which come back `ND` at 0.17–0.26: 35 *B. anthropi*, 15
*B. intermedia*, and most of the rest of the former *Ochrobactrum* clade. Four of those runs re-measured
against the vSNP set under the 2.1.0 copy threshold sit at 0.168–0.244 and share 3 to 15 of 1000 hashes, so
this band is one the threshold change leaves where it was. Nineteen are wrong-species calls
the extended set corrects, most of them marine: the vSNP set scored **0 of 11** on *B. pinnipedialis*,
sending six of them to *B. ceti*, where the extended set gets 8 of 11 with the other three not identified.

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
  re-measured under the 2.1.0 threshold it is 0.0029–0.0088, and one — a *B. pituitosa* run at 0.0029 —
  now falls inside `--loose-distance` and so no longer carries the note. That is the clearest single piece
  of evidence that 0.005 is due a re-measurement.

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
fixed. Past `--loose-distance`, 0.005 by default, the call carries a note giving the distance and the
shared hashes, which takes the ten limit cases from one flagged to eight.

### The loose-match distance

The 0.005 of `--loose-distance` came out of this sweep rather than out of a guess. Of the 4,020 identified
runs, 195 sit past 0.005 — 4.9%, so the note stays rare — and 185 of those 195 agreed with the submitter
anyway, which is why the note does not touch the species call. What those 185 had in common is the useful
part: **every one of them shared fewer than 900 of the 1000 hashes**, median 696, against the 900-plus a
clean library shares. A loose match is usually a diluted library rather than a wrong species, and saying so
is more useful than either hiding it or acting on it.

> **These figures were measured before 2.1.0**, under the k-mer copy threshold of one copy per 20x. That
> threshold left error k-mers in the sketch and so inflated every read-derived distance; the distances
> above, and the 4.9% rate that follows from them, are therefore an upper bound on what the note does now.
> The "fewer than 900 hashes" half of the finding is the part that carries the conclusion, and it is the
> part least likely to move. Re-measuring the calibration against the current threshold is pending; see
> [The k-mer copy threshold](K-mer-threshold).

## What this does not establish

The sweep scores agreement with submitters, not truth, and the 20 mislabelled read sets it found are a
reminder that the denominator has errors in it too. Every one is listed with its evidence in
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md).

It also establishes nothing about read sets that are not paired Illumina WGS of a single isolate: that is
what the run selection filtered to, so no long-read, metagenomic or deliberately mixed library is measured
anywhere here.
