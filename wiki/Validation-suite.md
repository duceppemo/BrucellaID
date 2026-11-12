# Validation suite

**22 samples whose answer is known**: reads simulated from 16 genomes, of which 13 are in the sketch under
test, plus two of those genomes again at depths where a fixed k-mer copy threshold breaks, and four real
public read sets that must keep giving the answer they are known to give.

This first suite measures **the vSNP reference set**, not the extended set that 2.0 identifies against by
default, so its expected answers are vSNP sketch IDs and the script passes `--reference-set vsnp`. Its
negative control is *B. anthropi* ATCC 49188, which has to come back `ND` precisely because the vSNP set
does not carry it. Sixteen genomes are simulated from, but they are not the sixteen of the sketch: **13 of
the sketch's 16 are tested**, covering 7 of its 8 species, and the other three samples are a *B. abortus*
genome deliberately left out of the sketch, the *B. anthropi* control and an *E. coli* control. The three
sketch genomes never tested are *B. ceti* Cudo and *B. suis* bv. 4 strain 40, which vSNP distributes under
its own names so the script cannot download them, and the *B. neotomae* reference.

What the default set does is measured separately, in [Validation at scale](Validation-at-scale).

Reproduce it with:

```bash
bash validation/run_validation.sh      # about 40 minutes and 1.5 GB, all of it in validation/data
REAL=0 bash validation/run_validation.sh   # simulated reads only: no 540 MB from the ENA
```

The script downloads the genomes listed in
[`validation/genomes.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/validation/genomes.tsv) from
NCBI, simulates 150 bp paired-end reads from each (random fragments, 0.2% substitutions, one seed per
sample; 20x unless the table fixes the depth), downloads the read sets listed in
[`validation/reads.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/validation/reads.tsv) from the
ENA, runs BrucellaID over the folder and checks the result.

**The sequence type is not listed anywhere** for the simulated samples: `check_results.py` reads it off
each genome by looking for an exact match of every allele of the package, on both strands, and compares
that with what BrucellaID called from the reads. The two answers share no code: one searches an assembly
for a string, the other counts k-mer votes in reads. For the real read sets there is no assembly to read,
so the species is what the submitter labelled the library and the reference and sequence type are the
answer the library gave once and must keep giving — a regression check rather than an independent truth.

It also checks that no sample comes back as `MASH (conflict)`: these samples are what they say they are, so
the distance and the sequence type must not name species that are nowhere near each other.

## What came back

```
sample                   kind            species              reference                      hashes     k-mer copies  ST called  ST expected    species from          result
-----------------------  --------------  -------------------  -----------------------------  ---------  ------------  ---------  -------------  --------------------  ------
abortus_9-941            simulated       Brucella abortus     NC_006932-NC_006933.fasta      1000/1000  4 (20x)       1          1              MASH+ST (agree)       ok
abortus_BER              simulated       Brucella abortus     NZ_CP007682-NZ_CP007683.fasta  1000/1000  4 (20x)       6          6              MASH+ST (agree)       ok
melitensis_16M           simulated       Brucella melitensis  NC_003317-NC_003318.fasta      1000/1000  4 (20x)       7          7              MASH+ST (agree)       ok
melitensis_ATCC23457     simulated       Brucella melitensis  NC_012441-NC_012442.fasta      1000/1000  4 (20x)       8          8              MASH+ST (agree)       ok
melitensis_Ether         simulated       Brucella melitensis  NZ_CP007760-NZ_CP007761.fasta  999/1000   4 (20x)       9          9              MASH+ST (agree)       ok
melitensis_BwIM_SOM_36b  simulated       Brucella melitensis  NZ_CP018508-NZ_CP018509.fasta  1000/1000  4 (20x)       12         12             MASH+ST (agree)       ok
suis_1330                simulated       Brucella suis        NC_017251-NC_017250.fasta      999/1000   4 (20x)       14         14             MASH+ST (agree)       ok
suis_ATCC23445           simulated       Brucella suis        NC_010169-NC_010167.fasta      1000/1000  4 (20x)       15         15             MASH+ST (agree)       ok
suis_686                 simulated       Brucella suis        NZ_CP007719-NZ_CP007718.fasta  1000/1000  4 (20x)       17         17             MASH                  ok
canis_ATCC23365          simulated       Brucella canis       NC_010103-NC_010104.fasta      1000/1000  4 (20x)       20         20             MASH+ST (arbitrated)  ok
ovis_ATCC25840           simulated       Brucella ovis        NC_009505-NC_009504.fasta      1000/1000  4 (20x)       13         13             MASH+ST (agree)       ok
ceti_TE10759-12          simulated       Brucella ceti        NC_022905-NC_022906.fasta      999/1000   4 (20x)       26         26             MASH+ST (agree)       ok
microti_CCM4915          simulated       Brucella microti     NC_013119_NC_013118.fasta      1000/1000  4 (20x)       53         53             MASH+ST (agree)       ok
abortus_2308             simulated       Brucella abortus     NC_006932-NC_006933.fasta      993/1000   4 (20x)       5          5              MASH+ST (agree)       ok
anthropi_ATCC49188       simulated       ND                   B-REF-BS4-40.fasta             4/1000     6 (29x)       ND         ND (no locus)                        ok
coli_MG1655              simulated       ND                   B-REF-BS4-40.fasta             0/1000     6 (28x)       ND         ND (no locus)                        ok
melitensis_16M_150x      simulated 150x  Brucella melitensis  NC_003317-NC_003318.fasta      1000/1000  30 (151x)     7          7              MASH+ST (agree)       ok
suis_1330_250x           simulated 250x  Brucella suis        NC_017251-NC_017250.fasta      1000/1000  48 (254x)     14         14             MASH+ST (agree)       ok
SRR19520329              real            Brucella melitensis  NZ_CP007760-NZ_CP007761.fasta  992/1000   11 (53x)      11         11             MASH+ST (agree)       ok
ERR3454415               real            Brucella abortus     NC_006932-NC_006933.fasta      995/1000   14 (70x)      1          1              MASH+ST (agree)       ok
SRR8550484               real            Brucella suis        NC_010169-NC_010167.fasta      989/1000   14 (69x)      15         15             MASH+ST (agree)       ok
ERR2136546               real            Brucella canis       NC_010103-NC_010104.fasta      996/1000   9 (47x)       20         20             MASH+ST (arbitrated)  ok

22 sample(s), all as expected.
```

The `k-mer copies` column is the threshold BrucellaID derived for each sample, with the coverage it
estimated from the files; the biovar of each simulated genome is in `genomes.tsv`. The two ~29x entries are
the two non-*Brucella* genomes: their 4.7 Mb at 20x of their own length is 29x of a *Brucella* genome,
which is the only unit the estimate is for.

Every simulated *Brucella* sample now keeps 999 or 1000 of its 1000 hashes, which is the right answer for
reads simulated from a genome the sketch holds, and the four real libraries keep 989 to 996 where they
kept 773 to 979 before. That is the re-derived copy threshold:
[The k-mer copy threshold](K-mer-threshold) is where it comes from.

## What that establishes

**The identification names the right genome.** Each of the genomes that are in the sketch came back as
itself, at every depth, with 999 or 1000 of its 1000 hashes shared, and with the species the reference
table gives it.

**And the right species for a genome that is not in the sketch.** *B. abortus* 2308 is not one of the 16.
It was identified as *Brucella abortus* through 9-941, sharing 993 hashes.

**The k-mer copy threshold rises with the depth, and the sketch stays clean.** The 150x and 250x samples
are simulated from genomes the sketch holds, so their true distance is zero — any distance they do
report is error k-mers crowding the sketch. With the derived threshold (30 at 150x, 48 at 250x) both keep
**1000/1000** hashes at distance 0.000000. With the fixed threshold of 2 that BrucellaID used before 1.2.0,
a sweep of four genomes at 150x and two at 250x kept 860–902 and 728–748 hashes, at distances of 0.0025 to
0.0082 — on references the reads were simulated from. `check_results.py` holds the line at 950 of 1000,
checks the estimated coverage against the depth the reads were simulated at, and checks the threshold
against the coverage, so a regression has to get past all three.

**Real reads agree with simulated reads.** The four ENA libraries are the ones the [Tutorial](Tutorial)
runs. Each gave the species its submitter gave it, the sequence type it is known to carry, and the same
reference genome as the simulated sample of that reference — real reads and simulated reads never
disagreed. Every one of them raised the threshold well above the floor -- 9 to 14 copies at 47x to 70x --
which is also checked.
*B. suis* SRR8550484 is the sample to watch across the two times this threshold changed: 422 of 1000 hashes
under a fixed 2, 775 once the threshold was first derived from the coverage, and 989 now.
[The k-mer copy threshold](K-mer-threshold) is that story.

**The sequence types agree, 16 of 16 and 4 of 4.** Every ST called from simulated reads is the one the
genome carries, read off the assembly by exact match: ST 1, 5, 6, 7, 8, 9, 12, 13, 14, 15, 17, 20, 26 and
53, with 7 and 14 called again at depth. The four real libraries keep giving ST 11, 1, 15 and 20.

**With room to spare, at 20x.** Across the 126 allele calls of the fourteen 20x *Brucella* samples, the
allele had 64 to 150 reads behind it (median 108) and beat its runner-up by 14 reads at the narrowest
(median 34). No locus came back as a tie or below the five-read minimum.

**What is not *Brucella* is not identified.** *Brucella anthropi* ATCC 49188, the nearest neighbour of the
genus (it was *Ochrobactrum anthropi* until the genus was redrawn), shares **4** of 1000 hashes with its
closest reference, and *E. coli* K-12 shares **0**. Both are reported as `ND`, far below the threshold of
150. There is a wide gap between a *Brucella* sample, which shares 993 or more at this depth, and anything
that is not one.

**How tight the genus is, in numbers.** The distance between the closest reference and the second closest
ranged from 0.00024 to 0.00352 on the simulated 20x samples. The narrowest was *B. canis* ATCC 23365,
which is 0.00024 away from *B. suis* bv. 4 strain 40 — ten hashes out of a thousand. The same measurement
taken within a species is one of the reasons there is no biovar in the reports: the two *B. melitensis*
biovar 1 references these runs measured against are 0.000756 apart, further than the biovar 2 and biovar 3
references are from each other ([more](How-it-works#no-biovar)).

**And the calls the distance did not make.** The simulated *B. canis* sample and the real one
(ERR2136546) are the two lines whose `species from` reads `arbitrated`: *B. suis* bv. 4 was inside the
margin both times, so MASH had not really chosen, and ST 20 did — it is *B. canis* in all 70 isolates of
the PubMLST isolate database that carry it, and never *B. suis*. The twenty other samples needed no
arbitration, and **not one of the twenty-two came back as a conflict** — the reference table and the ST to
species table never contradicted each other, which the check enforces.

*B. suis* 686 is the other instructive line: its `species from` is a plain `MASH`, because its ST 17 is
*B. suis* in only 18 of 20 isolates and does not clear the 95% bar. Nothing was lost — the distance had
named *B. suis* on its own, with no second species in reach — but it shows the threshold doing its job
rather than waving an impure ST through.

## What the validation caught

The first run of this validation, with the MLST database as vSNP shipped it, typed *B. abortus* 9-941 as a
**new** profile rather than ST 1. The reason: that copy of the database held a `gyrB_23` which is not a
variant of `gyrB_1` at all — it differs from it at 355 of 469 positions and occurs elsewhere in the
*B. abortus* genome — so reads from that other place voted for it and it won the locus. PubMLST has
deleted allele 23, and no profile uses that number any more. With
[PubMLST's data](https://github.com/duceppemo/BrucellaID/blob/main/scripts/make_mlst_data.py) the same
reads give ST 1. See the [Changelog](Changelog).

The depth problem travelled the other way: it was found on real libraries, because every simulated sample
here is 20x, and at 20x a fixed threshold of 2 is close enough to the right answer to hide the fault. The
150x and 250x samples, the hash floor of 950 and the threshold checks exist so that this validation would
have caught it — forcing the fixed threshold back on (`--min-kmer-copies 2`) now fails the run with four
problems on the two deep samples, which keep 902 and 728 of their 1000 hashes against the floor of 950.

**And it caught the second change to that threshold.** `check_results.py` restates the rule in its own
arithmetic rather than importing it from `brucellaid.mash`, deliberately, so that a drift in the shipped
formula shows up here as a failure instead of being matched silently. When the rule was re-derived from
one copy per 20x to one per 5x, this run duly failed: *the k-mer copy threshold is 30, expected about 8 at
150x*. The check was right to complain and the restatement was then updated by hand, which is the only
way it can work — a check that imports what it is checking cannot fail.

## What it does not establish

* **Messy reads.** The simulated reads have even coverage, 0.2% substitutions, no adapters, no
  contamination and no indels, and the four real libraries are clean runs of pure cultures. Nothing here
  is a metagenomic sample, a mixed culture or a degraded library.
* **Low depth.** The simulated default is 20x, and the shallowest real library is 47x. Run the script with
  `COVERAGE=5` to see where the ties begin.
* **Biovars.** Nothing at all: BrucellaID does not handle biovars, and this validation does not test for
  one. See [How it works](How-it-works#no-biovar).
* **Novel alleles.** Every genome here carries alleles the database holds.

The threshold those `k-mer copies` figures come from is derived in
[The k-mer copy threshold](K-mer-threshold).
