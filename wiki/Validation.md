# Validation

Sixteen genomes whose answer is known, reads simulated from each one, and the reports compared with what
the genome really is. Reproduce it with:

```bash
bash validation/run_validation.sh      # about 25 minutes and 1 GB, all of it in validation/data
```

The script downloads the genomes listed in
[`validation/genomes.tsv`](https://github.com/duceppemo/BrucellaID/blob/main/validation/genomes.tsv) from
NCBI, simulates 20x of 150 bp paired-end reads from each (random fragments, 0.2% substitutions, one seed
per sample), runs BrucellaID over the folder and checks the result.

**The sequence type is not listed anywhere.** `check_results.py` reads it off each genome by looking for an
exact match of every allele of the package, on both strands, and compares that with what BrucellaID called
from the reads. The two answers share no code: one searches an assembly for a string, the other counts
k-mer votes in reads.

It also checks that no sample comes back as `MASH (conflict)`: these genomes are what they say they are, so
the distance and the sequence type must not name species that are nowhere near each other.

## What came back

```
sample                   species              biovar of the genome  reference                      hashes     ST called  ST in the genome  species from          result
-----------------------  -------------------  --------------------  -----------------------------  ---------  ---------  ----------------  --------------------  ------
abortus_9-941            Brucella abortus     1                     NC_006932-NC_006933.fasta      999/1000   1          1                 MASH+ST (agree)       ok
abortus_BER              Brucella abortus     3                     NZ_CP007682-NZ_CP007683.fasta  998/1000   6          6                 MASH+ST (agree)       ok
melitensis_16M           Brucella melitensis  1                     NC_003317-NC_003318.fasta      995/1000   7          7                 MASH+ST (agree)       ok
melitensis_ATCC23457     Brucella melitensis  2                     NC_012441-NC_012442.fasta      999/1000   8          8                 MASH+ST (agree)       ok
melitensis_Ether         Brucella melitensis  3                     NZ_CP007760-NZ_CP007761.fasta  1000/1000  9          9                 MASH+ST (agree)       ok
melitensis_BwIM_SOM_36b  Brucella melitensis  1                     NZ_CP018508-NZ_CP018509.fasta  998/1000   12         12                MASH+ST (agree)       ok
suis_1330                Brucella suis        1                     NC_017251-NC_017250.fasta      995/1000   14         14                MASH+ST (agree)       ok
suis_ATCC23445           Brucella suis        2                     NC_010169-NC_010167.fasta      998/1000   15         15                MASH+ST (agree)       ok
suis_686                 Brucella suis        3                     NZ_CP007719-NZ_CP007718.fasta  999/1000   17         17                MASH                  ok
canis_ATCC23365          Brucella canis       -                     NC_010103-NC_010104.fasta      1000/1000  20         20                MASH+ST (arbitrated)  ok
ovis_ATCC25840           Brucella ovis        -                     NC_009505-NC_009504.fasta      996/1000   13         13                MASH+ST (agree)       ok
ceti_TE10759-12          Brucella ceti        -                     NC_022905-NC_022906.fasta      995/1000   26         26                MASH+ST (agree)       ok
microti_CCM4915          Brucella microti     -                     NC_013119_NC_013118.fasta      998/1000   53         53                MASH+ST (agree)       ok
abortus_2308             Brucella abortus     -                     NC_006932-NC_006933.fasta      992/1000   5          5                 MASH+ST (agree)       ok
anthropi_ATCC49188       ND                   -                     B-REF-BS4-40.fasta             4/1000     ND         ND (no locus)                           ok
coli_MG1655              ND                   -                     B-REF-BS4-40.fasta             0/1000     ND         ND (no locus)                           ok

16 sample(s), all as expected.
```

## What that establishes

**The identification names the right genome.** Each of the 13 genomes that are in the sketch came back as
itself, with 995 to 1000 of its 1000 hashes shared, and with the species the reference table gives it. The
"biovar of the genome" column is what `genomes.tsv` says each genome is, not something BrucellaID
reported: it reports no biovar (see [How it works](How-it-works#no-biovar)), and the column is there so
that a reader can see which biovars the check covered.

**And the right species for a genome that is not in the sketch.** *B. abortus* 2308 is not one of the 16.
It was identified as *Brucella abortus* through 9-941, sharing 992 hashes.

**The sequence types agree, 14 out of 14.** Every ST that BrucellaID called from the reads is the one the
genome carries, read off the assembly by exact match: ST 1, 5, 6, 7, 8, 9, 12, 13, 14, 15, 17, 20, 26 and
53. That is 126 allele calls, all correct.

**With room to spare, at 20x.** Across those 126 calls, the allele had 64 to 150 reads behind it (median
108) and beat its runner-up by 14 reads at the narrowest (median 34), a ratio of 1.19 at the narrowest and
1.45 at the median. No locus came back as a tie or below the five-read minimum.

**What is not *Brucella* is not identified.** *Brucella anthropi* ATCC 49188, the nearest neighbour of the
genus (it was *Ochrobactrum anthropi* until the genus was redrawn), shares **4** of 1000 hashes with its
closest reference, and *E. coli* K-12 shares **0**. Both are reported as `ND`, far below the threshold of
150. There is a wide gap between a *Brucella* sample, which shares 992 or more at this depth, and anything
that is not one.

**How tight the genus is, in numbers.** The distance between the closest reference and the second closest
ranged from 0.00024 to 0.00354. The narrowest was *B. canis* ATCC 23365, which is 0.00024 away from
*B. suis* bv. 4 strain 40 — ten hashes out of a thousand. The same measurement taken within a species is
why there is no biovar in the reports: the two *B. melitensis* biovar 1 references of this very sketch are
0.000756 apart, further than the biovar 2 and biovar 3 references are from each other.

**And the one call the distance did not make.** That *B. canis* sample is the only one of the sixteen whose
`species from` reads `arbitrated`: *B. suis* bv. 4 was inside the margin, so MASH had not really chosen, and
ST 20 did — it is *B. canis* in all 70 isolates of the PubMLST isolate database that carry it, and never
*B. suis*. The answer is the same as before this rule existed, which is the point: it was right by ten
hashes, and it is now right for a reason. The other fifteen samples needed no arbitration, and **not one of
them came back as a conflict** — on sixteen genomes whose species is known, the reference table and the ST
to species table never contradicted each other, which the check now enforces.

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

## What it does not establish

* **Real reads.** These are simulated: even coverage, 0.2% substitutions, no adapters, no contamination and
  no indels. The [Tutorial](Tutorial) runs four real public read sets, where the numbers are less tidy — one
  of them shares 422 hashes of 1000 and is still typed.
* **Low depth.** Everything here is 20x. Run the script with `COVERAGE=5` to see where the ties begin.
* **Biovars.** Nothing at all: BrucellaID does not report a biovar, and this validation does not test for
  one. See [How it works](How-it-works#no-biovar) for the three reasons.
* **Novel alleles.** Every genome here carries alleles the database holds.
