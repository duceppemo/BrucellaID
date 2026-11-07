# PCR evidence for the mislabelled records

The records in [`../mislabelled_records.md`](../mislabelled_records.md) were found with a MASH distance and
a nine-locus MLST. This folder holds a third and fourth line of evidence, from methods that share nothing
with either: *in silico* PCR of published assays.

| Tool | Assays | Run on |
| --- | --- | --- |
| [brucetyper](https://github.com/duceppemo/brucetyper) | Bruce-ladder, Suis-ladder, AMOS, RB51 junction | every assembly |
| [insilicoPCR](https://github.com/duceppemo/insilicoPCR) 0.10.1 | the two real-time assays of Hoover et al. 2026 | every assembly **and** every raw read set |

The Hoover assay is [doi:10.1128/spectrum.02974-25](https://doi.org/10.1128/spectrum.02974-25) (CFIA Ottawa
Laboratory-Fallowfield): a genus singleplex for the classical *Brucella* species, and a multiplex with one
assay each for *B. abortus*, *B. melitensis* and *B. suis*. `hoover2026_qpcr.fasta` is that assay in
insilicoPCR's naming convention, the probes carrying the paper's sequences with the dye and quencher
labels stripped. The amplification verdict is insilicoPCR's consolidated report. Primers were required to
match **exactly** (`-m 0`), which is what reproduces the assay's designed specificity — see the caveat
below.

## What was run

* **46 NCBI assemblies**, downloaded from the NCBI FTP for this analysis: every assembly in the dossier,
  plus one genome per species of the whole clade as a control panel, type strains where one exists. The
  downloads were byte-identical to the copies the curation was done on.
* **48 public read sets**, downloaded with their ENA md5 checksums verified and assembled with
  SKESA 2.5.1: the 21 mislabelled read sets, the 4 former *Ochrobactrum* candidates, the *B. papionis*
  read set, and all 22 read sets that appear to hold no *Brucella*. 47 produced an assembly;
  `SRR24420454` produced nothing, which is a result in itself.

`evidence_assemblies.tsv` and `evidence_reads.tsv` carry one row per subject with every assay verdict, and
`evidence_fastq_diagnostics.tsv` one row per assay that the read-mode mini-assembly missed, with why.

## What it establishes

**The control panel validates the assay and bounds it.** Across one genome per species, the genus
singleplex is positive on every classical species and negative on all 17 former *Ochrobactrum* species,
which is what the paper's "classical *Brucella* species" target implies. The three species assays are
perfectly specific: no cross-reaction anywhere in the clade.

**The mislabelled read sets are good *Brucella* genomes with the wrong name on them.** 20 of the 21 are
genus-positive, and the exception is the one whose evidence says *B. anthropi*, where a negative is the
correct answer. 19 of 21 assemble to 3.24–3.35 Mb, textbook *Brucella*, so these are not junk libraries
that confused a distance — they are real genomes carrying a wrong species. Of the 17 whose verdict species
the multiplex covers, **14 fired the matching species assay**; two of the three that did not are the two
worst assemblies in the set.

**The swap pair is confirmed from both directions**, which is the strongest single result here. The genome
labelled *B. inopinata* BO1 behaves like *B. suis*, the genome labelled *B. suis* behaves like
*B. inopinata*, and the genuine BO1 assembly behaves differently from the one claiming to be it.

**20 of the 22 "no *Brucella*" read sets are PCR-negative for a species on both substrates**,
independently confirming what the shared-hash counts said. The two exceptions are **mixed libraries**
rather than mistakes, and both are genus-positive, so each does hold some classical *Brucella*:

* `SRR952243` is negative on its assembly but genus- and *B. abortus*-positive on its reads, so a minority
  *B. abortus* population sits under a dominant other organism — the library assembles cleanly into
  5.4 Mb, far too large for *Brucella*. This is the one case where reads see more than the assembly, and
  0.10.1 confirms it.
* `SRR35797377` is genus-positive on both substrates and *B. melitensis*-positive on an assembly that is
  10,733 contigs of mixture; in the reads that *melitensis* amplicon is in a single read, which the
  threshold rejects. So: a trace of *B. melitensis* in a library that holds a little classical *Brucella*
  and a great deal of something else.

Neither is a counter-example to the shared-hash verdict — 26 and 146 of 1000 hashes is what a minority
population looks like — but "mixed library" describes them better than "not *Brucella*".

**For the former *Ochrobactrum* records this evidence is silent.** All four candidates are negative on
every assay on both substrates, and Bruce-ladder is equally silent, because neither assay was designed for
that half of the clade. Those records stand or fall on the distances alone.

## What the non-*Brucella* read sets of PRJEB4782 actually hold

The ten runs of that study which share one or two of 1000 hashes with the clade were identified, because
"not *Brucella*" is less useful to a submitter than a name. `prjeb4782_unidentified_libraries.tsv` has the
per-run detail. Four methods agree: mashID against proGenomes v4 and against a dereplicated RefSeq
bacterial sketch, fastANI against every *Camelimonas* genome in NCBI, and 16S rRNA extracted with barrnap
and compared with both the NCBI type-strain database and nt.

| Runs | Identification | Evidence |
| --- | --- | --- |
| `ERR471337`, `ERR554820`, `ERR775236` | ***Camelimonas fluminis*** | 98.35–98.39% ANI to the type strain, 16S 100% |
| `ERR473743`, `ERR775237` | ***Camelimonas***, unnamed species | 95.5% ANI to *Camelimonas* sp. ID_303_24, only 89.6% to the *C. fluminis* type strain |
| `ERR775244`–`ERR775248` (5) | **Alcaligenaceae**, unnamed genus | full-length 16S 94% to the *Advenella kashmirensis* type strain; **no hit in either genome database** |

*Camelimonas* explains the deposit rather than merely contradicting it: it sits in **Brucellaceae**, the
same family as *Brucella*, it is a slow-growing oxidase-positive alphaproteobacterium, and the genus
contains a *Camelimonas abortus*. An isolate typed biochemically in a brucellosis workup could reasonably
have been recorded as *Brucella*. The reference the unnamed-species pair matches is itself a clinical
isolate from human blood.

The five 1.9 Mb libraries are one organism — 99.3–100% ANI to each other, all 52.7% GC, assembled at
250–270×, so the genome is genuinely that small rather than half-covered. At 94% 16S identity to its
nearest named relatives it falls below genus-level similarity, and it has no match at all in 32,887
proGenomes representatives or in the RefSeq sketch. A reduced-genome member of a family whose relatives
run 3.5–4.5 Mb, with no sequenced representative, is worth flagging to whoever holds the isolate.

## Two caveats worth passing on

**Read mode used to lose targets silently; it no longer does.** The first pass of this analysis ran
insilicoPCR 0.9.0, where the reads disagreed with the whole-genome assembly on ten of the 48 read sets and
nine of those were the assembly seeing *more*. The cause was the Tadpole mini-assembly that read mode
searches against: it omitted the target region, or broke inside the amplicon so the two primers landed on
different contigs, while the reads carrying the oligos sat in the baited set all along. That was reported
as [insilicoPCR#7](https://github.com/duceppemo/insilicoPCR/issues/7) and fixed in 0.10.0, which calls
such targets directly from the reads and writes a diagnostics file explaining every miss.

Rerun on 0.10.1, **reads and assembly now agree on 44 of the 48 read sets, up from 38**. Ten read sets
gained an assay and none lost one, nine of them recovering exactly what their whole-genome assembly had
reported all along. The assembly-mode results are byte-for-byte identical between the two versions, so
nothing above rests on the change.

Two read sets gained more than a recovered negative, and both had poor assemblies — one of 2,165 contigs,
one of 4,251:

* `SRR867710`, the third library of the Broad BioSample, now supports *B. suis* where its own assembly
  gave only a genus signal, so **all three libraries of that BioSample now agree**.
* `ERR2938687` reports *two* species assays, *B. suis* in 54% of informative reads and *B. melitensis* in
  33%, just over the threshold. It is the weakest of the eleven PRJEB30030 findings on distance, and it
  has the worst assembly. **A mixture fits all three observations better than a simple mislabel**, and it
  should be reported that way.

The threshold earns its keep in the other direction too. Ten assays are flagged `WEAK_READ_SUPPORT` —
present in the reads but below the 2-read, 30% minimum — and they are the ones a looser rule would have
turned into false positives: 2% here, 15% there, and two single-read traces in libraries that hold no
*Brucella* to speak of.

**The *B. suis* assay's specificity rests on exact matching.** Tolerating two mismatches
(`-m 2 --3p-penalty 3`) makes it fire on nearly every classical *Brucella* genome, including *abortus*,
*melitensis*, *ceti*, *microti*, *neotomae*, *ovis*, *pinnipedialis* and *vulpis*. This does not affect
anything above, because the exact-match run is the one that reproduces the designed assay and a real
reaction enforces that stringency through its annealing temperature. It is noted because an *in silico*
screen run with a loose mismatch budget would reach the wrong conclusion, and because the assay's authors
may want to know.
