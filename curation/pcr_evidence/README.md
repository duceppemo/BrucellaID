# PCR evidence for the mislabelled records

The records in [`../mislabelled_records.md`](../mislabelled_records.md) were found with a MASH distance and
a nine-locus MLST. This folder holds a third and fourth line of evidence, from methods that share nothing
with either: *in silico* PCR of published assays.

| Tool | Assays | Run on |
| --- | --- | --- |
| [brucetyper](https://github.com/duceppemo/brucetyper) | Bruce-ladder, Suis-ladder, AMOS, RB51 junction | every assembly |
| [insilicoPCR](https://github.com/duceppemo/insilicoPCR) 0.9.0 | the two real-time assays of Hoover et al. 2026 | every assembly **and** every raw read set |

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

`evidence_assemblies.tsv` and `evidence_reads.tsv` carry one row per subject with every assay verdict.

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

**20 of the 22 "no *Brucella*" read sets are PCR-negative on both substrates**, independently confirming
what the shared-hash counts said. The two exceptions are more interesting than a clean negative would
have been, and both are **mixed libraries** rather than mistakes: `SRR952243` is negative on its assembly
but genus- and *abortus*-positive on its reads, so it holds a minority *B. abortus* population under a
dominant other organism — it assembles cleanly into 5.4 Mb, far too large for *Brucella*. `SRR35797377` is
*melitensis*-positive on an assembly that is 10,733 contigs of mixture, with no genus signal, so it is
equivocal and is reported that way.

**For the former *Ochrobactrum* records this evidence is silent.** All four candidates are negative on
every assay on both substrates, and Bruce-ladder is equally silent, because neither assay was designed for
that half of the clade. Those records stand or fall on the distances alone.

## Two caveats worth passing on

**A negative on raw reads is weak evidence of absence, and the reason is the mini-assembly step.** The
assembly and the reads agree on 38 of the 48 read sets. Of the ten that differ, nine have *more* signal
from the whole-genome assembly than from the reads.

In FASTQ mode insilicoPCR baits reads matching the oligos, recruits the reads overlapping those, assembles
the recruited set with Tadpole, and searches that mini-assembly. Tracing one of the nine through its
retained intermediates shows the loss is at the assembly step, not at baiting and not at the search: the
reads covering the missing target are present in both the baited and the recruited set — more of them than
for a target that did report in the same sample — and Tadpole's default thresholds simply build no contig
over them, so the target is absent from what the search sees. Re-running Tadpole on the identical recruited
reads with a lower seed count, coverage floor and minimum contig length recovers the missing target and
takes the mini-assembly from 5 contigs to 31.

So these negatives are an artefact of the default assembly settings rather than evidence that the target
is absent, and the remedy is to lower those thresholds or to read the verdict off the whole-genome
assembly. Worth passing to the tool's author.

The one case that runs the other way, `SRR952243`, is the one that matters biologically: a minority
population present in the reads and absent from its whole-genome assembly.

**The *B. suis* assay's specificity rests on exact matching.** Tolerating two mismatches
(`-m 2 --3p-penalty 3`) makes it fire on nearly every classical *Brucella* genome, including *abortus*,
*melitensis*, *ceti*, *microti*, *neotomae*, *ovis*, *pinnipedialis* and *vulpis*. This does not affect
anything above, because the exact-match run is the one that reproduces the designed assay and a real
reaction enforces that stringency through its annealing temperature. It is noted because an *in silico*
screen run with a loose mismatch budget would reach the wrong conclusion, and because the assay's authors
may want to know.
