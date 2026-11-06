# How *B. papionis* got into the extended reference set

*B. papionis* is the thirtieth named species of the *Brucella* clade, described from baboons in 2014. For
a while it looked as though no reference set could ever carry it, because **nothing in INSDC is labelled
*B. papionis***: NCBI's assembly database returns no genome for taxon 1462582, ENA's returns none, and
the species had exactly one public read set, `SRR4038994`.

That turned out to be a statement about labels, not about sequence. The genome is there:
**`GCA_000163135.1`, which NCBI calls "*Brucella* sp. NVSL 07-0026"**. It is a *B. papionis* genome filed
under an unnamed label, and it is now the species' reference in the extended set, relabelled on the
evidence below.

## How it was found

A read set is what led to it. `SRR4038994` (BioSample `SAMN05417907`, strain **B07-0026**, study
`PRJNA251693`) was submitted as *B. papionis* and came back from the Phase 2 validation as
*B. pinnipedialis* at 0.0028 with an `unresolved` flag — the only sensible answer available when the
species has no reference. It was assembled here with SKESA 2.5.1 to fill that gap:

```
skesa --reads SRR4038994_1.fastq.gz,SRR4038994_2.fastq.gz --cores 16 --memory 48 \
      --contigs_out SRR4038994_skesa.fa      # 50 contigs, 3,243,917 bp, N50 116,444
```

Measured against the 2,022 assemblies the reference set was curated from, that assembly turned out to be
**0.000048 from `GCA_000163135.1`** — 998 of 1000 hashes, the same genome. NCBI's strain name for it,
*NVSL 07-0026*, is the read set's strain **B07-0026** with the laboratory's prefix on it. The species had
an assembly all along, and the curation had **dropped** it: the rule that removes an unnamed *Brucella*
sp. genome sitting within 0.004 of a named species treated it as a redundant twin of *B. microti*, 0.0020
away. In the classical core 0.0020 is a species-level gap, not a within-species one, so the rule was
wrong here. That is recorded in `scripts/extended_references_curation.tsv`, where the genome is now a
`RELABEL` instead of a `DROP-SP`.

## Why we are confident it is *B. papionis*

| Evidence | Result |
| --- | --- |
| The read set's submitter | *B. papionis* strain B07-0026 |
| ST from the reads of `SRR4038994` | ST 54 |
| ST read off the SKESA assembly, by exact allele match | `9-12-20-10-1-3-5-8-2` = ST 54 |
| ST read off `GCA_000163135.1`, by exact allele match | `9-12-20-10-1-3-5-8-2` = ST 54 |
| What PubMLST records for ST 54 | *B. papionis*, both of the isolates that carry it |
| `GCA_000163135.1` vs the SKESA assembly | 0.000048, 998/1000 hashes |
| `GCA_000163135.1` vs its nearest named species | 0.0019 *pinnipedialis*, 0.0020 *microti*, 0.0020 *ceti* |

The allele checks are exact string matches against the contigs and share no code with the read-based
caller. The distances are what a species of the classical core looks like: near-equidistant from several
named species and identical to none. A mislabelled member of one of those species would sit at 0.0002
from it instead.

With the genome in the set, `SRR4038994` reports *Brucella papionis* at 0.00083, 966 of 1000 hashes.

**What this does not show.** One strain is all the species has, so this cannot demonstrate that a
*different* *B. papionis* isolate would be identified correctly — the reference and the only test sample
are the same strain. It shows that the clade's one known *papionis* lineage is now covered instead of
being reported as its nearest neighbour.

BrucellaID does not use ST 54 as corroboration, and reports this call as plain `MASH` rather than
`MASH+ST (agree)`: an ST must rest on at least three PubMLST isolates to be allowed to speak, and ST 54
has two. That is the intended conservatism, not a failure.

## What is in this folder

`SRR4038994_B_papionis_B07-0026_skesa.fna.gz` — the SKESA assembly. It is **not** part of the reference
set and is not sketched; `GCA_000163135.1` is, because a genome with an accession is citable, is what
everyone else can download, and needs no special case in
`scripts/make_extended_references.py`. The assembly is kept because it is the independent evidence that
`GCA_000163135.1` is *B. papionis*: two assemblies of one strain, from different data by different
methods, agreeing on the nine alleles and sitting 0.000048 apart.

## Worth reporting

`GCA_000163135.1` is a named species sitting in NCBI under *Brucella* sp. It is listed in
[`../mislabelled_records.md`](../mislabelled_records.md) so it can be reported, which would give
*B. papionis* a labelled genome in the public archive for everyone rather than only in this reference set.
