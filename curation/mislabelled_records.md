# Mislabelled public records found while curating the extended reference set

Records in NCBI (and, where noted, in the [mashID *Brucella* pack](https://doi.org/10.6084/m9.figshare.33970168))
whose sequence contradicts the species on the label. They are kept here so they can be reported to the
submitters and to NCBI with the evidence attached, rather than only silently excluded from our reference set.

**Nothing here is a claim about anyone's laboratory work.** A label can be wrong because of a strain
mix-up, a swapped submission, a taxonomy that moved under it, or a species boundary that genome data drew
differently from the phenotype it was named on. The evidence says what the sequence is, not how the label
came to be.

**How to read the evidence.** Every distance is a MASH distance (k = 21). Reference-to-reference distances
come from an all-against-all of 454 dereplicated genomes sketched at 10,000 hashes; read-to-reference
distances come from BrucellaID's packaged sketch at 1,000 hashes. A distance is only meaningful next to
another distance, so each row carries the comparison that makes it mean something: how far the genome is
from the species it claims against how far it is from the species the sequence puts it in. For context,
two genomes of the same classical *Brucella* species are typically 0.0001–0.002 apart, while the former
*Ochrobactrum* species are far more diverse — within *B. anthropi* the spread reaches 0.21.

**Status** is `to report` until it has been sent, then `reported` with the date, then `resolved` or
`disputed`. The reporting log at the end of this file is the record of what was sent where.

## Assemblies: the label is contradicted by the sequence

Each of these rests on two independent lines of evidence, or on type-strain geometry with a margin wide
enough that one line is decisive. All five are also in the figshare pack, so the pack carries the same
labels.

| Accession | Label in NCBI | The sequence is | Evidence | Status |
| --- | --- | --- | --- | --- |
| `GCA_036881765.1` | *B. abortus* | *B. melitensis* | 0.0007 from *melitensis*, 0.0016 from *abortus*. Carries ST 8, which PubMLST records as *melitensis* in 98.3% of 862 isolates. Two independent lines. | to report |
| `GCA_000662015.2` | *B. inopinata* BO1, **flagged as type material** | a *suis*/*canis* cluster genome | 0.00114 from *B. suis* `GCA_000741995.1` and 0.00118 from *B. canis* `GCA_000530495.1`, but **0.0237 from the other assembly of BO1**. It cannot be the same strain as the genome it claims to be. | to report |
| `GCA_000661995.2` | *B. suis* bv. 4 str. 40 | *B. inopinata* (BO1) | 0.00044 from the BO1 type assembly `GCA_000182725.1` and 0.00203 from a second *inopinata*, against 0.022 from every genome labelled *suis*. Its three nearest neighbours are all *inopinata*. | to report |
| `GCA_049065255.1` | *B. ciceri* (SA289) | *B. intermedia* | 0.00004 from `GCA_001917355.1`, which NCBI labels *intermedia*; 0.0066 from the nearest other *ciceri* and from the *ciceri* type strain DSM 22292. The SA148/IR068 cluster it sits in is labelled *intermedia* in 12 of 14 NCBI records. | to report |
| `GCA_059995225.1` | *B. cytisi* (WTB2) | *B. anthropi* | 0.0066 from *anthropi* against 0.0395 from the nearest other *cytisi*. All three deposits of the *cytisi* type strain (ESC1, LMG 22713, DSM 19778) agree with each other to 0.0012–0.0026 and sit 0.0381–0.0394 away, so the type cluster is tight and this genome is not in it. | to report |

`GCA_000662015.2` and `GCA_000661995.2` are **a swapped pair**: each one's sequence matches the other's
label. They are consecutive accessions from one submission, which is where the explanation probably lies.
The *inopinata* half is the one worth raising first, because it is flagged as **type material** — the
reference a future study would anchor on.

## Assemblies: candidates, one line of evidence only

These are flagged by geometry alone and need type-strain anchoring per species before anything is claimed.
The ex-*Ochrobactrum* species are diverse enough that "nearer to another species" is suggestive, not
conclusive. The ratio column is how much nearer the other species is than the labelled one.

| Accession | Label in NCBI | Nearer species | Distance | To its own label | Ratio | Note |
| --- | --- | --- | --- | --- | --- | --- |
| `GCA_001586735.1` | *B. anthropi* | *B. pseudogrignonensis* | 0.0129 | 0.1827 | 14× | |
| `GCA_000559065.1` | *B. rhizosphaerae* | *B. pituitosa* | 0.0141 | 0.1158 | 8× | only 4 *rhizosphaerae* genomes exist to compare with |
| `GCA_024127925.1` | *B. intermedia* | *B. pseudintermedia* | 0.0113 | 0.0835 | 7× | the two species are sisters |
| `GCA_002025625.1` | *B. pituitosa* | *B. pseudogrignonensis* | 0.0260 | 0.1240 | 5× | |
| `GCA_001575075.1` | *B. anthropi* | *B. tritici* | 0.0192 | 0.0717 | 3.7× | |
| `GCA_979502475.1` | *B. lupini* | *B. anthropi* | 0.0051 | 0.0111 | 2.2× | weakest of the set |
| `GCA_000292145.2` | *B. abortus* bv. 1 (BCB027) | an unnamed lineage | 0.0010 | 0.0016 | 1.6× | known atypical genome; margin is narrow |
| `GCA_000292205.2` | *B. abortus* bv. 1 (BCB013) | *B. suis* | 0.0011 | 0.0022 | 2× | ST 120 is a PubMLST singleton recorded as *abortus*, so the ST does **not** corroborate |
| `GCA_028621035.1` and `GCA_042679345.1` | both *B. gallinifaecis* | — | — | — | — | the only two *gallinifaecis* genomes are **0.167 apart**, so at most one of them can be that species. Needs the type strain to say which. |

### Checked and cleared

Kept here so they are not re-flagged by the next automated pass.

* `GCA_000182725.1` (*B. inopinata* BO1, type material) — **correctly labelled.** An automated rule flags
  it because its nearest neighbour at 0.00044 wears a *suis* label, but that neighbour is
  `GCA_000661995.2`, the mislabelled genome above. Its next neighbour is a genuine *inopinata* at 0.00167.
* `GCA_000741995.1` (*B. suis*) and `GCA_000530495.1` (*B. canis*) — 0.00019 apart, but that is real
  biology: *B. canis* sits inside the *B. suis* radius. Not a labelling error, and the reason BrucellaID
  uses the sequence type to separate those two species.

## Read sets: the submitted species is contradicted by two independent methods

From the Phase 2 validation over public ENA read sets. For each of these the MASH identification and the
nine-locus MLST — which share no code and no reference data — both point away from the submitted label and
at the same species. Distances are to the closest genome of the extended set. The `PubMLST species` column
gives the species that database records for the sequence type, with how consistently its isolates carry it.

Grouped by study, because several come from a single submission and are better raised as one item with the
submitter than as separate reports.

### PRJEB30030 — seven read sets, one submitter

National Institute of Health Dr. Ricardo Jorge. Seven distinct BioSamples, each of which the sequence puts
in a different species from its label.

| Run | BioSample | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `ERR2938660` | `SAMEA5145035` | *B. ovis* | *B. suis* | 0.00081 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2938669` | `SAMEA5145008` | *B. melitensis* | *B. suis* | 0.00181 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2938678` | `SAMEA5145000` | *B. melitensis* | *B. abortus* | 0.00111 | 5 | *B. abortus*, 98.6% of 71 isolates |
| `ERR2938687` | `SAMEA5145007` | *B. melitensis* | *B. suis* | 0.00591 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2993134` | `SAMEA5176151` | *B. melitensis* | *B. pinnipedialis* | 0.00026 | 25 | *B. pinnipedialis*, 95.0% of 80 isolates |
| `ERR2993152` | `SAMEA5176142` | *B. melitensis* | *B. abortus* | 0.00014 | 1 | *B. abortus*, 99.4% of 466 isolates |
| `ERR2993155` | `SAMEA5176143` | *B. melitensis* | *B. suis* | 0.00014 | 14 | *B. suis*, 96.4% of 140 isolates |

`ERR2993134` is the most clear-cut: a read set submitted as *B. melitensis* that is 0.00026 from
*B. pinnipedialis* and carries ST 25, a marine-mammal sequence type. `ERR2938687` is the weakest, because
its distance of 0.0059 is large enough that the species rests mainly on the ST.

### PRJNA173595 and PRJNA173601 — Broad Institute reference strains

| BioSample | Runs | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `SAMN00811170` | `SRR642185`, `SRR651210`, `SRR867710` | *B. abortus* 63/183 | *B. suis* | 0.0022–0.0023 | 14, 16 | *B. suis* (ST 16 is 100% of 15 isolates) |
| `SAMN00811205` | `SRR867711` | *B. melitensis* 63/227 | *B. suis* | 0.00205 | 14 | *B. suis*, 96.4% of 140 isolates |

The first row is **one BioSample sequenced in three runs**, not three records: all three runs agree, which
is a consistency check rather than three independent findings. These two carry strain designations
(63/183, 63/227), so a wrong label here would propagate further than an ordinary isolate's would.

### PRJEB4782 — Wellcome Sanger Institute

| Run | BioSample | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `ERR471314` | `SAMEA2266858` | *B. ceti* | *B. abortus* | 0.00017 | 5 | *B. abortus*, 98.6% of 71 isolates |
| `ERR471315` | `SAMEA2266859` | *B. ceti* | *B. abortus* | 0.00017 | 5 | *B. abortus*, 98.6% of 71 isolates |

### Single submissions

| Run | BioSample | Study | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `SRR31560207` | `SAMN15685105` | PRJNA649281 | *B. abortus* | *B. melitensis* | 0.00056 | 7 | *B. melitensis*, 98.6% of 72 isolates |
| `SRR11526759` | `SAMN14309819` | PRJNA610713 | *B. abortus* | *B. melitensis* | 0.00098 | 8 | *B. melitensis*, 98.3% of 862 isolates |

`SRR31357708` (`SAMN44729562`, PRJNA1183633), submitted as *B. cytisi*, belongs with these although no ST
is available: the former *Ochrobactrum* species do not carry the classical MLST alleles, and six of the
nine loci got no read at all. The distance evidence is decisive on its own — 0.0076 from the nearest
*B. anthropi* against 0.0381–0.0394 from the three deposits of the *cytisi* type strain, which rank 89th,
90th and 91st of the 428 references. Every one of its top eight hits is *anthropi*. Status: to report.

## Read sets: candidates

| Run | BioSample | Submitted as | Closest | Distance | Why it is only a candidate |
| --- | --- | --- | --- | --- | --- |
| `SRR6003475` | `SAMN07581394` | *B. thiophenivorans* | *B. pituitosa* | 0.0029 | sister species with one reference genome each, no ST available, no type-strain anchor measured yet |
| `SRR10948139` | `SAMN13812480` | *B. anthropi* | *B. intermedia* | 0.0087 | both species are diverse enough that 0.0087 is not conclusive, and no ST is available |

## Read sets where the label looks right and our reference set is what is missing

Kept here so they are not mistaken for mislabels. In each case the submitter is probably correct and
BrucellaID could not place the sample because the extended set has no suitable reference.

* `SRR4038994` (`SAMN05417907`, PRJNA251693), submitted as *B. papionis*, reported as
  *B. pinnipedialis* at 0.0028 with an `unresolved` flag. The label is almost certainly right: the run
  carries **ST 54, which PubMLST records as *B. papionis* in both of its isolates**, and the reason MASH
  cannot confirm it is that *B. papionis* has no assembly anywhere in INSDC, so the extended set has no
  genome of it to match. This is the one named species of the clade the set cannot cover, and this read
  set is the only public sequence data for it.
* `SRR12605656`, submitted as *B. haematophila*, reported as *Brucella* sp. at 0.0201. Only three
  *haematophila* genomes exist to compare with, and one of them is itself a candidate mislabel, so the
  reference side is too thin to call this either way.
* `SRR31481153`, submitted as *B. inopinata*, reported as *Brucella* sp. at 0.00053 with a `new` sequence
  type. It sits next to an unnamed lineage representative rather than to an *inopinata* genome, which may
  say more about how few *inopinata* genomes exist than about the label.

## Reproducing the evidence

```bash
# reference-to-reference distances (sketched at 10,000 hashes)
mash dist curation.msh curation.msh > allvsall.tsv

# a read set against the extended reference set, keeping every per-reference distance
brucellaid -r1 RUN_1.fastq.gz -r2 RUN_2.fastq.gz -o out --reference-set extended
#   out/brucellaid_mash.tsv holds the distance to all 428 references, ranked
#   out/brucellaid_mlst.tsv holds every allele call with the margin it had
```

The curation decisions that excluded these genomes from the shipped set are in
[`scripts/extended_references_curation.tsv`](../scripts/extended_references_curation.tsv).

## Reporting log

| Record | Reported to | Date | Outcome |
| --- | --- | --- | --- |
| — | — | — | nothing sent yet |

## Still to come

Phase 2 is still running at the time of writing: the read-set tables above cover the first 2,929 of the
4,044 runs and will grow. Re-run the extraction over the finished `results_default.tsv` and
`results_extended.tsv` before anything is sent, so that one report covers the whole sweep, and group any
new runs by study the way the tables above do -- six of the first seven findings came from one submission,
and raising them together is both fairer to the submitter and harder to dismiss.
