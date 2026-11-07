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
enough that one line is decisive. All six are also in the figshare pack, so the pack carries the same
labels.

| Accession | Label in NCBI | The sequence is | Evidence | Status |
| --- | --- | --- | --- | --- |
| `GCA_036881765.1` | *B. abortus* | *B. melitensis* | 0.0007 from *melitensis*, 0.0016 from *abortus*. Carries ST 8, which PubMLST records as *melitensis* in 98.3% of 862 isolates. Two independent lines. | to report |
| `GCA_000662015.2` | *B. inopinata* BO1, **flagged as type material** | a *suis*/*canis* cluster genome | 0.00114 from *B. suis* `GCA_000741995.1` and 0.00118 from *B. canis* `GCA_000530495.1`, but **0.0237 from the other assembly of BO1**. It cannot be the same strain as the genome it claims to be. | to report |
| `GCA_000661995.2` | *B. suis* bv. 4 str. 40 | *B. inopinata* (BO1) | 0.00044 from the BO1 type assembly `GCA_000182725.1` and 0.00203 from a second *inopinata*, against 0.022 from every genome labelled *suis*. Its three nearest neighbours are all *inopinata*. | to report |
| `GCA_049065255.1` | *B. ciceri* (SA289) | *B. intermedia* | 0.00004 from `GCA_001917355.1`, which NCBI labels *intermedia*; 0.0066 from the nearest other *ciceri* and from the *ciceri* type strain DSM 22292. The SA148/IR068 cluster it sits in is labelled *intermedia* in 12 of 14 NCBI records. | to report |
| `GCA_000163135.1` | *Brucella* sp. NVSL 07-0026 | *B. papionis* | Carries ST 54, which PubMLST records as *B. papionis*; the same nine alleles, by exact match, as an independent assembly of `SRR4038994` (strain B07-0026, submitted as *B. papionis*) sitting 0.000048 away. 0.0019–0.0020 from *pinnipedialis*, *microti* and *ceti*, so a species of its own. **This is the only *B. papionis* genome in INSDC and nothing there is labelled with that species.** | to report |
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

From the Phase 2 validation, which measured **4,042 of the 4,044 usable public *Brucella* read sets in ENA**
against both reference sets. For each of these the MASH identification and the nine-locus MLST — which share no
code and no reference data — both point away from the submitted label and at the same species. Distances
are to the closest genome of the extended set. The `PubMLST species` column gives the species that database
records for the sequence type, with how consistently its isolates carry it.

Grouped by study, because most come from a handful of submissions and are far better raised as one item
with each submitter than as twenty separate reports.

### PRJEB30030 — eleven read sets, one submitter

National Institute of Health Dr. Ricardo Jorge. Eleven distinct BioSamples, each of which the sequence
puts in a different species from its label. This is the largest single cluster in the whole sweep and the
first one to raise.

| Run | BioSample | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `ERR2938660` | `SAMEA5145035` | *B. ovis* | *B. suis* | 0.00081 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2938669` | `SAMEA5145008` | *B. melitensis* | *B. suis* | 0.00181 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2938678` | `SAMEA5145000` | *B. melitensis* | *B. abortus* | 0.00111 | 5 | *B. abortus*, 98.6% of 71 isolates |
| `ERR2938687` | `SAMEA5145007` | *B. melitensis* | *B. suis* | 0.00591 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2938695` | `SAMEA5145011` | *B. melitensis* | *B. suis* | 0.00244 | 15 | *B. suis*, 96.1% of 103 isolates |
| `ERR2993134` | `SAMEA5176151` | *B. melitensis* | *B. pinnipedialis* | 0.00026 | 25 | *B. pinnipedialis*, 95.0% of 80 isolates |
| `ERR2993138` | `SAMEA5176141` | *B. melitensis* | *B. suis* | 0.00014 | 19 | *B. suis*, 100% of 6 isolates |
| `ERR2993140` | `SAMEA5176147` | *B. melitensis* | *B. neotomae* | 0.00012 | 22 | *B. neotomae*, 100% of 9 isolates |
| `ERR2993145` | `SAMEA5176140` | *B. melitensis* | *B. suis* | 0.00007 | — | no ST; **both reference sets agree independently** |
| `ERR2993152` | `SAMEA5176142` | *B. melitensis* | *B. abortus* | 0.00014 | 1 | *B. abortus*, 99.4% of 466 isolates |
| `ERR2993155` | `SAMEA5176143` | *B. melitensis* | *B. suis* | 0.00014 | 14 | *B. suis*, 96.4% of 140 isolates |

Several of these are 0.0001–0.0003 from the species the sequence says they are, which is as close as two
genomes of one species ever get. `ERR2993134` is the most striking: submitted as *B. melitensis*, 0.00026
from *B. pinnipedialis* and carrying ST 25, a marine-mammal sequence type. `ERR2938687` is the weakest of
the eleven, because its distance of 0.0059 is large enough that the species rests mainly on the ST.

### PRJEB4782 — three read sets, Wellcome Sanger Institute

| Run | BioSample | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `ERR471314` | `SAMEA2266858` | *B. ceti* | *B. abortus* | 0.00017 | 5 | *B. abortus*, 98.6% of 71 isolates |
| `ERR471315` | `SAMEA2266859` | *B. ceti* | *B. abortus* | 0.00017 | 5 | *B. abortus*, 98.6% of 71 isolates |
| `ERR473742` | `SAMEA2266954` | *B. abortus* | *B. neotomae* | 0.00034 | 22 | *B. neotomae*, 100% of 9 isolates |

The same study also contributed ten read sets that contain no *Brucella* sequence at all — see the section
below. Both findings are worth raising together.

### PRJNA173595 and PRJNA173601 — Broad Institute reference strains

| BioSample | Runs | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- |
| `SAMN00811170` | `SRR642185`, `SRR651210`, `SRR867710` | *B. abortus* 63/183 | *B. suis* | 0.0022–0.0023 | 14, 16 | *B. suis* (ST 16 is 100% of 15 isolates) |
| `SAMN00811205` | `SRR867711` | *B. melitensis* 63/227 | *B. suis* | 0.00205 | 14 | *B. suis*, 96.4% of 140 isolates |

The first row is **one BioSample sequenced in three runs**, not three records: all three agree, which is a
consistency check rather than three independent findings. Both carry strain designations (63/183, 63/227),
so a wrong label here travels further than an ordinary isolate's would.

### Single submissions

| Run | BioSample | Study | Submitted as | Both methods say | Distance | ST | PubMLST species of that ST |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `SRR31560207` | `SAMN15685105` | PRJNA649281 | *B. abortus* | *B. melitensis* | 0.00056 | 7 | *B. melitensis*, 98.6% of 72 isolates |
| `SRR11526759` | `SAMN14309819` | PRJNA610713 | *B. abortus* | *B. melitensis* | 0.00098 | 8 | *B. melitensis*, 98.3% of 862 isolates |

`SRR31357708` (`SAMN44729562`, PRJNA1183633), submitted as *B. cytisi*, belongs with these although no ST
is available: the former *Ochrobactrum* species do not carry the classical MLST alleles, and six of the
nine loci got no read at all. The distance evidence is decisive on its own — 0.0076 from the nearest
*B. anthropi* against 0.0381–0.0394 from the three deposits of the *cytisi* type strain, which rank 89th,
90th and 91st of the 431 references. Every one of its top eight hits is *anthropi*. Status: to report.

## Read sets that contain no *Brucella* sequence

Twenty-two read sets carry a *Brucella* species on the label and share almost nothing with any genome of
the clade. These are not identification failures — BrucellaID refused to name them, which is what the
threshold of 150 shared hashes is for — but they are worth reporting, because a reader searching the
archive for *B. abortus* reads will find them.

A real *Brucella* read set at usable depth shares 900 or more of the 1000 hashes with its own species.

| Study | Runs | Submitted as | Shared hashes |
| --- | --- | --- | --- |
| PRJEB4782 (Sanger) | `ERR471337`, `ERR473743`, `ERR554820`, `ERR775236`, `ERR775237`, `ERR775244`, `ERR775245`, `ERR775246`, `ERR775247`, `ERR775248` | *B. abortus* | 1–2 of 1000 |
| PRJNA52097 (Broad) | `SRR825216`, `SRR825217`, `SRR825218` | *B. pinnipedialis* | **0 of 1000** |
| PRJNA197004 (Broad) | `SRR886241`, `SRR886242`, `SRR886243` | *B. ceti* | 1 of 1000 |
| PRJNA77813 (Broad) | `SRR848042` | *B. abortus* | **0 of 1000** |
| PRJNA213741 (Broad) | `SRR952243` | *B. abortus* | 146 of 1000 — just under the threshold, the one borderline case |
| PRJEB22377 (Sanger) | `ERR3218177` | *B. abortus* | 6 of 1000 |
| PRJNA1290416 | `SRR35797377` | *B. melitensis* | 26 of 1000 |
| PRJNA966239 | `SRR24420454` | *B. melitensis* | 105 of 1000 |
| PRJNA945449 | `SRR25828922` | *B. pseudogrignonensis* | 5 of 1000 |

The four runs sharing **zero** hashes are the clearest: whatever those files hold, it is not *Brucella*.
`SRR952243` at 146 of 1000 is the only one near the line and deserves a second look before it is reported.

## Read sets: candidates

Called a species other than the submitter's, but with only one line of evidence. All are former
*Ochrobactrum* species, where the reference side is thin and the within-species spread is wide, so a
distance of 0.005–0.02 is suggestive rather than conclusive.

| Run | Submitted as | Closest | Distance | Shared hashes |
| --- | --- | --- | --- | --- |
| `ERR11807265` | *B. anthropi* | *B. intermedia* | 0.0056 | 801/1000 |
| `SRR10948139` | *B. anthropi* | *B. intermedia* | 0.0087 | 713/1000 |
| `SRR1798201` | *B. lupini* | *B. anthropi* | 0.0089 | 708/1000 |
| `SRR6003475` | *B. thiophenivorans* | *B. pituitosa* | 0.0029 | 890/1000 |

## Read sets where the label looks right and our reference set is what is missing

Kept here so they are not mistaken for mislabels. In each case the submitter is probably correct and
BrucellaID could not place the sample because the extended set has no suitable reference.

* `SRR4038994` (`SAMN05417907`, PRJNA251693), submitted as *B. papionis* — **resolved.** The submitter was
  right, and chasing this run is what turned up `GCA_000163135.1` above: the species' genome was in NCBI
  under *Brucella* sp. all along. With that genome relabelled and added, the run now reports *B. papionis*
  at 0.00083. See [`papionis/README.md`](papionis/README.md).
* `SRR12605656` (*B. haematophila*, 0.0201, 488/1000) and `SRR25828932`, `SRR25828924` (*B. anthropi*,
  0.014 and 0.008) — reported as *Brucella* sp. because the nearest named genome is far away. With three
  *haematophila* genomes in existence, one of them itself a candidate mislabel, the reference side is too
  thin to call these either way.
* `SRR31481153`, submitted as *B. inopinata*, reported as *Brucella* sp. at 0.00053 with 978 of 1000 hashes
  and a `new` sequence type. It sits beside an unnamed lineage representative rather than an *inopinata*
  genome, which probably says more about how few *inopinata* genomes exist than about the label.

## PCR evidence

Every record above was found with a MASH distance and a nine-locus MLST. A third and fourth line of
evidence now sits in [`pcr_evidence/`](pcr_evidence/README.md): *in silico* Bruce-ladder, Suis-ladder, AMOS
and RB51 junction PCR with [brucetyper](https://github.com/duceppemo/brucetyper), and the two real-time
assays of Hoover et al. 2026 ([doi:10.1128/spectrum.02974-25](https://doi.org/10.1128/spectrum.02974-25))
with [insilicoPCR](https://github.com/duceppemo/insilicoPCR), run on every assembly **and** on every raw
read set. 46 assemblies were re-downloaded from NCBI and 48 read sets downloaded, md5-verified and
assembled with SKESA for the purpose, against a control panel of one genome per species of the clade.

In short: the mislabelled read sets are **good *Brucella* genomes with the wrong name on them**, not junk
that confused a distance — 19 of 21 assemble to 3.24–3.35 Mb and 20 of 21 are genus-positive, and of the
17 whose verdict species the multiplex covers, 14 fired the matching species assay. The swapped pair is
confirmed from both directions. 20 of the 22 read sets that appear to hold no *Brucella* are PCR-negative
on both substrates, while the other two turn out to be **mixed libraries** with a minority *Brucella*
population rather than errors. For the former *Ochrobactrum* records the PCR evidence is silent, because
neither assay was designed for that half of the clade; those records rest on the distances alone.

Two findings from the read side are worth carrying into any report. `ERR2938687`, the weakest of the
eleven PRJEB30030 records on distance and the one with the worst assembly, reports *B. suis* in 54% of
informative reads and *B. melitensis* in 33%: **a mixture fits it better than a simple mislabel**, and it
should be described that way rather than as a clean finding. And the two read sets that looked like
exceptions in the "no *Brucella*" group are mixed libraries holding a minority classical *Brucella*
population, not errors in that verdict.

`pcr_evidence/README.md` records the remaining caveat — the *B. suis* assay's specificity depends on exact
primer matching — and the history of a read-mode false negative that this analysis found and that
[insilicoPCR#7](https://github.com/duceppemo/insilicoPCR/issues/7) fixed in 0.10.0.

## Reproducing the evidence

```bash
# reference-to-reference distances (sketched at 10,000 hashes)
mash dist curation.msh curation.msh > allvsall.tsv

# a read set against the extended reference set, keeping every per-reference distance
brucellaid -r1 RUN_1.fastq.gz -r2 RUN_2.fastq.gz -o out --reference-set extended
#   out/brucellaid_mash.tsv holds the distance to all 431 references, ranked
#   out/brucellaid_mlst.tsv holds every allele call with the margin it had
```

The curation decisions that excluded these genomes from the shipped set are in
[`scripts/extended_references_curation.tsv`](../scripts/extended_references_curation.tsv).

## Reporting log

| Record | Reported to | Date | Outcome |
| --- | --- | --- | --- |
| — | — | — | nothing sent yet |

## A note on our own cutoffs

Two of the findings above were missed by the first pass of the curation, both for the same reason, and it
is worth writing down because it will bite again. The rules used **absolute** distance thresholds: a
genome was an intruder if another species was nearer *and* within 0.004, and an unnamed genome was a
redundant twin if a named species was within 0.004.

Those numbers were calibrated on the classical species, which are clonal — two *B. abortus* genomes are
0.0001–0.002 apart. They do not transfer:

* In the former *Ochrobactrum* species the within-species spread is enormous — up to 0.21 inside
  *B. anthropi* — so a mislabel sitting 0.0066 from another species slips past an 0.004 intruder rule.
  That is how `GCA_059995225.1` kept a *cytisi* label while being an *anthropi*.
* In the classical core 0.0020 is a *species-level* gap, not a within-species one, so the twin rule threw
  away a genome that was the sole representative of its species. That is how `GCA_000163135.1`, the only
  *B. papionis* genome anywhere, was dropped as a redundant *B. microti*.

A relative criterion catches both: compare a genome's distance to another species against its distance to
its own, rather than against a fixed number.

The other seven genomes the twin rule dropped at 0.0014–0.0028 have now been re-examined on that
criterion — their distance to the nearest member of the species they were called a twin of, against how
far that species spreads internally.

**Five were dropped correctly.** `GCA_040368495.1` (*pseudogrignonensis*), `GCA_045032515.1`
(*intermedia*), `GCA_058796175.1`, `GCA_047708795.1` and `GCA_040366835.1` (*anthropi*) sit 0.0014–0.0028
from a member of those species, far inside their internal spread, whose median runs 0.021–0.024. For these
species 0.0028 really is within-species variation, so they are redundant twins.

**Two were not twins at all**, and both are in the classical core, where the spread is a hundred times
tighter:

| Genome | NCBI label | Distance | Why it is not a twin | What it is |
| --- | --- | --- | --- | --- |
| `GCA_030157885.1` | *Brucella* sp. bbatCR03, complete genome, 3.28 Mb | 0.00211 to the one *B. microti* genome, 0.00221 to *B. pinnipedialis* | near-equidistant from two species, which is what a species of its own looks like — the *B. papionis* signature | a novel lineage: its MLST profile `26-2-41-2-22-3-23-?-33` matches **no known ST**, and a second assembly, `GCA_002994775.1`, sits 0.00019 from it |
| `GCA_000370945.1` | *Brucella* sp. 63/311, scaffold, 3.35 Mb | 0.00167 to *B. pinnipedialis* | **further than any two *pinnipedialis* genomes are from each other** (maximum 0.00121) | an unnamed lineage both databases agree on: ST 47, which PubMLST also records as *Brucella* sp. |

Unlike *B. papionis*, neither is a named species hiding under an unnamed label — they are genuinely
nameless lineages. Which is precisely what the curation policy says to keep, so that an atypical sample
lands on *Brucella* sp. instead of a confident wrong name. Dropping them worked against that intent: in
the simulated-read validation both were reported as a named species, *B. microti* and *B. pinnipedialis*,
each with a `MASH (unresolved)` flag, rather than as *Brucella* sp.

Both are now back in the set, which stands at 431 genomes. Measured against all 2,022 assemblies, exactly
three samples change, every one of them labelled *Brucella* sp. in NCBI: the two genomes themselves, and
`GCA_002994775.1`, which stops being called *B. microti* at 0.00224 and matches the bbatCR03 lineage at
0.00019 — reported as *Brucella* sp., which is the honest answer. **No sample of a named species moves at
all**, and a *B. microti* and a *B. pinnipedialis* read set return the same species, closest reference and
distance as before.

## Scope

The read-set tables cover the completed Phase 2 sweep: every usable paired Illumina WGS run under NCBI
taxon 234 in ENA as of 2026-11-03, **4,042 of 4,044 analysed**. The two that are missing,
`SRR9831875` and `SRR3096373`, could not be downloaded from ENA after two passes; neither failed to
analyse. Fifty-two runs failed on the first pass and a retry recovered fifty of them, which added one
finding: `SRR886241`, included above.

The assembly tables come from the curation of the extended reference set and cover the 454 dereplicated
representatives of the mashID pack source collection, not all 2,022 assemblies, so the assembly side is
less complete than the read side.
