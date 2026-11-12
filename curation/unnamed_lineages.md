# The read sets reported as *Brucella* sp.: which lineages they are, and whether any has a name

The 2.1.1 sweep of the extended set calls **43 of 4,024** public read sets *Brucella* sp. — inside the
clade, matching no named species. *B. papionis* was found by chasing exactly that kind of call
([`papionis/README.md`](papionis/README.md)), so all 43 were chased the same way: assembled, grouped,
and compared by ANI with the type strains, with the nearest named genomes and with each other.

**Nothing in the reference set was changed by this work.** Two findings would justify a change, and
both are relabels rather than additions; they are listed under [What was not changed](#what-was-not-changed)
for a decision, because a relabel changes what a sample is called.

| Verdict | Lineages | Runs |
| --- | --- | --- |
| A named species the set carries only under *Brucella* sp. | *Ochrobactrum teleogrylli* | 3 |
| A lineage with a published name that is not validly published | "*B. nosferati*", the vampire-bat lineage | 22 |
| Unnamed lineages, each already represented in the set | eleven | 18 |
| Within-species variation of a named *Brucella* species | none | 0 |

## The candidates

[`unnamed_lineages/candidates.tsv`](unnamed_lineages/candidates.tsv) has one row per run;
[`unnamed_lineages/make_candidates.py`](unnamed_lineages/make_candidates.py) rebuilds the list from the
sweep output, which is not shipped.

All 43 have a SKESA assembly. Their median is 3.26 Mb in 37 contigs, and 39 are "clean and right-sized"
(at most 150 contigs, within 0.5 Mb of 3.3 Mb). The four that fail are not four bad assemblies: three are
former *Ochrobactrum* genomes of 4.9–5.4 Mb, which is the right size for what they are, and only
`SRR33001687` (230 contigs, from 462,642 reads) is a shallow library.

**The premise that these lineages have no good reference does not hold for most of them.** 37 of the 43
sit within 0.0006 of a genome already in the set, and 40 within 0.0016. They are called *Brucella* sp.
because that reference is itself unnamed, which is what the set is built to do. Only the three former
*Ochrobactrum* runs are far from everything (0.0077, 0.0137, 0.0200).

By closest reference the 43 are: 22 on `GCA_030157885.1`; 3 each on `GCA_017797885.1` and
`GCA_024970695.1`; 2 each on `GCA_000177135.1`, `GCA_000370945.1`, `GCA_014084005.1` and
`GCA_040372265.1`; and one each on seven others. (The tally this work started from had 2 on
`GCA_030157885.1`. It is 22.)

**Independent libraries are fewer than runs.** Several clusters are one strain sequenced more than once,
or reads of the very strain that is the reference, and those corroborate nothing:

* `SRR642174` and `SRR647550` are one BioSample, strain 63/311, which is also the reference.
* `ERR775238` carries the sample alias `BCCN_84_3`, so it is the same dog isolate as `SRR5927194`.
* `SRR4038983` and `SRR13170268` are two laboratories' reads of the one strain BO2.
* `SRR9620238`, `SRR12630907`, `SRR28969979`, `SRR14074356` and `SRR19147302` are reads of strains 2280,
  6810, CVUAS_1139.3, 458 and bbatCR03, each of which is the reference it matches.

## What ANI can and cannot decide here

Every ANI below is fastANI 1.33 on assemblies; where two genomes were compared in both directions the
mean is given. Genome-to-genome Mash distances are at 10,000 hashes, read-to-reference distances are the
sweep's own at 1,000.

In the former *Ochrobactrum* half of the clade, ANI to the type strain decides species, as it did in
[`mislabelled_records.md`](mislabelled_records.md). **In the classical half it decides nothing**, and the
measurements say so directly ([`unnamed_lineages/lineage_matrix.tsv`](unnamed_lineages/lineage_matrix.tsv),
27 genomes all against all):

| Comparison | ANI |
| --- | --- |
| Named species of the classical core, type strain to type strain | 99.53–99.89% |
| *B. inopinata* BO1 to each of those | 97.57–97.70% |
| *B. vulpis* F60 to every other genome, BO1 included | 97.17–97.63% |
| The ten unnamed early-diverging lineages, lineage to lineage | 97.42–98.32% |

Every *Brucella* of the classical clade, named or not, clears 95% against every type strain in it. A
candidate there cannot be assigned to a species by ANI, and it cannot be excluded from one. What ANI and
distance can still show is **relative**: whether a lineage is as far from its nearest named species as
the named species are from each other. That is the criterion that separated *B. papionis*, and it is the
one used for the classical clade below.

## A named species: *Ochrobactrum teleogrylli*

`SRR12605656`, `SRR25828924` and `SRR25828932` are one species, and it has a name — but the name is in
genus *Ochrobactrum*, which NCBI files outside taxon 234, so the collection the reference set was curated
from never contained its type strain. All 126 GenBank assemblies under NCBI taxon 528 were downloaded on
2026-11-12 and compared ([`unnamed_lineages/teleogrylli_ani.tsv`](unnamed_lineages/teleogrylli_ani.tsv)):

| Run | BioSample | Study | Deposited | ANI to its own type | ANI to *O. teleogrylli* (CGMCC 1.13984) | ANI to *O. soli* (BO-7) | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `SRR12605656` | `SAMN16077189` | PRJNA662158 | *B. haematophila* | 89.86% (CIP 109452) | **97.99%** | 97.87% | *O. teleogrylli* |
| `SRR25828924` | `SAMN33797214` | PRJNA945449 | *B. anthropi* | 83.00% (ATCC 49188) | 98.14% | **99.05%** | *O. teleogrylli* |
| `SRR25828932` | `SAMN33797175` | PRJNA945449 | *B. anthropi* | 83.03% (ATCC 49188) | 98.14% | **98.95%** | *O. teleogrylli* |

Neither deposited label survives: all three are 8 to 16 ANI points nearer these type strains than their
own. Status: **to report**; nothing has been sent.

*O. teleogrylli* and *O. soli* are one species by the same criterion — their type strains are **98.16%**
apart, and the two deposits of the *teleogrylli* type (LCB8, CGMCC 1.13984) agree at 99.99%. LPSN, read on
2026-11-12, gives *O. teleogrylli* Hu *et al.* 2020 as validly published in IJSEM in 2020 and *O. soli*
Choi *et al.* as effectively published in 2020 but validated only in 2024 (Validation List no. 216), so
*teleogrylli* has priority. LPSN lists both as orphaned species and has no *Brucella* combination for
either.

The set already holds three genomes of this species, all labelled *Brucella* sp.: `GCA_030063905.1` (NM4,
99.17% to the *soli* type), `GCA_040372265.1` (NBRC 113783, 98.99%) and `GCA_040391195.1` (98.90%, though
it is a 2.1 Mb partial genome with 636 of 686 fragments aligned). That is why these runs are reported as
*Brucella* sp. and not as a wrong named species.

Their read distances are large because the species is diverse, not because a reference is missing: its
members here are 97.6–99.5% ANI apart. Adding the type strain would give the species a name in the set
but would bring `SRR12605656` no closer — it is 97.99% from the type against 97.85% from NM4.

This replaces the earlier reading in `mislabelled_records.md`, which held that the submitters of these
three were probably right and the reference side was too thin to say.

## A lineage with a published but not validly published name: "*B. nosferati*"

Twenty-two runs land on `GCA_030157885.1` (*Brucella* sp. bbatCR03): twenty isolates from one colony of
common vampire bats in Costa Rica (PRJNA836284, 2020–2021) and two libraries of BCCN 84-3, from a dog in
Costa Rica in 1984. So there are **two independent sources 36 years apart**, not 22.

| Evidence | Result |
| --- | --- |
| Bat isolates to each other, assemblies | 0.00000–0.00003 |
| Bat isolates to the 1984 dog isolate | 0.00018–0.00030 |
| All 22 runs, nine-locus profile | `26-2-41-2-22-3-23-21-33`, no known ST |
| ANI to the nearest type strain | 99.73% *B. microti* CCM 4915; 99.51–99.71% to the whole classical core |
| bbatCR03 to its nearest named species, genome to genome | 0.00211 *B. microti*, 0.00235 *B. ceti*, 0.00238 *B. amazoniensis* |
| The same measure for named species | 0.00131 *ceti*–*microti*, 0.00147 *pinnipedialis*–*ceti*, 0.00153 *amazoniensis*–*microti*, 0.00200 *papionis*–*microti* |

It is further from every named species than *B. papionis*, *B. ceti*, *B. pinnipedialis* and
*B. amazoniensis* each are from their nearest neighbour, so on the relative criterion it is a lineage of
species rank.

It has been given a name. Hernández-Mora *et al.*, mSphere 2023
([doi:10.1128/msphere.00061-23](https://doi.org/10.1128/msphere.00061-23)) describe *Brucella nosferati*
sp. nov. from these bats, and their abstract places it as the agent of "a reported canine brucellosis
case", which fits BCCN 84-3. **The name is not validly published**: LPSN lists it in quotation marks with
the note that the taxon "lacks a discernible description (protologue)", NCBI Taxonomy has no entry for
it, and `GCA_030157885.1` is still *Brucella* sp. bbatCR03 there.

So this is not the *B. papionis* situation. That was a validly published species whose only genome sat
under an unnamed label. This is an unnamed-in-the-archives lineage whose reference is already in the set
and already catches all 22 runs, 18 of them at a distance of zero.

## Unnamed lineages, each already represented

None of these matches a named species on the relative criterion, each is as distinct from the others as
*B. inopinata* and *B. vulpis* are, and each already has a representative in the set. Nothing needs
adding.

| Lineage, and its reference in the set | Runs | Independent isolates among them | ST | Read distance | Nearest named species, genome to genome |
| --- | --- | --- | --- | --- | --- |
| 63/311, `GCA_000370945.1` | 2 | none: one BioSample, and it is the reference strain | 47 | 0.00017 | 0.00172 *B. microti*, 0.00175 *B. pinnipedialis* |
| 1315 / 2594, `GCA_024970695.1` | 3 | 3, amphibians from California, Massachusetts and Florida, 2018–2022; the reference is from a frog in Vienna | new | 0.00005–0.00012 | 0.0178 *B. pinnipedialis* |
| 458, `GCA_017797885.1` | 3 | 1 beyond the reference strain: two Wisconsin amphibian isolates 0.00001 apart, counted as one | 76 | 0.00010; 0.00144–0.00150 | 0.0169 a genome labelled *inopinata* |
| BO3 / B13-0095, `GCA_014084005.1` | 2 | 2: a patient in France, 2019, and a frog in Texas, 2010 | 80 | 0.00005–0.00029 | 0.0189 a genome labelled *inopinata* |
| BO2, `GCA_000177135.1` | 2 | none: one strain, two laboratories | new | 0.00007–0.00019 | 0.0161 *B. inopinata* |
| 09RB8918, `GCA_009664935.1` | 1 | 1, an amphibian in California, 2017 | 65 | 0.00106 | 0.0176 a genome labelled *inopinata* |
| 191011898, `GCA_903131625.1` | 1 | 1, a marine toad in Michigan, 2020, submitted as *B. inopinata* | new | 0.00043 | 0.0176 a genome labelled *inopinata* |
| NF 2653 / 83/13, `GCA_050410765.1` | 1 | none: strain NF 2653, of which an assembly exists | 50 | 0.00007 | 0.0207 *B. pinnipedialis* |
| 2280, `GCA_009601725.1` | 1 | none: the reference strain | 85 | 0.00014 | 0.0182 *B. inopinata* |
| 6810, `GCA_014495905.1` | 1 | none: the reference strain | 100 | 0.00007 | 0.0201 *B. inopinata* |
| CVUAS_1139.3, `GCA_050600715.1` | 1 | none: the reference strain | new | 0.00005 | 0.0190 *B. inopinata* |

PubMLST records STs 47, 50, 65, 76, 80, 85 and 100 as *Brucella* sp. in every isolate that carries them,
so the second method agrees these are unnamed. All 43 calls are reported as plain `MASH`.

Three of these are worth a sentence each.

**1315 / 2594 is the best-corroborated unnamed lineage in the sweep.** Five isolates from two continents
over seven years — two Austrian assemblies and three American read sets — sit within 0.0002 of each other
with one new nine-locus profile, `28-29-42-22-24-25-25-24-36`. It is 97.53–97.99% from every other genome
in the matrix and nearer the classical core than the other early-diverging lineages are. It is a
candidate for a species description, by someone with the strains.

**The 458 lineage has a second branch.** The two Wisconsin isolates share ST 76 with strain 458 and are
99.92–99.94% ANI from it, but their assemblies are 0.0016 from the reference and their reads match at
0.0014–0.0015, past the 0.001 where the sweep attaches a loose-match note. See below for why that is left
alone.

**09RB8918 and 191011898 are two branches of one group**, 99.43% ANI and 0.0059 apart. Both are in the
set, and each candidate lands on its own branch.

## Two things found on the way

Neither concerns the 43 candidates directly; both came out of the same measurements.

**A genome the set labels *B. inopinata* is not nearer the type strain than the unnamed lineages are.**
`GCA_900095155.1` (strain 141012304) is 97.87% ANI and 0.0212 from BO1. The unnamed BO2, 2280 and BO3
references are 98.24%, 98.08% and 98.01% from BO1, so three genomes the set calls *Brucella* sp. are
closer to the *inopinata* type than a genome it calls *inopinata*. The third *inopinata* genome,
`GCA_030270765.1`, is a real one: 99.88% and 0.0017 from BO1. The label has a consequence in the sweep:
`SRR12224742`, submitted as *Brucella* sp. DSM 103976, is reported as *B. inopinata* at 0.00012 because
`GCA_900095155.1` is what it matches. `SRR31481153`, submitted as *B. inopinata*, is the mirror image —
97.84% from BO1 and 99.98% from the unnamed 191011898. "*inopinata*" is evidently being used in the
archive for early-diverging *Brucella* generally, and the set has inherited one instance of that.

**19 of the set's 25 unnamed former-*Ochrobactrum* references are at or above 95% ANI to a type strain**
([`unnamed_lineages/unnamed_references_vs_types.tsv`](unnamed_lineages/unnamed_references_vs_types.tsv)):
14 to a named *Brucella* species, the three *O. teleogrylli* genomes above, and `GCA_039674265.1` and
`GCA_039674445.1` at 98.5% to the type of *Ochrobactrum chromiisoli*. A read set landing on one of them
is reported as *Brucella* sp. although its species has a name. None of the 4,024 public runs does, so
this cost nothing in the sweep; it is recorded because it is the same cause as the two cutoffs described
under "A note on our own cutoffs" in `mislabelled_records.md` — an absolute distance of 0.004 does not
describe species whose members are 0.02 apart.

## What the evidence does not settle

* **Whether the bat lineage should be called *B. nosferati*.** The genomes say it is a distinct lineage
  of species rank by the clade's own spacing. Whether the name stands is nomenclature, and today it does
  not. The paper itself could not be read from here (the publisher refused the request), so which strain
  it designates as type, and whether bbatCR03 is that strain, is not verified; only the abstract and the
  LPSN record were.
* **Whether any early-diverging lineage is a species.** They are as far apart as *B. inopinata* and
  *B. vulpis*, both of which are validly published, but LPSN gives those two the taxonomic status
  "synonym", and the whole classical clade is one species by the 95% criterion. "As distinct as a named
  species" is the most the data supports.
* **Whether the 458 lineage wants a second representative.** 0.0016 is a tenth of the smallest gap
  between early-diverging lineages (0.0154), and it is also wider than the gap between *B. ceti* and
  *B. microti* (0.0013). Which yardstick applies inside a lineage of three isolates cannot be said from
  three isolates. There is also nothing to add: BioProject PRJNA903886 has no assembly in GenBank, and
  the set takes accessioned genomes only.
* **Whether *O. soli* is a synonym of *O. teleogrylli*.** 98.16% between the type strains says one
  species; nobody has published that, and LPSN lists no synonymy. The three runs are assigned to the
  senior name on that ANI, which is this note's inference.
* **Independence inside PRJNA903886.** Six candidates come from that one submission, and its two
  Wisconsin isolates are 0.00001 apart. They are counted as one source here; the metadata does not say
  whether they are one animal.
* **The assemblies were built from reads subsampled to 100x and are unpolished.** That is ample for ANI
  and for the distances quoted, and none of the verdicts rests on a margin an assembly error could close.

## What was not changed

A reference's label decides what a sample is called, so each of these is a release-level decision and
none was made here.

| Change | Evidence | What it would do |
| --- | --- | --- |
| Label `GCA_030063905.1`, `GCA_040372265.1` and `GCA_040391195.1` as *O. teleogrylli*, with or without adding the type strain `GCA_042681285.1` | 98.9–99.2% ANI to the *soli* type, 98.1–98.5% to the *teleogrylli* type | Three runs of the sweep would get a species. It would be the set's first name outside genus *Brucella*, since no *Brucella* combination exists. |
| Label `GCA_030157885.1` as "*B. nosferati*" | A lineage of species rank with a published name | 22 runs would get that name. The set would then carry a name NCBI does not have. *B. amazoniensis* is a partial precedent: LPSN has no entry for it either, but NCBI has a taxon. |
| Relabel `GCA_900095155.1` from *B. inopinata* to *Brucella* sp. | 97.87% from BO1, further than three unnamed references | `SRR12224742` would be reported as *Brucella* sp., as submitted. |
| Name the 16 other unnamed references that clear 95% to a type strain | `unnamed_lineages/unnamed_references_vs_types.tsv` | No run of the sweep changes. |

## Reproducing the evidence

```bash
export PATH=/home/bioinfo/miniconda3/envs/brucellaid/bin:$PATH

# the candidate list, from the sweep output and the assemblies
python3 -I curation/unnamed_lineages/make_candidates.py > candidates_raw.tsv

# distances between the assemblies, and to all 2,022 assemblies of the source collection
mash sketch -k 21 -s 10000 -l assemblies.txt -o cand
mash dist cand.msh cand.msh > cand_vs_cand.tsv
mash dist all2022.msh cand.msh > cand_vs_2022.tsv

# ANI, in the fastani conda environment (fastANI 1.33)
fastANI --ql assemblies.txt --rl panel.txt -t 48 -o cand_vs_panel.tsv
```

The raw fastANI and Mash output, the sketches, the 126 *Ochrobactrum* assemblies and the scripts that
turned them into the tables below were kept at `/media/8tb/analyses/brucellaid_unnamed_lineages/`, which
is not shipped.

The panel was the type strain of every species of the clade, every genome of *B. pinnipedialis*,
*B. microti*, *B. inopinata*, *B. vulpis*, *B. amazoniensis* and *B. papionis*, all 94 genomes NCBI labels
*Brucella* sp., and the candidates themselves. `GCA_000662015.2` was left out: it is the swapped record
that carries the BO1 label on a *suis*/*canis* genome.

| File in [`unnamed_lineages/`](unnamed_lineages/) | Holds |
| --- | --- |
| `candidates.tsv` | the 43 runs: lineage, BioSample, study, host, the sweep's call, assembly size, ANI to the reference and to the best type strain |
| `candidates_ani.tsv` | each candidate against every type strain and named genome of the panel, and every unnamed genome at 99% or above |
| `lineage_matrix.tsv` | ANI and Mash distance, all against all, for 27 type strains and lineage representatives |
| `teleogrylli_ani.tsv` | the three former *Ochrobactrum* runs and the set's three genomes of that species against the *Ochrobactrum* type strains |
| `unnamed_references_vs_types.tsv` | the set's 25 unnamed former-*Ochrobactrum* references against the type strains |
| `make_candidates.py` | rebuilds the candidate list |
