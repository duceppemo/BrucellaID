# Reference data

Nine files, about 4 MB — 3.5 MB of it the extended sketch. `brucellaid --sketch/--references/--alleles/
--profiles/--st-species` take files of your own instead of the first five, and `--reference-set extended`
switches to the extended pair; the logo is only used by the reports.

| File | What it is | Where it comes from |
| --- | --- | --- |
| `references.msh` | MASH sketch of 16 *Brucella* reference genomes (k = 21, 1000 hashes) | The *Brucella* part of [vSNP](https://github.com/USDA-VS/vsnp)'s reference set, as `vsnp_reference.msh` |
| `references.tsv` | What each genome of the sketch is: vSNP species code, species, biovar, strain, NCBI accessions | Written for BrucellaID, from the sketch and vSNP's `species_accessions.csv` |
| `references_extended.msh` | MASH sketch of 429 curated genomes covering the whole *Brucella* clade (k = 21, 1000 hashes), used by `--reference-set extended` | The assemblies behind the [mashID *Brucella* pack](https://doi.org/10.6084/m9.figshare.33970168), dereplicated at 0.0003 and curated: `scripts/extended_references_curation.tsv` logs every decision |
| `references_extended.tsv` | What each genome of the extended sketch is. The code is the species epithet (plus the biovar when NCBI records one); only 7 genomes have a biovar and 349 a strain, because the NCBI records of the rest carry none | Generated from the NCBI assembly data report of the 429 accessions |
| `mlst_alleles.fasta` | 309 alleles of the nine MLST loci | [PubMLST](https://pubmlst.org/organisms/brucella-spp), "MLST (9 locus)", downloaded 2026-10-31 |
| `mlst_profiles.tsv` | 122 sequence types | The same scheme, the same day |
| `mlst_st_species.tsv` | The species of each of 118 of those STs, with how consistently its isolates carry it | The [PubMLST *Brucella* isolate database](https://pubmlst.org/organisms/brucella-spp), 3426 isolates, read 2026-10-31 |
| `logo.svg` | The logo the HTML report inlines, with its four colours as CSS variables so that one drawing serves the light and the dark page | `assets/make_logo.py` |
| `logo.png` | The same logo for the PDF report, which draws a png without needing anything that reads svg | `assets/make_logo.py` |

The sketch ID of a genome is the name of the file it was built from, which is what `mash dist` reports and
what `references.tsv` keys on.

`mlst_st_species.tsv` is what decides the one species pair the sketch cannot: *B. canis* ATCC 23365 and
*B. suis* bv. 4 strain 40 are 0.00024 apart, and the two species share no ST.

Two genomes of the sketch, `Bceti1Cudo.fasta` (*B. ceti* strain Cudo) and `B-REF-BS4-40.fasta` (*B. suis*
bv. 4 strain 40), are distributed by vSNP under its own names rather than as accessions, so their
`accessions` column is empty and `scripts/make_reference_data.py` cannot download them.

## Rebuilding
```bash
python3 scripts/make_mlst_data.py                                   # the two MLST files, from PubMLST
python3 scripts/make_reference_data.py --genomes /tmp/refs --out /tmp/references.msh   # the sketch
python3 scripts/make_extended_references.py --genomes /tmp/extended_refs   # the extended sketch
python3 scripts/make_st_species.py                                  # the ST to species table
```
PubMLST serves the data submitted up to a cut-off date unless you authenticate; the script prints the date
it was given. The MLST data here was downloaded on 2026-10-31, with data submitted up to 2024-12-31.

The isolate database has no bulk export, so `make_st_species.py` reads its 3426 records one at a time
with a pause between them — about forty minutes — and caches them, so a second run is free. The
species and biovar labels in it are the submitters', not the curators', which is why the table
reports a purity column instead of a single name: 105 of the 118 STs are seen in one species only, and
the thirteen that are not disagree by one to three isolates out of dozens.

The scheme is Whatmore, Perrett and MacMillan (2007),
[doi:10.1186/1471-2180-7-34](https://doi.org/10.1186/1471-2180-7-34). PubMLST is maintained at the
University of Oxford; publications using its data should cite Jolley, Bray and Maiden (2018),
[doi:10.12688/wellcomeopenres.14826.1](https://doi.org/10.12688/wellcomeopenres.14826.1).

A sketch must keep MASH's defaults (k = 21, 1000 hashes): the threshold of 150 shared hashes
(`--min-hashes`) is read against that sketch size.
