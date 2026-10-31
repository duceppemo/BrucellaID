# Reference data

Five files, about 300 kB. `brucellaid --sketch/--references/--alleles/--profiles` take files of your own
instead of these.

| File | What it is | Where it comes from |
| --- | --- | --- |
| `references.msh` | MASH sketch of 16 *Brucella* reference genomes (k = 21, 1000 hashes) | The *Brucella* part of [vSNP](https://github.com/USDA-VS/vsnp)'s reference set, as `vsnp_reference.msh` |
| `references.tsv` | What each genome of the sketch is: vSNP species code, species, biovar, strain, NCBI accessions | Written for BrucellaID, from the sketch and vSNP's `species_accessions.csv` |
| `mlst_alleles.fasta` | 309 alleles of the nine MLST loci | [PubMLST](https://pubmlst.org/organisms/brucella-spp), "MLST (9 locus)", downloaded 2026-10-31 |
| `mlst_profiles.tsv` | 122 sequence types | The same scheme, the same day |
| `mlst_st_species.tsv` | The species of each of 118 of those STs, with how consistently its isolates carry it | The [PubMLST *Brucella* isolate database](https://pubmlst.org/organisms/brucella-spp), 3426 isolates, read 2026-10-31 |

The sketch ID of a genome is the name of the file it was built from, which is what `mash dist` reports and
what `references.tsv` keys on.

Two genomes of the sketch, `Bceti1Cudo.fasta` (*B. ceti* strain Cudo) and `B-REF-BS4-40.fasta` (*B. suis*
bv. 4 strain 40), are distributed by vSNP under its own names rather than as accessions, so their
`accessions` column is empty and `scripts/make_reference_data.py` cannot download them.

## Rebuilding
```bash
python3 scripts/make_mlst_data.py                                   # the two MLST files, from PubMLST
python3 scripts/make_reference_data.py --genomes /tmp/refs --out /tmp/references.msh   # the sketch
python3 scripts/make_st_species.py                                  # the ST to species table
```
PubMLST serves the data submitted up to a cut-off date unless you authenticate; the script prints the date
it was given. The MLST data here was downloaded on 2026-10-31, with data submitted up to 2024-12-31.

The isolate database has no bulk export, so `make_st_species.py` reads its 3426 records one at a time
with a pause between them -- about forty minutes -- and caches them, so a second run is free. The
species and biovar labels in it are the submitters', not the curators', which is why the table
reports a purity column instead of a single name: 105 of the 118 STs are seen in one species only, and
the thirteen that are not disagree by one to three isolates out of dozens.

The scheme is Whatmore, Perrett and MacMillan (2007),
[doi:10.1186/1471-2180-7-34](https://doi.org/10.1186/1471-2180-7-34). PubMLST is maintained at the
University of Oxford; publications using its data should cite Jolley, Bray and Maiden (2018),
[doi:10.12688/wellcomeopenres.14826.1](https://doi.org/10.12688/wellcomeopenres.14826.1).

A sketch must keep MASH's defaults (k = 21, 1000 hashes): the threshold of 150 shared hashes
(`--min-hashes`) is read against that sketch size.
