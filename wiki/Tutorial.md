# Tutorial

Four public *Brucella* read sets, from downloading them to reading the reports: about ten minutes, most of
it the download. Every command and every output below is a real run.

## 1. Install
```bash
conda create -n brucellaid -c conda-forge -c bioconda mash "gsl=2.7.0" bbmap python">=3.10"
conda activate brucellaid
pip install brucellaid
brucellaid --version
```
```
brucellaid 2.0.0
```
One `conda install brucellaid` will do all of it once the
[bioconda recipe](https://github.com/bioconda/bioconda-recipes/pull/71040) is merged. See
[Installation](Installation).

## 2. Get four samples
Four runs from the ENA, one per species, about 540 MB in all:

| Run | Species, as the submitter labelled it | Size |
| --- | --- | --- |
| [SRR19520329](https://www.ebi.ac.uk/ena/browser/view/SRR19520329) | *Brucella melitensis* | 80 MB |
| [ERR3454415](https://www.ebi.ac.uk/ena/browser/view/ERR3454415) | *Brucella abortus* | 179 MB |
| [SRR8550484](https://www.ebi.ac.uk/ena/browser/view/SRR8550484) | *Brucella suis* | 151 MB |
| [ERR2136546](https://www.ebi.ac.uk/ena/browser/view/ERR2136546) | *Brucella canis* | 127 MB |

```bash
mkdir -p brucella_tutorial/reads && cd brucella_tutorial
for run in SRR19520329 ERR3454415 SRR8550484 ERR2136546; do
  urls=$(curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" \
         | tail -1 | cut -f2 | tr ';' ' ')
  for url in $urls; do curl -fsS "https://$url" -o "reads/$(basename "$url")"; done
done
ls reads/
```
```
ERR2136546_1.fastq.gz  ERR3454415_1.fastq.gz  SRR19520329_1.fastq.gz  SRR8550484_1.fastq.gz
ERR2136546_2.fastq.gz  ERR3454415_2.fastq.gz  SRR19520329_2.fastq.gz  SRR8550484_2.fastq.gz
```
The `_1` and `_2` files are paired automatically, and the sample name is what is left of the file name.

## 3. Run it
```bash
brucellaid -i reads -o results -t 8 -p 2
```
```
432 reference genome(s) in .../brucellaid/data/references_extended.msh
9 MLST locus/loci: gap, aroA, glk, dnaK, gyrB, trpE, cobQ, int_hyp, omp25
118 sequence type(s) with a species, to decide the species MASH cannot
4 sample(s) found in reads
4 sample(s), 2 at a time, 4 thread(s) each
ERR2136546: identifying the species with MASH (2 file(s))
ERR3454415: identifying the species with MASH (2 file(s))
ERR3454415: about 230,479,484 bases, roughly 70x coverage: sketching the k-mers seen at least 14 time(s)
ERR2136546: about 153,926,473 bases, roughly 47x coverage: sketching the k-mers seen at least 9 time(s)
ERR2136546: Brucella canis str. Oliveri (GCA_000530495.1_Oliveri_genomic.fna.gz, MASH distance 0.00012, 995/1000 hashes)
ERR2136546: calling the MLST alleles with Seal
ERR2136546: ST 20 (profile 1-6-4-1-5-3-5-4-6)
ERR2136546: Brucella canis (MASH+ST, arbitrated)
ERR3454415: Brucella abortus str. 2308 (GCA_000742195.1_ASM74219v1_genomic.fna.gz, MASH distance 0.00014, 994/1000 hashes)
ERR3454415: ST 1 (profile 2-1-1-2-1-3-1-1-1)
ERR3454415: Brucella abortus (MASH+ST, agree)
SRR19520329: about 174,716,055 bases, roughly 53x coverage: sketching the k-mers seen at least 11 time(s)
SRR19520329: Brucella melitensis str. 48_S36 (GCA_061325495.1_ASM6132549v1_genomic.fna.gz, MASH distance 0.00012, 995/1000 hashes)
SRR19520329: ST 11 (profile 3-2-3-2-1-5-3-2-10)
SRR8550484: about 229,332,903 bases, roughly 69x coverage: sketching the k-mers seen at least 14 time(s)
SRR8550484: Brucella suis str. Thompsen (GCA_000742055.1_ASM74205v1_genomic.fna.gz, MASH distance 0.00026, 989/1000 hashes)
SRR8550484: ST 15 (profile 1-2-7-1-3-3-5-3-2)
Done in 33s
```
33 seconds on 8 threads, for 4.6 million read pairs. `-p 2` runs two samples at a time and gives each of
them four of the eight threads. The `about ... coverage` line of each sample is the k-mer copy threshold being derived
from each sample’s depth — 14 copies for the ~70x libraries, 11 at 53x and 9 at 47x — which is what keeps
a deep read set’s sequencing errors out of the sketch. Section 6 shows what that is worth, and
[The k-mer copy threshold](K-mer-threshold) is where those numbers come from.

## 4. Read the main table
`results/brucellaid.tsv`, also printed to the terminal:

```
Sample       Species              SpeciesSource         ClosestReference                             SharedHashes  SequenceType  AlleleProfile       STSpecies            Reads
ERR2136546   Brucella canis       MASH+ST (arbitrated)  GCA_000530495.1_Oliveri_genomic.fna.gz       995/1000      20            1-6-4-1-5-3-5-4-6   Brucella canis       625576
ERR3454415   Brucella abortus     MASH+ST (agree)       GCA_000742195.1_ASM74219v1_genomic.fna.gz    994/1000      1             2-1-1-2-1-3-1-1-1   Brucella abortus     1175082
SRR19520329  Brucella melitensis  MASH+ST (agree)       GCA_061325495.1_ASM6132549v1_genomic.fna.gz  995/1000      11            3-2-3-2-1-5-3-2-10  Brucella melitensis  1136924
SRR8550484   Brucella suis        MASH+ST (agree)       GCA_000742055.1_ASM74205v1_genomic.fna.gz    989/1000      15            1-2-7-1-3-3-5-3-2   Brucella suis        1693672
```
Four species, four sequence types, and each species is the one the submitter labelled the run with.
`MashDistance`, left out above to fit, reads 0.000119, 0.000144, 0.000119 and 0.000264. There is no biovar
column: see [How it works](How-it-works#no-biovar).

The reference names are GenBank accessions because a run identifies against the 432-genome clade set by
default from 2.0. `--reference-set vsnp` gives vSNP's 16 genomes and its own file names instead, and on
these four samples the same four species.

**The first sample is the interesting one.** Its `SpeciesSource` says `arbitrated`, and the `Notes` column,
also left out above, says why:

```
nearly as close: GCA_000741995.1_ASM74199v1_genomic.fna.gz (Brucella suis, 986/1000 hashes), 3 other
Brucella canis genome(s); MASH cannot separate Brucella canis and Brucella suis: ST 20 is Brucella canis
in 70 of 70 isolate(s), which agrees with the closest reference
```
Two things are being said there, and both matter.

*B. canis* came first by **nine hashes out of a thousand** — 0.000119 against 0.000337, a margin that is
noise, not a measurement. So the species was not decided by that margin: ST 20 was, and ST 20 is *B. canis*
in all 70 isolates of the PubMLST isolate database that carry it, and never *B. suis*. Had the distances
come out the other way round, as they easily could have, the answer would have been the same. See
[How it works](How-it-works).

And **three other *B. canis* genomes were nearly as close as well**, counted rather than listed. In a set
of 432 a cluster of same-species genomes crowds the top of the table and listing them all would bury the
one line that matters, which is the *B. suis* genome nine hashes behind.

None of these four samples carries a **loose match** note, which is the other note this column can hold:
it appears when the closest reference is further than `--loose-distance`, 0.005 by default, and says that
the reads are not close to anything in the set rather than that the species is wrong. Earlier versions of
this page showed two of these libraries carrying one, at 0.0065. That was the k-mer copy threshold being
too low rather than anything about the libraries, and
[The k-mer copy threshold](K-mer-threshold) is what came of it.

## 5. Read what the calls rest on
`results/brucellaid_mash.tsv` has every distance. For the *B. canis* sample, the first four:

```
Sample      Reference                                    Species          Biovar  MashDistance  SharedHashes
ERR2136546  GCA_000530495.1_Oliveri_genomic.fna.gz       Brucella canis           0.000119      995
ERR2136546  GCA_000741995.1_ASM74199v1_genomic.fna.gz    Brucella suis            0.000337      986
ERR2136546  GCA_056582615.1_ASM5658261v1_genomic.fna.gz  Brucella canis           0.000361      985
ERR2136546  GCA_000367265.1_Bruc_cani_CNGB_1172_V1...    Brucella canis           0.000483      980
```
The *B. suis* genome in second place is the whole reason this sample needs the sequence type. Rows three
and four are the other *B. canis* genomes the note counted.

`results/brucellaid_mlst.tsv` has the nine allele calls and the margin of each. For the *B. abortus*
sample:

```
Sample      Locus    Allele  Reads  RunnerUp  RunnerUpReads
ERR3454415  gap      2       432    6         300
ERR3454415  aroA     1       322    30        210
ERR3454415  glk      1       222    21        158
ERR3454415  dnaK     2       340    26        222
ERR3454415  gyrB     1       334    22        224
ERR3454415  trpE     3       320    4         242
ERR3454415  cobQ     1       308    12        196
ERR3454415  int_hyp  1       368    4         256
ERR3454415  omp25    1       292    12        186
```
Each allele was the best match for more reads than any other allele of its locus, by about 100 reads. That
is what a solid call looks like; a call one or two reads ahead of its runner-up is one to look at twice.

## 6. One sample that shows what the depth handling is worth
`SRR8550484` is the deepest library here: 1.7 million read pairs, about 69x. At that depth sequencing
errors turn up often enough to pass a fixed two-copy filter, and the same reads against the same reference
have been reported three different ways as that filter was corrected:

| what the threshold was | copies | shared hashes | distance |
| --- | --- | --- | --- |
| a fixed 2, up to 1.1.0 | 2 | 422/1000 | 0.0248 |
| one more copy per 20x, 1.2.0 to 2.0.0 | 3 | 775/1000 | 0.006455 |
| one more copy per 5x, from 2.1.0 | 14 | 989/1000 | 0.000264 |

An earlier version of this page said the 0.0065 was real distance — that this isolate genuinely sat that
far from *B. suis* str. Thompsen. It was not. It was still mostly error k-mers, and the threshold was
re-derived because of it: see [The k-mer copy threshold](K-mer-threshold). The species and ST 15 were the
same at every one of the three, which is the point worth keeping: the identification was never in doubt,
but the *number beside it* measured how deeply the library happened to be sequenced rather than how far
the genome is.

Read a distance from this table the way the report means it — it is now a property of the sample. The
report says which threshold each sample got and the coverage it came from, and `--min-kmer-copies 2`
reproduces a 1.1.0 report exactly.

## 7. The quality assurance report
`results/brucellaid_report.html` is the run with its evidence attached: open it in a browser, or print it
&mdash; it lays out one sample to a page.
**[This is what it looks like](https://duceppemo.github.io/BrucellaID/example_report.html)** for the four samples of this tutorial, with the machine it ran on
removed.
The *B. canis* sample is the one to look at — its block lists the five closest references, with
`GCA_000530495.1_Oliveri_genomic.fna.gz` highlighted as the closest and
`GCA_000741995.1_ASM74199v1_genomic.fna.gz` shaded, because at
0.000119 against 0.000337 the distance did not really choose. Beside them, **The species call** shows
what decided it:

```
Species reported                          Brucella canis
Decided by                                MASH+ST
Species within reach of the closest       Brucella canis, Brucella suis
Sequence type                             20
What that ST is, in the isolate database  ST 20 is Brucella canis in 70 of 70 isolate(s)
```

The file fetches nothing from the network, so it can be emailed, archived or opened on an offline machine
and still look the same. `--pdf` writes [the same report as a PDF](https://duceppemo.github.io/BrucellaID/example_report.pdf) (it needs reportlab),
`--operator NAME` puts a name on it, and `--no-html` skips it.

`results/brucellaid_mqc.json` needs no options at all:

```bash
multiqc results/
```
gives a MultiQC report with a **BrucellaID** table of the four samples.

## 8. The JSON report
`results/brucellaid.json` holds everything, including what produced it:

```json
"software": {
  "brucellaid": "2.0.0",
  "Python": "3.12.16",
  "MASH": "2.3",
  "BBTools": "40.05",
  "Java": "openjdk version \"25.0.4.1\" 2026-08-18"
}
```
and, per sample, the species with its accessions, the close matches, every allele with its reads and its
runner-up, and the error if the sample failed.

## 9. The legacy report
If something downstream reads the xlsx of the first version:

```bash
brucellaid -i reads -o results --xlsx
```
writes `results/stat_alignment_summary_<date>.xlsx` next to the other reports, with the same four columns
as before and the reference table's code (`canis`, `abortus`, `melitensis`, `suis` with the default set;
vSNP's `ab1`, `mel3`, `suis2` with `--reference-set vsnp`) in the `species` column.

## Where to go next
* [Output files](Output-files) — every column of every report
* [How it works](How-it-works) — what the two calls rest on, and where each of them stops
* [Validation suite](Validation-suite) — the same pipeline against genomes whose answer is known, and
  reproducible from the repository
