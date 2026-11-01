# Tutorial

Four public *Brucella* read sets, from downloading them to reading the reports: about ten minutes, most of
it the download. Every command and every output below is a real run.

## 1. Install
```bash
conda create -n brucellaid -c conda-forge -c bioconda brucellaid
conda activate brucellaid
brucellaid --version
```
```
brucellaid 1.0.0
```

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
16 reference genome(s) in .../brucellaid/data/references.msh
9 MLST locus/loci: gap, aroA, glk, dnaK, gyrB, trpE, cobQ, int_hyp, omp25
118 sequence type(s) with a species, to decide the species MASH cannot
4 sample(s) found in reads
4 sample(s), 2 at a time, 4 thread(s) each
ERR2136546: identifying the species with MASH (2 file(s))
ERR2136546: Brucella canis str. ATCC 23365 (NC_010103-NC_010104.fasta, MASH distance 0.00652, 773/1000 hashes)
ERR2136546: calling the MLST alleles with Seal
ERR2136546: ST 20 (profile 1-6-4-1-5-3-5-4-6)
ERR2136546: Brucella canis (MASH+ST, arbitrated)
ERR3454415: Brucella abortus bv. 1 str. 9-941 (NC_006932-NC_006933.fasta, MASH distance 0.00455, 833/1000 hashes)
ERR3454415: ST 1 (profile 2-1-1-2-1-3-1-1-1)
SRR19520329: Brucella melitensis bv. 3 str. Ether (NZ_CP007760-NZ_CP007761.fasta, MASH distance 0.00329, 875/1000 hashes)
SRR19520329: ST 11 (profile 3-2-3-2-1-5-3-2-10)
SRR8550484: Brucella suis bv. 2 str. ATCC 23445 (NC_010169-NC_010167.fasta, MASH distance 0.02484, 422/1000 hashes)
SRR8550484: ST 15 (profile 1-2-7-1-3-3-5-3-2)
Done in 30s
```
30 seconds on 8 threads, for 4.6 million read pairs. `-p 2` runs two samples at a time and gives each of
them four of the eight threads.

## 4. Read the main table
`results/brucellaid.tsv`, also printed to the terminal:

```
Sample       Species             Biovar  SpeciesSource         ClosestReference               SharedHashes  SequenceType  AlleleProfile       STSpecies            Reads
ERR2136546   Brucella canis              MASH+ST (arbitrated)  NC_010103-NC_010104.fasta      773/1000      20            1-6-4-1-5-3-5-4-6   Brucella canis       625576
ERR3454415   Brucella abortus    1       MASH+ST (agree)       NC_006932-NC_006933.fasta      833/1000      1             2-1-1-2-1-3-1-1-1   Brucella abortus     1175082
SRR19520329  Brucella melitensis 3       MASH+ST (agree)       NZ_CP007760-NZ_CP007761.fasta  875/1000      11            3-2-3-2-1-5-3-2-10  Brucella melitensis  1136924
SRR8550484   Brucella suis       2       MASH+ST (agree)       NC_010169-NC_010167.fasta      422/1000      15            1-2-7-1-3-3-5-3-2   Brucella suis        1693672
```
Four species, four sequence types, and each species is the one the submitter labelled the run with.
`MashDistance`, left out above to fit, reads 0.006524, 0.004549, 0.003285 and 0.024841.

**The first sample is the interesting one.** Its `SpeciesSource` says `arbitrated`, and the `Notes` column,
also left out above, says why:

```
nearly as close: B-REF-BS4-40.fasta (Brucella suis bv. 4, 768/1000 hashes); MASH cannot separate
Brucella canis and Brucella suis: ST 20 is Brucella canis in 70 of 70 isolate(s), which agrees with
the closest reference
```
*B. canis* came first by **five hashes out of a thousand** — 0.006524 against 0.006698, a margin that is
noise, not a measurement. So the species was not decided by that margin: ST 20 was, and ST 20 is *B. canis*
in all 70 isolates of the PubMLST isolate database that carry it, and never *B. suis*. Had the distances
come out the other way round, as they easily could have, the answer would have been the same. See
[How it works](How-it-works).

## 5. Read what the calls rest on
`results/brucellaid_mash.tsv` has every distance. For the *B. canis* sample, the first four:

```
Sample      Reference                      Species           Biovar  MashDistance  SharedHashes
ERR2136546  NC_010103-NC_010104.fasta      Brucella canis            0.006524      773
ERR2136546  B-REF-BS4-40.fasta             Brucella suis     4       0.006698      768
ERR2136546  NC_017251-NC_017250.fasta      Brucella suis     1       0.007196      754
ERR2136546  NZ_CP007719-NZ_CP007718.fasta  Brucella suis     3       0.007268      752
```

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

## 6. One sample that is worth a second look
`SRR8550484` is identified as *B. suis* and typed as ST 15, with all nine loci called and every one of them
ahead of its runner-up by a comfortable margin. But its MASH distance is 0.0248 where the other three are
between 0.0033 and 0.0065, and it shares 422 of 1000 hashes where they share 773 to 875. Something in that read set is not *Brucella*, or not clean: at 1.7 million
read pairs, sequencing errors are deep enough to survive MASH's `-m 2` filter and fill the sketch with
k-mers no reference has.

It is still typed, and the MLST is unbothered — the allele calls only use the reads that match an allele —
but it is the kind of number worth carrying into the report. BrucellaID prints it rather than rounding it
away.

## 7. The quality assurance report
`results/brucellaid_report.html` is the run with its evidence attached: open it in a browser, or print it.
The *B. canis* sample is the one to look at — its block holds all sixteen distances, with
`NC_010103-NC_010104.fasta` highlighted as the closest and `B-REF-BS4-40.fasta` shaded, because at
0.006524 against 0.006698 the distance did not really choose. Under them, **The species call** shows what
decided it:

```
Species reported                          Brucella canis
Decided by                                MASH+ST
Species within reach of the closest       Brucella canis, Brucella suis
Sequence type                             20
What that ST is, in the isolate database  ST 20 is Brucella canis in 70 of 70 isolate(s)
```

The file fetches nothing from the network, so it can be emailed, archived or opened on an offline machine
and still look the same. `--pdf` writes the same report as a PDF (it needs reportlab), `--operator NAME`
puts a name on it, and `--no-html` skips it.

`results/brucellaid_mqc.json` needs no options at all:

```bash
multiqc results/
```
gives a MultiQC report with a **BrucellaID** table of the four samples.

## 8. The JSON report
`results/brucellaid.json` holds everything, including what produced it:

```json
"software": {
  "brucellaid": "1.0.0",
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
as before and the vSNP species code (`canis`, `ab1`, `mel3`, `suis2`) in the `species` column.

## Where to go next
* [Output files](Output-files) — every column of every report
* [How it works](How-it-works) — what the two calls rest on, and where each of them stops
* [Validation](Validation) — the same pipeline against genomes whose answer is known
