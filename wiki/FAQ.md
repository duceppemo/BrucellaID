# FAQ

### Can it read assemblies?
No, and it refuses them rather than typing them badly. Both steps count reads: MASH sketches a read set and
drops the k-mers seen fewer times than a depth-derived copy threshold, and the allele calling works because
a locus is covered by many reads that
each vote. A contig votes once, and a contig holding several loci would be assigned to one of them. For an
assembly, [brucetyper](https://github.com/duceppemo/brucetyper) does *in silico* PCR, and `mash dist`
against `brucellaid/data/references.msh` identifies a genome directly.

### Nanopore or PacBio reads?
They work. The reads are long enough to cover a whole locus, which helps, but their error rate costs k-mers:
`hdist=0` means a read with an error in a 25-mer loses that k-mer for every allele at once, so the ranking
holds while the counts drop. Use reads that have been corrected or basecalled with a recent model, and check
the margins in `brucellaid_mlst.tsv`.

### How much coverage do I need?
About 10x is comfortable, 20x is plenty. Below 5x, loci start coming back as ties, which BrucellaID reports
rather than guesses. The species identification is the more forgiving of the two: at 5x the sketch still
shares most of its hashes with the right reference.

### What does `SpeciesSource` mean, and can the ST overrule MASH?
It says what the species rests on. `MASH` is the distance alone, `MASH+ST (agree)` is both saying the same
thing, and `(arbitrated)` means a reference of a second species was within 0.0005 of the closest one — so
MASH had not really chosen — and the sequence type chose among those species. In the sketch of the package
that happens for *B. canis* against *B. suis* bv. 4, and for no other pair.

The ST cannot overrule a species MASH picked cleanly. An ST pointing at a species that was nowhere near the
closest reference is reported as `MASH (conflict)`, with the MASH call kept: a thousand min-hashes against
nine loci. `--no-st-species` leaves the ST out of it entirely. See
[How it works](How-it-works).

### Why not just add more reference genomes instead?
More references fix a different problem than the one the ST fixes. *B. canis* and *B. suis* bv. 4 are
0.00024 apart, which is inside the variation within either species: in a dereplicated 428-genome *Brucella*
database the nearest pair of *canis* and *suis* genomes is 0.000194 apart. More references move that
boundary around; they do not create one, and the MLST loci, which do separate the two, are what settles it.
What more references *do* fix is coverage: the 16 genomes carry 8 species, and a sample of any other —
*B. pinnipedialis*, the whole former *Ochrobactrum* clade — can only land on a wrong neighbour or on `ND`.
That is what `--reference-set extended` is for.

### What is the extended reference set?
428 curated genomes covering 29 of the 30 named species of the *Brucella* clade (the 2020 taxonomy moved
*Ochrobactrum* into it) plus the lineages that have no name yet, selected with `--reference-set extended`.
The thirtieth is *B. papionis*, and no reference set can carry it: INSDC holds no assembly of that species
at all, only a single read set, so there is no genome to sketch.
[How it works](How-it-works#the-extended-reference-set) says where it comes from and
[Validation](Validation#at-scale-the-extended-reference-set) what it does on simulated reads of 2,022
assemblies and [on every public *Brucella* read set in ENA](Validation#on-real-public-read-sets-the-whole-archive).
The 16 vSNP genomes stay the default.

### Where did the biovar column go?
It was removed in 1.1, because it was never a biovar call and was read as one. It held the biovar of the
closest reference genome. Three things make a real biovar call impossible here: the sketch has references
for only *B. abortus* biovars 1 and 3, so 43% of biovar-labelled *B. abortus* isolates in PubMLST are of a
biovar it cannot name; the two *B. melitensis* biovar 1 references are further apart (0.000756) than the
biovar 2 and 3 references are from each other (0.000706), so the distance carries no biovar signal; and a
biovar is a phenotype that is not monophyletic — the best possible ST → biovar rule reaches 66% for
*B. abortus* and 63% for *B. melitensis*. See [How it works](How-it-works#no-biovar).

What is still there: `brucellaid_mash.tsv` gives the biovar of each *reference genome*, and the quality
assurance report shows the biovars recorded for the sample's ST with their counts. For a biovar, use
classical biotyping or [brucetyper](https://github.com/duceppemo/brucetyper).

### Why 150 shared hashes?
It is the threshold the first version of BrucellaID used, and it separates "this read set is *Brucella*"
from "the closest of 16 genomes happens to be this one" with room to spare: a real *Brucella* sample at 20x
shares 900 or more of the 1000 hashes, and a sample of the nearest neighbouring species shares almost none.
`--min-hashes` changes it.

### Can I update the MLST database?
Yes: `python3 scripts/make_mlst_data.py` downloads the scheme from PubMLST and rewrites
`brucellaid/data/mlst_alleles.fasta` and `mlst_profiles.tsv`. Without an API key PubMLST serves the data
submitted up to a cut-off date, which the script prints. `--alleles` and `--profiles` take files of your own
without touching the package.

### Can I add a reference genome?
Yes: `scripts/make_reference_data.py` rebuilds the sketch from the accessions of `references.tsv`, so add a
row and run it, or build a sketch yourself with `mash sketch` (MASH's defaults, k = 21 and 1000 hashes) and
pass `--sketch` and `--references`. BrucellaID warns when the sketch and the table do not name the same
genomes.

### What does `new` mean in the SequenceType column?
All nine alleles were called, but their combination is not in the profile table. See
[Troubleshooting](Troubleshooting).

### Why is the species column of the xlsx report `ab1` and not *Brucella abortus*?
Because that is what the first version wrote, and `--xlsx` exists for the pipelines that read it. The
species name is in `brucellaid.tsv` and in the JSON report, with the vSNP code as `species_code`.

### What is the HTML report for, and can I send it to someone?
It is the evidence behind a run, in one file: the closest reference genomes rather than the winning one
alone, every allele call with the margin it had, the checksums of the inputs and of the reference data,
the exact command, and a box for a reviewer to sign. It fetches nothing from the network — no script, no stylesheet, no font — so it
renders the same in an email attachment, on an offline machine, or out of an archive years later, and it
prints to clean pages. `--pdf` writes the same report as a PDF, which needs reportlab. `--no-html` skips it.

### Can I get the results into MultiQC?
Yes, with nothing to configure: every run writes `brucellaid_mqc.json`, so `multiqc results/` shows a
**BrucellaID** table with the species, what decided it, the sequence type and its species, the shared
hashes, the closest reference and the status of each sample.

### Does it detect mixed samples?
Not as such. A mixture shows up as a lower share of hashes, a close second reference, or ties in the allele
calls — all of which are reported — but BrucellaID will not tell you "two strains". Read
`brucellaid_mash.tsv` and the per-locus margins.

### Is it the same answer as vSNP?
For the species, almost: the same 16 reference genomes, the same sketch and the same threshold of 150 shared
hashes, so the MASH identification is vSNP's. It can differ for one pair — a *B. canis* or *B. suis* bv. 4
read set, where the sequence type now chooses between the two instead of the distance deciding by a margin
of noise. `--no-st-species` gives vSNP's answer exactly. For the MLST, the alleles are now called from the reads with Seal rather than with
MLSTsippr, and the database comes from PubMLST rather than from a copy that had drifted: vSNP's copy still
held a `gyrB_23` that PubMLST has since deleted, and that allele matched the *B. abortus* 9-941 genome
better than the real one, which turned ST 1 into `new`.

### What does it not do?
Assembly, SNP calling, phylogeny, AMR, biovar typing and MLVA. It names the species and the sequence type,
and it writes down what both answers rest on.
