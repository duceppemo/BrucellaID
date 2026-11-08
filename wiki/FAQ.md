# FAQ

### Can it read assemblies?
No, and it refuses them rather than typing them badly. Both steps count reads: MASH sketches a read set and
drops the k-mers seen fewer times than a depth-derived copy threshold, and the allele calling works because
a locus is covered by many reads that
each vote. A contig votes once, and a contig holding several loci would be assigned to one of them. For an
assembly, [brucetyper](https://github.com/duceppemo/brucetyper) does *in silico* PCR, and `mash dist`
against `brucellaid/data/references_extended.msh` identifies a genome directly.

### Nanopore or PacBio reads?
**Untested.** Nothing in the validation is a long read: the suite simulates 150 bp paired-end, and the
at-scale sweep is filtered to paired Illumina WGS runs, so there is no measurement behind an answer either
way. In principle they should work, and the reads are long enough to cover a whole locus, which helps — but
their error rate costs k-mers:
`hdist=0` means a read with an error in a 25-mer loses that k-mer for every allele at once, so the ranking
holds while the counts drop. Use reads that have been corrected or basecalled with a recent model, and check
the margins in `brucellaid_mlst.tsv`.

### How much coverage do I need?
20x is plenty and 10x is comfortable, and the two answers give out at different depths. The species
identification is the more forgiving: striding three real libraries down depth by depth, the sketch still
shares 834 to 950 of its 1000 hashes with the right reference at 5x to 7x, and 983 or more from 14x up — the
figures are on [The k-mer copy threshold](K-mer-threshold). The allele calls go first, as loci start coming
back as ties, which BrucellaID reports rather than guesses. How shallow that begins is not measured here:
the suite simulates 20x and the shallowest real library in it is 47x, so run
`COVERAGE=5 bash validation/run_validation.sh` to see where the ties start on data like yours.

### What does `SpeciesSource` mean, and can the ST overrule MASH?
It says what the species rests on. `MASH` is the distance alone, `MASH+ST (agree)` is both saying the same
thing, and `(arbitrated)` means a reference of a second species was within 0.0005 of the closest one — so
MASH had not really chosen — and the sequence type chose among those species. Comparing each packaged
genome with all the others, a second species falls inside the margin for 38 of the 432 default references
over 14 species pairs — *canis*/*suis*, *ceti*/*pinnipedialis* and the former *Ochrobactrum* ones among
them — and for 7 of vSNP's 16 over 8 pairs. The classical *canis*/*suis* case is the one that comes up in
practice, and it comes up the same way under either set.

`MASH (unresolved)` is the same situation as `(arbitrated)` with the arbitration missing: two species were
within reach and the ST could not choose, so the species is the closest reference's and the sample is
flagged for review. [Troubleshooting](Troubleshooting#speciessource-says-mash-unresolved) lists the seven
reasons an ST cannot be used.

The ST cannot overrule a species MASH picked cleanly. An ST pointing at a species that was nowhere near the
closest reference is reported as `MASH (conflict)`, with the MASH call kept: a thousand min-hashes against
nine loci. `--no-st-species` leaves the ST out of it entirely. See
[How it works](How-it-works#the-species-call).

### Why not just add more reference genomes instead?
More references fix a different problem than the one the ST fixes. *B. canis* and *B. suis* bv. 4 are
0.00024 apart, which is inside the variation within either species: in the 432-genome default set the
nearest pair of *canis* and *suis* genomes is 0.000264 apart, still the closest pair of different species
in the sketch. More references move that boundary around; they do not create one, and the MLST loci, which
do separate the two, are what settles it.
What more references *do* fix is coverage, which is why the extended set is the default from 2.0: vSNP's
16 genomes carry 8 species, and a sample of any other —
*B. pinnipedialis*, the whole former *Ochrobactrum* clade — can only land on a wrong neighbour or on `ND`.
`--reference-set vsnp` is still there for the runs that want the old 16, and nothing else.

### What is the extended reference set?
432 curated genomes covering all 30 named species of the *Brucella* clade (the 2020 taxonomy moved
*Ochrobactrum* into it) plus the lineages that have no name yet. It is what a run uses from 2.0.
One of the 30 took some finding. Nothing in INSDC is labelled *B. papionis*, so it looked as though no
set could carry it — but the genome is there under *Brucella* sp., and
[`curation/papionis/README.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/papionis/README.md)
records how it was identified and relabelled.
[How it works](How-it-works#the-reference-sets) says where it comes from and
[Validation](Validation-at-scale#simulated-reads-of-2022-assemblies) what it does on simulated reads of 2,022
assemblies and [on every public *Brucella* read set in ENA](Validation-at-scale#every-usable-public-read-set-in-ena).
The 16 vSNP genomes are still a run away, as `--reference-set vsnp`.

### Where did the biovar column go?
It was removed in 1.1, because it was never a biovar call and was read as one. It held the biovar of the
closest reference genome. Three things make a real biovar call impossible here: hardly any reference
carries a biovar at all, 7 of the 432 genomes of the default set, so 93% of the *B. abortus* isolates in
PubMLST that have a numbered biovar are of a biovar none of them can name; the two *B. suis* biovar 2
genomes are further apart (0.00163) than the biovar 3 genome is from one of them (0.00150), so the
distance carries no biovar signal; and a
biovar is a phenotype that is not monophyletic — the best possible ST → biovar rule reaches 66% for
*B. abortus* and 63% for *B. melitensis*. See [How it works](How-it-works#no-biovar).

What is still there: `brucellaid_mash.tsv` gives the biovar of each *reference genome*, and the quality
assurance report shows the biovars recorded for the sample's ST with their counts. For a biovar, use
classical biotyping or [brucetyper](https://github.com/duceppemo/brucetyper).

### What does "loose match" in the notes mean?
That the closest reference is further than `--loose-distance` (0.001 by default) — the sample matched
something, but not closely. It is **not** a claim that the species is wrong, and BrucellaID does not change
the call: over every public *Brucella* read set in ENA, agreement with the submitter barely moves with
distance — 99.8% of calls within 0.0005 agree, and so do 95% to 100% of those in every band beyond it. The
note flags about 4% of identified runs, and 98% of those still agree with their submitter, so it marks a
library worth a second look rather than a species worth doubting;
[the derivation](Validation-at-scale#the-loose-match-distance) has the figures. The note names three
causes: the reads
hold more than *Brucella*, which is the usual one; the sample is a lineage with no near relative in the
reference set; or a minority *Brucella* population is being filtered out along with the sequencing errors,
which `--min-kmer-copies` either side of the derived value will show — see
[The k-mer copy threshold](K-mer-threshold). Read the note's distance and
shared hashes, and `brucellaid_mash.tsv` for the runners-up.

### Why 150 shared hashes?
It is the threshold the first version of BrucellaID used, and it separates "this read set is *Brucella*"
from "the closest genome of the set happens to be this one" with room to spare: a real *Brucella* sample at 20x
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

### Why is the species column of the xlsx report a code and not *Brucella abortus*?
Because that is what the first version wrote, and `--xlsx` exists for the pipelines that read it. The
species name is in `brucellaid.tsv` and in the JSON report, where the same code is `species_code`.
Which code depends on the reference set: the default set spells the species epithet (`abortus`), and
`--reference-set vsnp` gives vSNP's own `ab1`, which is what 1.x wrote. The run warns when it is writing
epithets into that column — see [Output files](Output-files).

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
For the species, `--reference-set vsnp` is almost exactly vSNP's answer: the same 16 reference genomes,
the same sketch and the same threshold of 150 shared
hashes, so the MASH identification is vSNP's. It can differ for one pair — a *B. canis* or *B. suis* bv. 4
read set, where the sequence type now chooses between the two instead of the distance deciding by a margin
of noise. `--no-st-species` gives vSNP's answer exactly. For the MLST, the alleles are now called from the reads with Seal rather than with
MLSTsippr, and the database comes from PubMLST rather than from a copy that had drifted: vSNP's copy still
held a `gyrB_23` that PubMLST has since deleted, and that allele matched the *B. abortus* 9-941 genome
better than the real one, which turned ST 1 into `new`.

### What does it not do?
Assembly, SNP calling, phylogeny, AMR, biovar typing and MLVA. It names the species and the sequence type,
and it writes down what both answers rest on.
