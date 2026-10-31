# How it works

A run does two independent things to the same reads. Neither needs the other, and neither corrects the
other: when they disagree, both numbers are in the reports and the disagreement is the finding.

## 1. Which species: the MASH distance to 16 reference genomes

The reads of a sample — R1 and R2 together, as one read set — are piped into
`mash dist -r -m 2 references.msh -`. MASH reduces the read set to the 1000 smallest hashes of its 21-mers
and compares that sketch to the sketch of each reference genome, which gives a distance and the number of
hashes the two share.

* `-m 2` drops the k-mers seen once in the read set. Most of them are sequencing errors, and they would
  otherwise fill the sketch with noise.
* `-r` tells MASH the input is reads, so the genome size is estimated from k-mer multiplicity.
* The k-mer size and the sketch size are inherited from the reference sketch, which is why a custom sketch
  has to keep MASH's defaults.

The closest reference names the species, and the reference table also gives its biovar, strain and NCBI
accessions. A reference is only reported when the reads share at least **150 of its 1000 hashes**
(`--min-hashes`). Below that the sample is reported as `ND`: there is some *Brucella* in it, or there is
not, but the closest of 16 genomes is not an answer. The first version of BrucellaID used the same
threshold.

**What this tells you, and what it does not.** The *Brucella* species are extraordinarily close — they are
over 99% identical — so a MASH distance separates them by a small margin, and the biovar references by a
smaller one still. Read the species as the answer and **the biovar as the nearest reference, not as a
biovar call**: `brucellaid_mash.tsv` lists every distance, and a note on the sample says when another
reference was nearly as close. For a biovar, use the classical biotyping, or
[brucetyper](https://github.com/duceppemo/brucetyper), which reads the Bruce-ladder and Suis-ladder assays
*in silico*.

Of the 16 references, two (*B. ceti* Cudo and *B. suis* bv. 4 strain 40) come from vSNP under their own
names rather than from an accession; the other 14 are NCBI genomes, listed in `references.tsv`.

## 2. Which sequence type: the nine MLST alleles, called with Seal

The same reads go to `seal.sh` with the 309 alleles of the nine loci as references, and these settings:

```
k=25 hdist=0 rcomp=t maskmiddle=f ambiguous=all clearzone=0
```

Seal counts, for each read, which reference sequences share the **highest number of 25-mers** with it, and
counts the read for each of those and for no other. The consequence is what the call rests on:

* the allele the sample really carries shares every k-mer of every read of that locus, so it is counted for
  **all** of them;
* another allele of the same locus is only counted for the reads that fall where the two agree, and loses
  every read that covers a position where they differ.

So for each locus the allele with the highest count wins, and `brucellaid_mlst.tsv` prints the runner-up and
its count next to it, which is the margin. Two alleles with the same count are reported as a tie
(`ND` for that locus) rather than broken by a coin toss: at a low depth, two alleles can agree everywhere
the reads happen to cover. A locus with fewer than `--min-reads` reads (5) is not called either.

`hdist=0` is deliberate. The alleles of a locus differ by a few bases — `trpE_3` and `trpE_5` differ at one
position of 486 — and a k-mer allowed a mismatch would match both.

The nine allele numbers are then looked up in the profile table:

| Result | What it means |
| --- | --- |
| A number | That profile is in the table: the sample is that ST |
| `new` | All nine alleles were called, but the combination is not in the table |
| `ND` | At least one locus has no call. The `Notes` column says which, and why |

**Novel alleles.** An allele that is not in the database is called as its nearest neighbour, which usually
shows up as a lower count and a close runner-up, and often as a `new` profile. BrucellaID does not assemble,
so it cannot tell you the new sequence: assemble the sample and submit the locus to
[PubMLST](https://pubmlst.org/organisms/brucella-spp) if you suspect one.

**Depth.** The alleles are 422 to 589 bp long. At 20x, a locus is covered by roughly 80 reads and the margin
over the runner-up is comfortable; below about 5x the ties start, and they are reported as ties rather than
guessed. Coverage at the very ends of a locus matters as much as in the middle: two alleles differing at a
single position near one end can only be told apart by a read that reaches it.

## What the reference data is

| File | Source |
| --- | --- |
| `references.msh`, `references.tsv` | The *Brucella* part of [vSNP](https://github.com/USDA-VS/vsnp)'s reference set, 16 genomes |
| `mlst_alleles.fasta`, `mlst_profiles.tsv` | The MLST (9 locus) scheme of [PubMLST](https://pubmlst.org/organisms/brucella-spp), downloaded with `scripts/make_mlst_data.py` |

The MLST scheme itself is [Whatmore, Perrett and MacMillan (2007)](https://doi.org/10.1186/1471-2180-7-34):
nine loci, 4,396 bp in all.

## Why not an assembly?

Both steps count reads. MASH sketches a read set and filters the k-mers seen once, and the allele calling
works because a locus is covered by many reads that each vote for an allele; a contig votes once, and a
contig holding several loci would be assigned to one of them. A fasta file is refused rather than quietly
mis-typed. For assemblies, [brucetyper](https://github.com/duceppemo/brucetyper) does *in silico* PCR
instead.
