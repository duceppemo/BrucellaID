# Changelog

## Unreleased

**One genome of the extended set changes label: `GCA_900095155.1` is now *Brucella* sp., not
*B. inopinata*.** Strain 141012304 is 97.87% ANI and 0.0212 from the *B. inopinata* type strain BO1, which
is further than three references the set already calls *Brucella* sp. (BO2 at 98.24%, 2280 at 98.08%, BO3
at 98.01%). The other *inopinata* genome beside the type, `GCA_030270765.1`, is 99.88% from it. The set
was calling one early-diverging lineage *inopinata* and three nearer ones unnamed, and now does not.

A sample that matches this genome is reported as *Brucella* sp. where it was reported as *B. inopinata*.
That is one of the 4,024 public read sets, `SRR12224742`, whose submitter named no species, and none of
the 2,022 assemblies but the genome itself; its distance, shared hashes and sequence type are unchanged,
and so is every figure on [Validation at scale](Validation-at-scale). The sketch is not rebuilt — the
label is a row of `references_extended.tsv`.

**Four genomes of the extended set take a name outside the genus: *Ochrobactrum teleogrylli* and
*Ochrobactrum chromiisoli*.** They were *Brucella* sp. Both species were described after *Ochrobactrum*
was moved into *Brucella*, neither has a *Brucella* combination, and NCBI files both outside taxon 234,
so the collection the set was curated from held some of their genomes and none of their type strains.
By ANI to those type strains:

* `GCA_030063905.1` and `GCA_040372265.1` are *O. teleogrylli*: 99.0–99.2% from the *O. soli* type
  strain and 98.2% from the *O. teleogrylli* one. The two types are 98.16% apart, so
  one species, and *teleogrylli* (2020) has priority over *soli* (validated 2024).
* `GCA_039674265.1` and `GCA_039674445.1` are *O. chromiisoli*: 98.5% from its type strain.

**A sample can now be reported under a species name that does not begin with *Brucella*.** Anything that
parses the `Species` column on that assumption should be checked. Three of the 4,024 public read sets
change: `SRR12605656`, `SRR25828924` and `SRR25828932`, deposited as *B. haematophila* and *B. anthropi*,
go from *Brucella* sp. to *O. teleogrylli* with every measurement unchanged, and their assemblies are
83–90% ANI from the type strains of the species on their labels. No run matches a *chromiisoli* genome.
Removing the species from the set was considered and rejected: without its genomes those libraries share
93 to 95 of 1000 hashes with their best match and would be reported as not *Brucella*.

**The extended set is 429 genomes, not 432: three partial genomes are removed and the sketch rebuilt.**
`GCA_040391195.1` is 2.1 Mb of *O. teleogrylli*, and `GCA_032085995.1` and `GCA_032104425.1` are 3.2 and
3.7 Mb of *B. pseudogrignonensis*, all species of about 5 Mb that the set holds complete genomes of. A
reference holding part of its genome cannot be matched closely by a complete library of its own species,
and none was: not one of the three is the best match of any identified run among the 4,024 public read
sets, or within the close-match margin of one. The other 429 genomes are in the rebuilt sketch with
exactly the hashes they had, so every distance to them is what it was. A custom `--references` table
written against the 432 will name three genomes the sketch no longer holds, which BrucellaID reports at
start-up.

**Eleven more unnamed genomes of the set take the name of the *Brucella* species they belong to.** They
were *Brucella* sp. and are 96.2–99.3% ANI from a type strain, measured in both directions: one
*B. tritici*, three *B. anthropi*, two *B. intermedia*, one *B. haematophila* and four
*B. pseudogrignonensis*. Four of them take the senior name of a group ANI cannot separate, not the name
of their nearest type strain, which is the rule the set already used. A sample matching one of the
eleven was reported as *Brucella* sp. and is now reported under its species. No public read set matches
any of them, so no call in the sweep changes. `GCA_054871705.1` was left unnamed: at 95.2% from the
*B. pecoris* type strain it is too near the boundary to call.

**Arbitration reaches fewer references: a second species is within the close-match margin of 33 of the
429, across 10 species pairs, where it was 38 across 14.** Four of the fourteen pairs were a named species
against one of its own unnamed genomes, so a sample between them would have been reported as
unresolved between *B. anthropi* and *Brucella* sp. when both genomes were *B. anthropi*. The densest case is still 15.

With all of this the set is 30 named *Brucella* species, the two *Ochrobactrum* ones in four genomes, and
31 unnamed genomes. No held-out named sample of the simulated-read validation moves, and no count on
[Validation at scale](Validation-at-scale) that `scripts/score_public_sweep.py` derives changes.

**One sentence on that page was stale and is corrected.** It gave the 95 unnamed simulated samples as 53
reported *Brucella* sp. and 40 named, 36 of them with a note. The full run against the 432-genome set
gave 56 and 37, which is the old figure plus the three samples that moved when two unnamed lineages were
restored to the set, so it had probably not been recomputed since. After all the relabels above it is
**37 and 56**, ten of the 56 with a note.

All of this came out of chasing the 43 public read sets reported as *Brucella* sp., written up in
`curation/unnamed_lineages.md`. That note also records what was left alone: 22 of those runs are the
vampire-bat lineage published as "*B. nosferati*", a name that is not validly published, and it stays
*Brucella* sp.

**The at-scale validation is re-measured under the shipped k-mer copy threshold.** Every figure on
[Validation at scale](Validation-at-scale) came from the sweep that ran before 2.1.0 re-derived the
threshold, and was carried with a caveat saying so. All 4,024 usable public ENA read sets were swept again
against both reference sets, with MLST, and the page now states measurements rather than caveats.

**The headline figures are unchanged**: the extended set matches the submitter on 98.7% and the vSNP set
on 96.2%, the same percentages as before. Not one of the extended set's 52 non-matches is a case where it
held the right reference and named the wrong species.

**`--loose-distance 0.001` is confirmed.** It had been set from a weighted 207-run sample that put the
95th percentile of agreeing calls at 0.00066; over all 4,001 identified runs it is **0.00068**, so the
constant needs no change. `scripts/loose_distance_calibration.tsv` now holds every identified run instead
of the sample, so the percentiles carry no sampling error.

**One claim is withdrawn.** The sample supported "agreement with the submitter is flat in distance", and
the full sweep does not: it falls from 99.5% on close calls to **81.2%** between 0.005 and 0.02. A loose
match does carry some signal that the species is wrong. The signal is weak — the worst band is still
mostly right, and part of the fall is the submitters' own mislabels, since the thin-reference former
*Ochrobactrum* runs cluster there — so the note is still reported rather than acted on, which is what the
claim was supporting. Corrected in `brucellaid/mash.py`, [How it works](How-it-works), the
[FAQ](FAQ) and on the page.

* Head-to-head figures re-measured: the two sets disagree on **168** runs rather than 146, the extended set
  is right on **95** of them and the vSNP set on none, and **89** are runs the vSNP set cannot name at all.
* **23 runs share too little with the clade to be named**, holding 0 to 8 of 1000 hashes, where the
  pre-2.1 sweep saw 22 runs spanning 0 to 146. The difference is the re-derived threshold removing a trace
  population from deep mixed libraries, which is what the one borderline case at 146 hashes turned out to
  be.
* **Eighteen runs publish three fastq files**, which is how ENA presents a run holding both paired and
  unpaired reads: `_1` is the unpaired ones and `_2`/`_3` are the mates. The sweep's staging paired `_1`
  with `_2`, feeding the orphan file in as R1, so Seal refused the pair and all eighteen came back with no
  sequence type. Re-run on `_2`/`_3` they type cleanly — seventeen *B. melitensis* ST 8 and one
  *B. microti* ST 53, each agreeing with the species MASH called, whose distances moved by at most one
  shared hash. `find_samples()` now warns when a `_3` file sits beside a pair it has formed, since that is
  the only point at which the layout is visible.
* **A false `new` sequence type, found by the independent check.** `SRR27874758` is a mixture — 230 of
  1000 shared hashes, distance 0.047, 7.8 Mb of assembly against 3.3 Mb of genome — and its reads give a
  profile that is not in the table, reported as `new`. The type read off its own assembly is ST 8, a
  common *B. melitensis* type agreeing with the species call. Every locus had a runner-up 20 to 54 reads
  behind a winner with hundreds, which is what a mixture looks like and the opposite of what a novel
  allele looks like. Written up under
  [The sequence type is `new`](Troubleshooting#the-sequence-type-is-new); no code changed, since the
  sample was already flagged by the loose-match note and its shared-hash count.
* 20 of the 4,044 runs were not analysed: 18 never finished downloading a second mate, and two were
  already unavailable. A run with one mate is left out rather than analysed half-depth, since the copy
  threshold is derived from coverage.

## 2.1.1 (2026-11-08)

**The reports no longer show the biovars recorded for a sample's sequence type.** The quality assurance
report carried a `Biovars of its isolates` row in each sample's evidence block, holding the biovars PubMLST
records for that ST with their counts. It is gone from both the HTML and the PDF.

The reason is the shape of the data rather than the idea. Of the 118 sequence types in the shipped table,
**42 have a single biovar on record** and 59 have none, so on about a third of samples the row rendered one
bare value — `2:60` — beside that sample's own species and sequence type. That reads as a biovar call, which
is precisely how the reference-genome biovar column came to be misread before 1.1.0 removed it. The 17 STs
that record several biovars did convey the ambiguity well; they were the minority.

Nothing measured is lost. The counts remain in the packaged `mlst_st_species.tsv`, in a custom
`--st-species` table, and in `brucellaid.json`; the reference genomes' own biovars remain in
`brucellaid_mash.tsv`, where they are facts about the references. The report's methods section still sets
out why no biovar can be measured — 7 of the 432 references carry one, the two *B. suis* biovar 2 genomes
are further apart than the biovar 3 genome is from one of them, and the best any ST-to-biovar rule manages
is 66% — and now says plainly that this is why it shows none at all.

* The wiki's biovar explanation is cut to one short treatment under
  [How it works](How-it-works#no-biovar), where it had grown to seven places and about forty lines of
  argument for something the tool does not do. The pages that mentioned it now point there in a clause.
* **393 tests**, 100% statement and branch coverage. One is new, and asserts the row is absent from both
  reports while the fixture still carries biovars, so it cannot come back unnoticed.

## 2.1.0 (2026-11-08)

**The k-mer copy threshold was re-derived, and reported distances change because of it.** Before MASH
sketches a read set it drops every k-mer seen fewer than *m* times, to throw away the ones that exist only
because of a sequencing error. *m* was one more copy per 20x of coverage, capped at 20; it is now one more
copy per **5x**, capped at **48**. The floor of 2 is unchanged.

No species call and no sequence type changes: the 2,022-sample simulated validation returns exactly the
same 1,538 of 1,543 held-out named samples, with the same five mislabelled genomes, and the repository's
own suite still passes 22 of 22. What changes is `MashDistance` and `SharedHashes`, on real libraries
especially, because the old threshold was leaving error k-mers in the sketch:

| library | depth | before | after |
| --- | --- | --- | --- |
| ERR2136546 | 47x | 0.00652, 773/1000 | **0.000119, 995/1000** |
| ERR3454415 | 70x | 0.00068, 972/1000 | **0.000144, 994/1000** |
| SRR19520329 | 53x | 0.00051, 979/1000 | **0.000119, 995/1000** |
| SRR8550484 | 69x | 0.00645, 775/1000 | **0.000264, 989/1000** |

Simulated samples of a genome the sketch holds now return 999 or 1000 of 1000 hashes, which is the correct
answer for reads simulated from that genome.

**If you compare distances across versions, they are not comparable.** A distance from 2.0.0 or earlier is
larger than the same measurement today, by a factor that depends on the depth of the library — roughly 5
at 50x and 25 at 70x on the libraries above. `--min-kmer-copies` pins any fixed value if you need to
reproduce an old report: 2 gives a 1.1.0 report exactly.

**What it was fitted on.** 81 threshold sweeps over 63 real libraries from 5x to 1903x: three tutorial
libraries strided down to each depth, which separates depth from library quality, and 60 public ENA runs
stratified across the range, because the tutorial libraries stop at 70x and the median public *Brucella*
run is 141x. The new rule lands inside the measured plateau for 80 of the 81 curves; the old one for 26.
The measurements are in `scripts/kmer_threshold_sweep.tsv` and the derivation is on the wiki under
[The k-mer copy threshold](K-mer-threshold).

* **`--loose-distance` is recalibrated, 0.005 to 0.001.** The threshold is, as before, the 95th
  percentile of the MASH distance of calls that agree with the ENA submitter, taken to a round number —
  but that percentile was measured on distances the old copy threshold had inflated. Re-measured on 207
  public read sets under the new one it is 0.00066 rather than 0.00464, so the threshold follows it down.
  The note stays about as rare as it was, flagging 4.0% of the identified runs against 4.9% before, and
  98% of what it flags still agreeing with its submitter. (Both from the 207-run sample. The full sweep,
  recorded under Unreleased, gives 3.4% and 92%.) Left at 0.005 it would have missed the
  thin-reference cases it exists to catch, a *B. pituitosa* run at 0.0029 among them.
  `scripts/loose_distance_calibration.tsv` holds the measurements and
  `scripts/recalibrate_loose_distance.py` derives the number from them. (Both files were replaced when
  the sweep was re-run under this threshold: the calibration file now holds every identified run rather
  than a weighted sample, and the script no longer weights anything. The sample had put the percentile at
  0.00066 against the 0.00068 the full sweep gives, which is the evidence that the sampling was sound.)
* **A claim behind the old calibration was circular, and is withdrawn.** The derivation offered, as
  corroboration, that every agreeing call past the threshold shared fewer than 900 of its 1000 hashes.
  MASH derives the distance from the shared-hash count — over the 4,042 runs of the sweep, 317 distinct
  shared-hash values give 317 distinct distances and not one gives two — so that is the threshold restated
  in another unit, not support for it. What does support reporting the note rather than acting on it is how
  often a call at a given distance is still the submitter's species. Both tables are on
  [Validation at scale](Validation-at-scale#the-loose-match-distance). (This entry first said that
  agreement was flat in distance. It is not: measured over all 4,001 identified runs it falls from 99.5%
  on close calls to 81.2% between 0.005 and 0.02. The flat reading came from a 207-run sample. The
  conclusion is unchanged, because the worst band is still mostly right and some of the fall is the
  submitters' own mislabels, but the signal is not absent.)
* **Two of the four tutorial libraries no longer carry a loose-match note.** They were flagged only
  because their distance was inflated past `--loose-distance`; at 0.00012 and 0.00026 they are nowhere
  near even the new threshold.
* **The loose-match note gives different advice**, because its old advice no longer fits the rule. It used
  to suggest `--min-kmer-copies` *above* the derived value, which made sense when the derived value was
  three times too low. Now that it sits on the measured plateau, the remaining cause of a loose match on a
  deep library is a minority *Brucella* population being filtered out along with the errors — and for that
  the threshold has to come *down*. `SRR24420453` is the worked example: at the derived 48 it shares 160 of
  1000 hashes, at 20 it shares 396. The note now says to try either side of the derived value.
* **Three head-to-head figures of the public sweep are corrected**, and the evidence for them is now
  shipped. They had been computed before the sweep's retry pass filled in its last runs and were never
  recomputed: the extended set wins 94 of the 146 disagreements rather than 91, 89 of them are runs the
  vSNP set cannot name rather than 72, and 21 of those are *B. intermedia* rather than 15. Every figure
  moved in the extended set's favour. `scripts/public_sweep_disagreements.tsv` lists all 146 and
  `scripts/score_public_sweep.py` derives the figures from it.
* **The simulated-reads comparison is shipped too**, as `scripts/kmer_threshold_simulated.tsv`, and its
  claim is corrected. It read "their plateau runs 3-20 at every depth above 50x"; that sweep stopped at
  20 copies, so the upper edge was the edge of the measurement rather than a result. The uncensored point
  is the lower edge, which stays at 3 or 4 from 50x to 400x where a real library of 420-535x needs 32.
* **392 tests**, 100% statement and branch coverage.
* `validation/check_results.py` restates the rule in its own arithmetic rather than importing it, so that
  a drift in the shipped formula fails the validation instead of being matched silently. It duly failed
  when the rule changed, and was updated by hand.
* **The validation pages are split by topic**: [Validation](Validation) is now a short hub over
  [Validation suite](Validation-suite), [Validation at scale](Validation-at-scale) and
  [The k-mer copy threshold](K-mer-threshold).

## 2.0.0 (2026-11-07)

**The extended reference set is now the default, and that is a breaking change.** A run with no
`--reference-set` identifies against 432 curated genomes covering all 30 named species of the *Brucella*
clade instead of vSNP's 16 genomes and 8 species. Over every usable public *Brucella* read set in ENA —
4,042 runs — the two sets disagree on 146, and the extended set is right on 94 of them and the 16 on
none; across the 2,018 *B. melitensis* and 1,118 *B. abortus* runs that carry most real traffic they agree,
so the change adds coverage without disturbing the calls the smaller set already made well.

What this means in practice:

* **`--reference-set` takes `extended` (the default) or `vsnp`.** The option is new in 2.0 — 1.2.1 had no
  `--reference-set` at all — and the development builds that carried it briefly used `default` for the 16
  genomes. That value is gone: passing it gets an argparse error naming the two valid choices and exit 2,
  rather than silently changing behaviour.
* `--reference-set vsnp` gives back 1.x's reference data exactly: the same 16 genomes, the same sketch,
  the same threshold, and therefore the same species, closest reference, distance, shared hashes and
  sequence type. (**The distance and shared-hash half of that is no longer true from 2.1.0**, which
  re-derived the k-mer copy threshold: the species, closest reference and sequence type still match, but
  the numbers behind them moved. See the 2.1.0 entry.) It is not a byte-for-byte 1.x run, because the 2.0 features that are independent of the
  reference set still apply: a library past `--loose-distance` gains a loose-match note, which is inserted
  first and so shifts the other notes down, and `Loose match beyond` is added to the run parameters.
* **The `code` column changes with the set**, which matters for the legacy `--xlsx` report and for
  `species_code` in the JSON. With the extended set the codes are species epithets (`abortus5`,
  `anthropi`). **The break is partial, which is worse than total:** 7 of vSNP's 16 codes — `canis`,
  `microti`, `neotomae`, `ovis`, `suis1`, `suis2`, `suis3` — happen to spell the same thing in both
  tables and keep matching, while the other 9 — `ab1`, `ab3`, `mel1`, `mel1b`, `mel2`, `mel3`, `ceti1`,
  `ceti2`, `suis4` — silently stop. A code-to-species table will half work, per species. A pipeline
  reading those codes needs `--reference-set vsnp`, and `--xlsx` now warns when it is writing epithets
  instead, because that report exists for the pipelines that read the old codes and nothing else.
* **Arbitration now spans more species pairs, but is not triggered more often.** Comparing each packaged
  genome with all the others, a second species falls inside the close-match margin for 38 of the 432 and
  across 14 species pairs, against 7 of 16 and 8 pairs for vSNP's set — the new pairs are mostly in the
  former *Ochrobactrum* clade. The proportion is lower, not higher, because in a dense set a sample's
  nearest neighbours are usually other genomes of its own species, and on the four tutorial read sets the
  candidates and the verdict are identical under both sets. `--no-mlst` and `--no-st-species` switch the
  arbitrator off and report such a sample `MASH (unresolved)` with a warning, as they always did.
* **The distance no longer lists references it cannot separate.** Past 40 references inside the
  close-match margin the margin is measuring nothing: mash derives the distance from the shared-hash
  count, so at a low count every reference sharing the same few hashes gets exactly the same distance and
  ties with the best. A sample sharing 4 of 1000 hashes was seen to tie 127 references across twelve
  species, which put 10 kB of genome names in one `Notes` cell and 128 rows in the report. The note now
  gives the count and says the distance cannot separate them; the densest real case over the 432 packaged
  genomes is 15, so nothing legitimate is suppressed. The species call says so too, rather than claiming
  an agreement the distance never expressed: where the references cannot be separated the verdict is
  `arbitrated` when the sequence type can name the species and `unresolved` when it cannot, never
  `agree`, which is the verdict of a clean match.
* **A sample that matches nothing no longer carries 431 close matches.** A read set that shares no hashes
  sits at distance 1 from every reference, so every one of them tied with the closest and was reported as
  nearly as close — with 16 genomes that was invisible, with 432 it put a whole reference table in the
  sample's block of the quality assurance report and in `close_matches` in the JSON. A sample that was not
  identified now has no close matches, which is what the report was always meant to say. The caption over
  that table no longer claims a ranking the distances do not support either: where every reference is the
  same distance away it says so, instead of "closest first" over an arbitrary five of a 432-way tie and
  "the other 427 are further away" about references that are not.
* **One half of a custom pair no longer passes in silence.** `--sketch` and `--references` go
  together, and a 1.x command line that passed only `--sketch` fell back on the 16-genome table, which
  matched it. It now falls back on the extended table, whose IDs are GCA file names, so nothing in such a
  sketch could be named and every sample came back `ND` from a run that exited 0. A sketch none of whose
  genomes the reference table names is refused, by name, before any sample is read — but only when the
  sequence type is switched off as well, because otherwise the ST can still name the species and a sketch
  of your own genomes against the packaged ST table is a supported way to run. With the ST available the
  run warns instead and carries on, and the sample is reported `ST (arbitrated)` with no reference behind
  it, which now counts as a **warning** rather than `ok`: half the evidence the report exists to show is
  missing. A custom sketch that is a subset of a packaged set still works on its own, as before.
* **Callers of the library get the extended set too**, not just the command line: `reference_sketch()` and
  `load_references()` resolve to it when they are given no file, so a script that reaches
  `analyse_samples()` directly measures against the 432 genomes without having to know about
  `--reference-set`. There is deliberately no constant meaning "whichever set is the default" --
  `EXTENDED_*` and `VSNP_*` are the only names, so code that cares which genomes it used has to say which.
* **The no-biovar figures are now the default set's**, and they are starker: a biovar is recorded for only
  7 of the 432 genomes, because the NCBI records of the rest carry none, so 93% of the 240 *B. abortus*
  isolates of the PubMLST isolate database that have a numbered biovar are of a biovar no reference
  has, against 42% for
  vSNP's 16. The resolution argument is unchanged in kind -- the two *B. suis* biovar 2 genomes are 0.00163
  apart, further than the biovar 3 genome is from one of them -- and the reports say so in those terms.
* **The `Notes` column is now populated on nearly every sample**, where in 1.x it was usually empty: the
  432-genome set almost always has other genomes of the same species inside the close-match margin, so
  most samples carry a "nearly as close: N other Brucella x genome(s)" line. On the four tutorial read
  sets it goes from 1 of 4 samples to 4 of 4, and the longest cell from 218 to 691 characters. No
  column was added or removed anywhere, so a parser reading by position is unaffected — but a rule of
  "a non-empty Notes means look at this one" now fires on everything and needs rewriting.
* **The per-run output is bigger**, almost all of it `brucellaid_mash.tsv`, which holds every distance
  and so grows with the reference set: for four samples it goes from 5 kB to 170 kB, and
  `brucellaid.json` from 24 kB to 132 kB because its reference table has 432 entries instead of 16.
  The HTML report barely moves, because it still shows the closest few.
* A mismatch between a sketch and a reference table names the first eight genomes and counts the rest,
  instead of listing all of them: against the 432-genome table that warning was a single 19 kB line.
* Nothing about the MLST changes, the 150-shared-hash threshold is unchanged, and there is still no biovar.
* `validation/run_validation.sh` pins `--reference-set vsnp`. It reads those 16 genomes genome for genome,
  so its expected answers are their sketch IDs, and its negative control is *B. anthropi* ATCC 49188, which
  comes back `ND` only because the vSNP set does not carry it. The default set's evidence is the two
  at-scale validations instead.

### Also fixed in 2.0
None of these follows from the change of default; they were found reviewing it.

* **The legacy `--xlsx` report no longer writes a formula.** xlsxwriter sends any value starting with
  `=` to `write_formula()`, so a sample named `=HYPERLINK(...)` became a live formula and the cell
  displayed its result instead of the sample it belonged to. Every value is written as text now.
* **One bad file in an input folder can no longer take the run down.** A named pipe blocked the run for
  ever on `open()` and a link to a device read without end, because `-i` accepted anything the filename
  pattern matched while `-r1` had always required a regular file; such a file is now skipped by name.
  And an unexpected error in one sample is recorded against that sample, with its type named, instead of
  escaping the thread pool and ending the run before any report was written.
* **Smaller report fixes.** The `#` command line of each table goes through the same cell cleaning as the
  data, so a path holding a line break no longer splits it in two; a sample name that is only whitespace
  falls back to the name of its folder, and surrounding spaces are dropped, so the TSV, JSON, MultiQC
  file and report no longer disagree about what a sample is called — the files a name is built from still
  decide which reads pair with which, so two file names that differ only in spacing are a collision and
  not a merge; and the HTML report replaces bidirectional and control characters, which `html.escape`
  passes through and which can make a name render as something the machine-readable files do not contain.

### Also new in 2.0
* **The extended reference set**, which is what the paragraphs above are about. It identifies against
  432 curated genomes covering the whole *Brucella* clade -- all 30 named species, the former
  *Ochrobactrum* included
  (*B. papionis* is labelled nowhere in INSDC: its genome was found under "Brucella sp." and relabelled on
  sequence-type and distance evidence, see curation/papionis/README.md), plus the
  lineages with no name yet, which are reported as *Brucella sp.*. Derived from the 2,022
  assemblies behind the mashID *Brucella* pack and curated genome by genome
  (`scripts/extended_references_curation.tsv` logs every exclusion and why); validated on simulated reads
  of all 2,022, every held-out sample of a named species right once five genomes mislabelled in NCBI --
  caught by the curation itself -- are accounted for, and on all 4,042 usable public Brucella read sets in
  ENA, where it matches the submitter on 98.7% against 96.2% and wins every disagreement that has a right
  answer. Of its 52 non-matches, 20 are mislabelled records and 22 are read sets holding little or no
  Brucella; none is a wrong species named off a reference that was there.
  `scripts/make_extended_references.py` rebuilds the extended sketch from NCBI.
* **A loose match is now said to be loose.** When the closest reference is further than
  `--loose-distance` (0.005 by default), the call carries a note giving the distance and the shared hashes.
  It does not change the species: a sample can be the right species and still match loosely, because the
  library holds more than *Brucella*, because its lineage has no near relative in the sketch, or because
  the derived k-mer copy threshold is low for its depth — which is what two of the four tutorial
  libraries turned out to be under the 2.0 threshold, see
  [The k-mer copy threshold](K-mer-threshold). The
  threshold comes from every public *Brucella* read set in ENA, where nineteen in twenty calls that agreed
  with the submitter sat within 0.005, and all 185 agreeing calls beyond it shared fewer than 900 of the
  1000 hashes.
* A nearly-as-close reference of the **same species** as the closest one is counted, not listed: "nearly as
  close: 4 other Brucella melitensis genome(s)" instead of four file names that change nothing about the
  call. A reference of another species is still named in full, with its shared hashes.

## 1.2.1 (2026-11-02)
* **Concatenated gzip files are sized correctly.** The coverage estimate read the gzip trailer, which names
  only the *last* member of a file, so lane files joined with `cat` and the 64 KB blocks `bgzip` writes
  were estimated at a fraction of their size -- a deep bgzipped library could read as under 1x, keep the
  k-mer copy threshold at 2, and let the error k-mers back into the sketch. The head of the file is now
  decompressed member by member: several members mean the trailer is ignored outright, and a trailer that
  contradicts the measured compression ratio is ignored too.
* The reports are written in UTF-8 whatever the locale: a sample name with an accent no longer loses the
  reports of a finished run under `LC_ALL=C`.
* A gzip file with a corrupt stream inside makes the sample fail with the real error instead of crashing
  the run while its content was being sniffed for Seal.
* A coverage below 10x is written with one decimal, and a trace of reads reads "under 0.1x" -- never
  "roughly 0x coverage" beside a real base count.
* **Report tables are only as wide as their content**, in the HTML and the PDF alike, every column is
  left-aligned, and a rule opens and closes each table. A table stretched to the full page put unrelated
  amounts of air between its columns, and a right-aligned number sat away from its label. In the PDF, a
  column now gets the width its longest value asks for, and when a table is wider than the page only the
  columns with the most to say wrap.
* The `About ... coverage: sketching the k-mers seen at least N time(s)` log line now begins with the name
  of the sample it describes. With `--parallel`, two samples identify at once and their lines interleave
  with no way to tell whose coverage was whose.
* 345 tests. 1.2.0 claimed 100% coverage while four lines of the depth change had none; three tests close
  that, and the number is true again.

## 1.2.0 (2026-11-02)
**The k-mer copy threshold rises with the sequencing depth.** `mash dist -m 2` keeps every k-mer the read
set holds at least twice, and an error k-mer turns up about depth × error rate times, so past roughly 50x
the errors reach two copies on their own and crowd the true k-mers out of the bottom-1000 sketch. The
distance then grows with the depth instead of with the genetic distance: on 244 real libraries, the median
distance to the correct reference was 0.0219 at 100–150x where the assembly says 0.0001, and above about
240x the shared hashes fell under the 150 that `--min-hashes` asks for and the sample got **no call at
all**. The threshold is now derived from each sample — one more required copy per 20x of coverage, never
below 2, never above 20 — with the coverage estimated from the read files in a tenth of a second, without
decompressing them (the gzip trailer gives the size, the first thousand records give the bases per byte).

* Below about 50x the derived threshold is the same 2 as before and **nothing changes**, which is what the
  low end needs: a flat 10 would destroy a 5x sketch.
* At depth, everything does: on real libraries the medians moved from 666/1000 shared hashes at distance
  0.010637 to 966/1000 at 0.000843, and the 0.0005 close-match margin — which the sequence-type
  arbitration sits inside — means something again. The verdict changed on 14 of the 244 libraries, 6 of
  which had been reporting close matches off a margin of exactly zero.
* The reports say which threshold was used and the coverage it came from, and `brucellaid.json` carries
  `min_kmer_copies` and `estimated_coverage` under `closest_reference`.
* `--min-kmer-copies N` pins any fixed value, including the old 2.

**The validation covers depth now, and real reads.** Everything it simulated was 20x, where the derived
threshold returns the same 2 the old fixed value gave: the bug above would have sailed through it. Two
genomes are now simulated deep — *B. melitensis* 16M at 150x and *B. suis* 1330 at 250x — where the
derived threshold keeps 1000/1000 sketch hashes and the fixed 2 kept 860–902 and 728–748 in a six-sample
sweep, so `check_results.py` holds a floor of 950 and also checks the estimated coverage against the depth
the reads were simulated at and the threshold against the coverage. And the four public ENA libraries of
the [Tutorial](Tutorial) run next to the simulated samples: each must give the species its submitter gave
it, the reference and sequence type it is known to give, and the same species and reference as the
simulated sample of its own reference genome. All 22 samples pass, and the fix shows on real data:
SRR8550484 (~69x) moved from 422/1000 hashes at distance 0.0248 to 775/1000 at 0.0065, while ERR2136546
(47x) stayed exactly where it was. See [The k-mer copy threshold](K-mer-threshold).

**Upgrading:** no report or table changed shape. Results change only on read sets above about 50x, where
they were the ones being distorted: distances shrink toward what an assembly would give, shared hashes
rise, and a sample deep enough to have been `ND` can get its call back.

* 331 tests, 100% statement and branch coverage.

## 1.1.0 (2026-10-31)
**The biovar column is gone.** It never was a biovar call — it held the biovar of the closest
reference genome — and it was read as one. Three things rule out a real biovar from this data, and
none of them is a matter of trying harder:

* **Coverage.** The sketch holds *B. abortus* biovars 1 and 3 and no others, so **43%** of the
  biovar-labelled *B. abortus* isolates of the PubMLST isolate database are of a biovar it has no
  reference for. *B. suis* bv 5 is missing too.
* **Resolution.** The two *B. melitensis* **biovar 1** references of the sketch are **0.000756** apart,
  *further* than the biovar 2 and biovar 3 references are from each other (0.000706). Two genomes of one
  biovar are further apart than genomes of two biovars, so the ranking carries no biovar signal.
* **Biology.** A biovar is a phenotype and is not monophyletic: the best any rule from sequence type to
  biovar could manage is **66%** for *B. abortus* and **63%** for *B. melitensis*. ST 8 alone is biovar 3
  in 47 isolates, biovar 1 in 34 and biovar 2 in 22. Only *B. suis* is tidy, at 96%.

So `Biovar` is out of `brucellaid.tsv`, out of the MultiQC table, out of `brucellaid.json` and out of the
reports. What stays, because it is a fact about a reference rather than a claim about a sample: the biovar
of each reference genome in `brucellaid_mash.tsv`, and the biovars recorded for the sample's sequence type,
with their counts, in the quality assurance report. For a biovar, use classical biotyping or
[brucetyper](https://github.com/duceppemo/brucetyper). See [How it works](How-it-works#no-biovar).

**The report is easier to read.**

* A sample block lists the **five closest** reference genomes, plus any other near enough to the closest
  that the distance cannot separate them, instead of all sixteen. `brucellaid_mash.tsv` still holds every
  distance, and the report says how many it left out.
* **One sample to a printed page.** The fact lists are laid out in two columns and the tables across the
  page, which is what makes a sample fit.
* The summary is one table with a labelled row per group rather than one table per group, statuses are
  shown as coloured pills, sections are numbered, and the type, spacing and tables were reworked in both
  the light and the dark theme.
* The allele table merges the runner-up and its read count into one column.
* **The PDF was rebuilt to match**, rather than merely holding the same words: the same palette,
  numbered sections, one grouped summary table with coloured statuses, two columns of facts above each
  table, quieter tables with a rule under the header instead of a grid, a tinted box to sign, a rule
  above the footer, and one sample to a page there too.

**A quality assurance report.** Every run now writes `brucellaid_report.html`: one self-contained file —
no script, no stylesheet, no font, nothing fetched from anywhere — holding what a result rests on rather
than only the result.

* **`brucellaid_report.html`**, written on every run (`--no-html` skips it): a summary with the samples
  grouped by what they were called and a box for a reviewer's name, date and signature; the notes,
  warnings and errors one row each; a block per sample with the five closest references, and any other
  the distance cannot separate from the first (the closest highlighted, the near ones shaded;
  `brucellaid_mash.tsv` holds every distance), how the species was decided and what the sequence type
  contributed, the nine allele calls with the margin each had over the allele it beat, and the input files
  with their size, date and MD5; what the numbers mean; and the run itself — operator, computer, command,
  versions, and **the MD5 of every reference data file**. It prints to clean pages as it stands.
* **`brucellaid_mqc.json`**, also on every run: a MultiQC custom content table, so `multiqc results/` shows
  a BrucellaID section with the species, what decided it, the sequence type and the status of each sample.
* **`--pdf`** writes the same report as a PDF, with "Page 3 of 9" in every footer. It needs reportlab,
  which is an optional extra: `pip install 'brucellaid[pdf]'`. BrucellaID itself still installs with
  nothing but the standard library, and the HTML report holds the same information.
* **`--no-md5`** skips the checksums of the input files, which cost a second or two per gigabyte, and
  **`--operator NAME`** puts a name other than the user account on the report.
* The report carries the BrucellaID logo, inlined as svg in the HTML (one drawing, themed with CSS
  variables, so it follows the light and the dark page) and drawn from a png in the PDF. Both ship in
  the package, written by `assets/make_logo.py`.
* An example report, [HTML](https://duceppemo.github.io/BrucellaID/example_report.html) and
  [PDF](https://duceppemo.github.io/BrucellaID/example_report.pdf): the four read sets of the
  [Tutorial](Tutorial), rendered with the computer they ran on left out by
  `scripts/make_example_report.py`, and published with GitHub Pages from `docs/` so that a browser
  renders them instead of showing their source.
* Each sample now has a status — `ok`, `warning` (not identified, no ST, a `new` profile, or the two kinds
  of evidence disagreeing) or `failed` — in the reports and in the MultiQC table.
* `brucellaid.json` gains, per sample, the `status`, the input files with their checksums and, in
  `species_call`, `st_not_used`: why the sequence type was left out of the species call when it was. Its
  run block gains the operator, the user, the computer, the operating system, the working directory and
  the checksummed `reference_data`.

**The sequence type now decides the species MASH cannot.** *B. canis* ATCC 23365 and *B. suis* bv. 4 strain
40 are 0.00024 apart in the sketch — the closest pair of the sixteen, closer than two *B. suis* references
are to each other — so which one a *canis* or a *suis* bv. 4 read set lands on first is noise. The two
species share no sequence type at all, so the ST settles it.

* New data file `brucellaid/data/mlst_st_species.tsv`: the species of each ST of the scheme, counted over
  the 3,426 isolates of the PubMLST *Brucella* isolate database, with the share of them that agree and the
  biovars recorded. 118 STs, 105 of them seen in one species only. Built by the new
  `scripts/make_st_species.py`.
* Where a reference of another species is within 0.0005 of the closest one, the species is no longer the
  closest reference's by default: a sequence type whose isolates agree chooses between the species MASH left
  within reach, and the biovar and the vSNP code follow the species chosen. An ST that contradicts a species
  MASH picked cleanly is reported but never obeyed. See [How it works](How-it-works).
* Two new columns in `brucellaid.tsv`: **`SpeciesSource`** (`MASH`, `MASH+ST (agree)`,
  `MASH+ST (arbitrated)`, `ST (arbitrated)`, `MASH (unresolved)`, `MASH (conflict)`) and **`STSpecies`**.
  `brucellaid.json` gains a `species_call` object. Anything reading the main table by column *number* has to
  be updated; by column name, nothing changes.
* New options: `--no-st-species` (decide the species as 1.0.0 did), `--min-st-purity` (0.95),
  `--min-st-isolates` (3) and `--st-species FILE`.
* An ST is only used when at least `--min-st-purity` of at least `--min-st-isolates` isolates agree on a
  species, and never when that species is "*Brucella* sp.". The labels are the submitters' and some are
  wrong, so a split ST is left out of the decision rather than settled by a majority.
* `CONTRIBUTING.md` listed the allele and sequence type counts of vSNP's old copy of the database, 263 and
  81; they are PubMLST's 309 and 122.
* The wiki front page carries the logo, in the light and the dark variant, as the README does.
* The Zenodo DOI of the archived releases is recorded in `CITATION.cff`, in the README badge and in the
  citation on this wiki: [10.5281/zenodo.23845956](https://doi.org/10.5281/zenodo.23845956) resolves to
  the latest version.

## 1.0.0 (2026-10-31)
BrucellaID becomes an installable, tested package. It answers the same two questions as the first version —
which *Brucella* species a read set is, and which sequence type — and it answers the first one exactly as
before. What changed is everything around that, and how the MLST is called.

**Packaging and interface**
* Installable with `pip install brucellaid` (and from bioconda once the recipe is merged): a
  `brucellaid` command, and `python -m brucellaid`.
* **The two OLC dependencies are gone**: `accessoryFunctions` (OLCTools), which every function imported, and
  `MLSTsippr` (sipprverse), which ran the MLST. Neither installs from PyPI, so neither could be packaged.
  BrucellaID now needs nothing but the Python standard library, plus MASH and BBTools.
* **The MLST database is in the package**, so `-d` is gone. The input arguments are `-r1`/`-r2` for one
  sample or `-i` for a folder, instead of `-s`.
* The main table is **printed to standard output**, with the messages on standard error, so
  `brucellaid ... > results.tsv` works.
* A folder is searched recursively, R1/R2 files are paired, and samples run `-p` at a time.
* A sample that fails no longer stops the run: its error is a column of the report, and the exit code says
  that something failed.
* fasta input is **refused with a reason** instead of being typed badly: both steps count reads.

**Reports**
* `brucellaid.tsv` (one row per sample), `brucellaid_mlst.tsv` (one row per locus, with the reads behind
  each call and the runner-up), `brucellaid_mash.tsv` (every distance to every reference) and
  `brucellaid.json` (everything, including the versions of MASH, BBTools and Java).
* `--xlsx` still writes `stat_alignment_summary_<date>.xlsx`, with the same four columns and the vSNP
  species code in the `species` column, for the pipelines that read it.

**Identification**
* Unchanged in substance: the same 16 *Brucella* reference genomes from vSNP, the same MASH sketch, the same
  threshold of 150 shared hashes out of 1000. The reads of a sample are still sketched as one read set.
* What is new: the species, biovar, strain and NCBI accessions of each reference are in a table that ships
  with the package, the distance to **every** reference is reported, and a reference that is nearly as close
  as the closest one is called out in a note. The sketch and the table are checked against each other at
  startup.

**MLST**
* The alleles are called with **Seal** from BBTools instead of MLSTsippr's bowtie2 and samtools. Seal counts
  a read for every allele of a locus it shares the highest number of 25-mers with, and for no other
  (`ambiguous=all`, `clearzone=0`, no mismatch), so the allele the sample carries is counted for every read
  of the locus and the others only where they agree with it. The allele with the highest count wins, and a
  tie is reported rather than broken. See [How it works](How-it-works).
* Every call now comes with the reads behind it and the runner-up, so the margin is visible.
* **The database is PubMLST's**, downloaded by `scripts/make_mlst_data.py` on 2026-10-31 (the MLST (9 locus)
  scheme, data submitted up to 2024-12-31): 309 alleles and 122 sequence types, where the copy shipped with
  the old dependencies folder had 263 alleles and 81. That copy had also drifted: it still held a `gyrB_23`
  that PubMLST has deleted, and that sequence matched the *B. abortus* 9-941 genome better than the real
  gyrB_1, which turned ST 1 into a new profile. The 24 other gyrB alleles of the old copy are identical to
  PubMLST's.

**Around the code**
* 194 tests, 100% statement and branch coverage, run on Python 3.10, 3.12 and 3.14 by GitHub Actions.
* A [Validation suite](Validation-suite) page, and `validation/run_validation.sh`, which reproduces it from public
  genomes.
* A [wiki](Home), a [tutorial](Tutorial) on public read sets, a bioconda recipe, `CITATION.cff`, and
  `examples/run_example.sh`, which runs in a few seconds on data it makes itself.
* The 193 MB of per-reference vSNP `script_dependents` folders are no longer in the repository: no code read
  them. They remain in the history and in the
  [v0.1 release](https://github.com/duceppemo/BrucellaID/releases/tag/v0.1), which is the last state of the
  clone-and-run version.
