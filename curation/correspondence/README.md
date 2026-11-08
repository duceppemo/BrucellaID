# Correspondence

One message per recipient, plain text, ready to send. None carries a `To:` line: neither ENA nor NCBI
publishes submitter email addresses, so the recipient has to be chosen per record — the submitting
organisation, the corresponding author of an associated paper, or the archive itself.

There are two kinds of file here. The `archive_*` letters are addressed to the **curation teams** and
consolidate every record in that archive into one request — use these if you would rather go through the
archive than chase individual submitters. The rest are addressed to the **submitters**, one per
submitting group, and say "your records"; they are not interchangeable with the archive letters.

| To the archives | Records | Status |
| --- | --- | --- |
| `archive_ena_2026-11.txt` | 13 BioSamples to change, 11 to annotate, 1 flagged as a probable mixture | **sent 2026-11-07** via the [ENA support form](https://www.ebi.ac.uk/ena/browser/support) |
| `archive_ncbi_2026-11.txt` | 5 BioSamples to change, 10 to annotate, 3 explicitly needing no action | **sent 2026-11-07** to `biosamplehelp@ncbi.nlm.nih.gov` |
| `archive_ena_followup_2026-11.txt` | `SAMEA114240991`, resolved by type-strain ANI after the first message | **sent 2026-11-07** via the ENA support form |
| `archive_ncbi_followup_2026-11.txt` | `SAMN07581394`, `SAMN13812480`, same | **sent 2026-11-07** to `biosamplehelp@ncbi.nlm.nih.gov` |

| To the submitters | Records | Submitter | Status |
| --- | --- | --- | --- |
| `genbank_2026-11.txt` | 5 assemblies | GenBank | **sent 2026-11-07** |
| `prjeb30030_2026-11.txt` | 11 BioSamples | Instituto Nacional de Saúde Dr. Ricardo Jorge | **sent 2026-11-07** |
| `prjeb4782_2026-11.txt` | 3 BioSamples + 10 runs | Wellcome Sanger Institute | **sent 2026-11-07** |
| `broad_2026-11.txt` | 2 BioSamples + 8 runs | Broad Institute | **sent 2026-11-07** |
| `single_prjeb22377.txt` | `ERR3218177` | Wellcome Sanger Institute | **sent 2026-11-07** |
| `single_prjna649281.txt` | `SRR31560207` | University of Pretoria | **sent 2026-11-07** |
| `single_prjna610713.txt` | `SRR11526759` | Sardar Vallabhbhai Patel Univ. of Agriculture and Technology | **sent 2026-11-07** |
| `single_prjna1183633.txt` | `SRR31357708` | Institute of Genetics and Developmental Biology | **sent 2026-11-07** |
| `single_prjna1290416.txt` | `SRR35797377` | Shandong Center for Disease Control and Prevention | **sent 2026-11-07** — **not** a correction request |
| `single_prjna945449.txt` | `SRR25828922` | NIAID | **sent 2026-11-07** |

`SRR24420454` (PRJNA966239) gets no message: its library holds a third of a *Brucella* genome at 0.5× and
nothing else, so there is nothing to report beyond insufficient depth.

## How a correction actually gets made

Neither archive lets a third party edit someone else's record, so each message offers the route rather
than assuming the recipient knows it.

* **NCBI** (SRR/SAMN records): organism changes on a BioSample go through `biosamplehelp@ncbi.nlm.nih.gov`,
  with `pd-help@ncbi.nlm.nih.gov` for pathogen data. NCBI notes it often cannot reach the original
  submitter, so the help desk may act directly.
* **ENA** (ERR/SAMEA records): no third-party error channel. The routes are the
  [support form](https://www.ebi.ac.uk/ena/browser/support) and the ELIXIR Clearinghouse, which exists for
  third-party curation of records one did not submit.

## Two messages that are not what they look like

* `single_prjna610713.txt` — the submitter's own project title already says *B. melitensis*; only the
  run's organism field says *abortus*, and our reading agrees with their title. It also flags that
  "strain 2308" is conventionally *B. abortus* 2308, which only their records can settle.
* `single_prjna1290416.txt` — their species is **correct**. The library is mixed, which is why
  sketch-based tools report the run as not identified, and the message explains that rather than asking
  for anything.

## Every letter here has been sent, and they are kept verbatim

All fourteen went out on 2026-11-07: the five archive and GenBank letters first, then the nine
per-submitter ones. Three of the first five describe a "431-genome reference set", because the set was 431
when they went and is 432 now; the nine that went later were corrected to 432 before they were sent. The
earlier three are left exactly as sent, because that is what their recipients hold.

Do not "fix" a figure in any of these files. They are a record of what was sent, not a description of the
current reference set, and `../mislabelled_records.md` records the date each one went.
