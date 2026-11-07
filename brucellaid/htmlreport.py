"""
The quality assurance report: one self-contained HTML file holding everything needed to check a run.

Self-contained is the whole point. No stylesheet, no script, no font and no image is fetched from anywhere,
so the file can be attached to an email, dropped into a laboratory information system or archived for seven
years and still render exactly as it did on the day it was written. It prints to clean pages as it stands.

What it holds, in the order a reviewer reads it: a summary of the run with the samples grouped by what they
were called and a box to sign; every warning and error, one row each; then, per sample, the identification
with the five closest reference genomes (and any other the distance cannot separate from the first), how
the species was decided, the nine
allele calls with the margin each one had, and the input files with their checksums; then what the numbers
mean; and last, the run itself -- operator, machine, command, versions and the checksum of every reference
data file it measured against.
"""

import pathlib
import re
from functools import cache
from html import escape as _escape

from . import __version__
from .errors import BrucellaIDError
from .mash import kmer_copies_note
from .mlst import UNKNOWN
from .qa import (
    FAILED,
    OK,
    WARNING,
    close_hits,
    counted,
    groups,
    hits_caption,
    margin,
    status,
    top_hits,
    warnings,
)
from .references import data_path
from .species import AGREE, ARBITRATED, CONFLICT, UNRESOLVED

REPORT_FILE = 'brucellaid_report.html'
LOGO_FILE = 'logo.svg'
TITLE = 'BrucellaID quality assurance report'
HEADING = 'Quality assurance report'

# One stylesheet, inline, with the colours as tokens so the dark variant only has to redefine them. The
# print rules matter as much as the screen ones: this is a document that gets signed on paper.
CSS = """
:root {
  color-scheme: light dark;
  --bg: #ffffff; --panel: #f6f7f9; --raised: #ffffff; --ink: #14171a; --soft-ink: #3d454d;
  --muted: #6b747d; --line: #e3e7ec; --hair: #eef1f4;
  --accent: #15639f; --accent-soft: #eaf2fa; --accent-line: #bcd7ee;
  --ok: #16723c; --ok-soft: #e9f5ed; --warn: #8a5000; --warn-soft: #fdf1df;
  --bad: #a51f1f; --bad-soft: #fcebeb;
  --shadow: 0 1px 2px rgba(16, 24, 32, 0.05), 0 1px 10px rgba(16, 24, 32, 0.04);
  --mono: ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace;
  --logo-ink: #2f2a3d; --logo-accent: #7d4a9e; --logo-reference: #9aa5ad; --logo-muted: #5d6a72;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #111418; --panel: #191d22; --raised: #171b20; --ink: #e9ecef; --soft-ink: #c7ced5;
    --muted: #949da6; --line: #2a3037; --hair: #22272d;
    --accent: #6cb2ef; --accent-soft: #15242f; --accent-line: #2b4559;
    --ok: #5ec488; --ok-soft: #12241a; --warn: #e0a95c; --warn-soft: #2a2114;
    --bad: #ef8d8d; --bad-soft: #2b1818;
    --shadow: 0 1px 2px rgba(0, 0, 0, 0.4);
    --logo-ink: #f0edf5; --logo-accent: #b185d8; --logo-reference: #6d787f; --logo-muted: #a9b6bd;
  }
}
:root[data-theme="dark"] {
  --bg: #111418; --panel: #191d22; --raised: #171b20; --ink: #e9ecef; --soft-ink: #c7ced5;
  --muted: #949da6; --line: #2a3037; --hair: #22272d;
  --accent: #6cb2ef; --accent-soft: #15242f; --accent-line: #2b4559;
  --ok: #5ec488; --ok-soft: #12241a; --warn: #e0a95c; --warn-soft: #2a2114;
  --bad: #ef8d8d; --bad-soft: #2b1818;
  --shadow: 0 1px 2px rgba(0, 0, 0, 0.4);
  --logo-ink: #f0edf5; --logo-accent: #b185d8; --logo-reference: #6d787f; --logo-muted: #a9b6bd;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; padding: 0 20px 72px; background: var(--bg); color: var(--ink);
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  font-variant-numeric: tabular-nums;
}
main { max-width: 1120px; margin: 0 auto; }
.masthead { padding: 34px 0 6px; }
.logo svg { width: 326px; max-width: 100%; height: auto; display: block; }
h1 { font-size: 1.55rem; font-weight: 650; margin: 14px 0 2px; letter-spacing: -0.015em; }
.sub { color: var(--muted); font-size: 0.92rem; margin: 0; }
h2 {
  font-size: 1.22rem; font-weight: 650; letter-spacing: -0.01em; margin: 46px 0 14px;
  padding-bottom: 8px; border-bottom: 1px solid var(--line); display: flex; align-items: baseline; gap: 10px;
}
h2 .n {
  font-size: 0.72rem; font-weight: 700; color: var(--accent); background: var(--accent-soft);
  border: 1px solid var(--accent-line); border-radius: 5px; padding: 1px 7px; letter-spacing: 0.02em;
}
h3 { font-size: 1.02rem; font-weight: 650; margin: 24px 0 8px; }
h4 {
  font-size: 0.72rem; font-weight: 700; margin: 18px 0 7px; color: var(--muted);
  text-transform: uppercase; letter-spacing: 0.08em;
}
p, ul, ol { margin: 9px 0; }
a { color: var(--accent); text-underline-offset: 2px; }
code, .mono { font-family: var(--mono); font-size: 0.86em; }

/* The run at a glance */
.chips { display: flex; flex-wrap: wrap; gap: 8px; margin: 0 0 16px; }
.chip {
  display: inline-flex; align-items: baseline; gap: 6px; padding: 5px 12px; border-radius: 999px;
  background: var(--panel); border: 1px solid var(--line); font-size: 0.84rem; color: var(--soft-ink);
}
.chip b { font-size: 1rem; font-weight: 650; color: var(--ink); }
.chip.ok { background: var(--ok-soft); border-color: transparent; color: var(--ok); }
.chip.ok b { color: var(--ok); }
.chip.warning { background: var(--warn-soft); border-color: transparent; color: var(--warn); }
.chip.warning b { color: var(--warn); }
.chip.failed { background: var(--bad-soft); border-color: transparent; color: var(--bad); }
.chip.failed b { color: var(--bad); }

.pill {
  display: inline-flex; align-items: center; gap: 5px; padding: 2px 9px; border-radius: 999px;
  font-size: 0.76rem; font-weight: 600; border: 1px solid var(--line); background: var(--panel);
  color: var(--soft-ink); white-space: nowrap;
}
.pill::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: currentColor; }
.pill.ok { background: var(--ok-soft); border-color: transparent; color: var(--ok); }
.pill.warning { background: var(--warn-soft); border-color: transparent; color: var(--warn); }
.pill.failed { background: var(--bad-soft); border-color: transparent; color: var(--bad); }
.pill.decided { background: var(--accent-soft); border-color: transparent; color: var(--accent); }
.pill.plain::before { display: none; }

table { border-collapse: collapse; margin: 8px 0 2px; font-size: 0.86rem; border-top: 1px solid var(--line); }
caption { text-align: left; color: var(--muted); font-size: 0.82rem; padding: 0 0 7px; line-height: 1.45; }
th {
  text-align: left; padding: 7px 10px; background: transparent; color: var(--muted); font-weight: 650;
  font-size: 0.74rem; text-transform: uppercase; letter-spacing: 0.05em; white-space: nowrap;
  border-bottom: 1px solid var(--line);
}
td { padding: 7px 10px; border-bottom: 1px solid var(--hair); vertical-align: top; color: var(--soft-ink); }
/* The closing rule: the last row's own border, because at equal width a cell border wins the collapse */
tbody tr:last-child td { border-bottom: 1px solid var(--line); }
td:first-child { color: var(--ink); }
td.nw { white-space: nowrap; }
tr.best td { background: var(--accent-soft); color: var(--ink); }
tr.best td:first-child { box-shadow: inset 3px 0 0 var(--accent); font-weight: 600; }
tr.close td { background: var(--panel); }
tr.warning td:first-child { box-shadow: inset 3px 0 0 var(--warn); }
tr.failed td:first-child { box-shadow: inset 3px 0 0 var(--bad); }

tr.grouprow td {
  background: var(--panel); font-weight: 650; color: var(--ink); font-size: 0.85rem;
  border-bottom: 1px solid var(--line); padding: 7px 10px;
}
tr.grouprow .count { font-weight: 400; color: var(--muted); font-size: 0.8rem; margin-left: 10px; }

.panel {
  background: var(--panel); border: 1px solid var(--line); border-radius: 10px; padding: 14px 16px;
  margin: 14px 0;
}
.note { background: var(--accent-soft); border-color: var(--accent-line); color: var(--soft-ink); }

.sample {
  border: 1px solid var(--line); border-radius: 12px; padding: 2px 18px 18px; margin: 16px 0;
  background: var(--raised); box-shadow: var(--shadow); border-left: 3px solid var(--line);
}
.sample.ok { border-left-color: var(--ok); }
.sample.warning { border-left-color: var(--warn); }
.sample.failed { border-left-color: var(--bad); }
.sample h3 { display: flex; flex-wrap: wrap; gap: 9px; align-items: center; margin-top: 16px; }
.sample h3 .name { font-size: 1.08rem; }

/* Two columns of facts: the page is wide, and a sample should stay on one printed page */
.cols { display: grid; grid-template-columns: 1fr 1fr; gap: 4px 34px; align-items: start; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: 2px 14px; margin: 6px 0 2px; }
dl.facts dt { color: var(--muted); font-size: 0.82rem; line-height: 1.5; }
dl.facts dd { margin: 0; font-size: 0.86rem; line-height: 1.5; }

.sign { display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin-top: 10px; }
.sign .line { height: 44px; border-bottom: 1px solid var(--soft-ink); }
.sign .what { padding-top: 5px; font-size: 0.78rem; color: var(--muted); }

nav.toc { font-size: 0.9rem; margin: 18px 0 0; }
nav.toc ol { margin: 6px 0 0; padding-left: 20px; color: var(--muted); }
nav.toc li { margin: 2px 0; }
footer {
  margin-top: 48px; padding-top: 12px; border-top: 1px solid var(--line); color: var(--muted);
  font-size: 0.78rem;
}
.wrap { overflow-wrap: anywhere; }

/* Screen only: an A4 page is about 750 CSS pixels wide, so a plain max-width rule would collapse the
   columns when printing, which is exactly where they are needed. */
@media screen and (max-width: 760px) {
  .cols, .sign { grid-template-columns: 1fr; }
  dl.facts { grid-template-columns: 1fr; }
  dl.facts dd { margin-bottom: 7px; }
  table { font-size: 0.78rem; }
  th, td { padding: 6px 7px; }
}

@media print {
  @page { size: A4; margin: 13mm 11mm 15mm; }
  body { padding: 0; font-size: 9.2pt; line-height: 1.42; background: #fff; }
  :root {
    --bg: #fff; --panel: #f4f5f7; --raised: #fff; --ink: #000; --soft-ink: #1f2429; --muted: #4a5157;
    --line: #c2c8ce; --hair: #e2e6ea; --accent: #0d4f82; --accent-soft: #eaf2fa; --shadow: none;
  }
  nav.toc, .noprint { display: none; }
  .logo svg { width: 196px; }
  h1 { font-size: 15pt; }
  h2 { font-size: 12pt; margin: 0 0 8pt; break-before: page; }
  h2.first { break-before: avoid; }
  h2, h3, h4 { break-after: avoid; }
  /* One sample to a page: with the five closest references it fits, and a reviewer reads one at a time.
     Everything below is what buys the room for that -- a sample that spills makes two half-empty pages. */
  .sample {
    break-before: page; break-inside: avoid; box-shadow: none; margin: 0;
    border: 0; border-left: 2pt solid; border-radius: 0; padding: 0 0 0 8pt;
  }
  .sample:first-of-type { break-before: avoid; }
  .sample h3 { margin: 0 0 2pt; gap: 7px; }
  table { font-size: 8.3pt; margin: 3pt 0 1pt; }
  th { padding: 3pt 5pt 2.5pt; font-size: 7pt; }
  td { padding: 2.6pt 5pt; }
  caption { font-size: 7.4pt; padding-bottom: 3pt; line-height: 1.35; }
  h4 { margin: 8pt 0 2pt; font-size: 7.2pt; }
  h3 { margin: 10pt 0 3pt; }
  p { margin: 4pt 0; }
  dl.facts { gap: 0 12px; margin: 3pt 0 0; }
  dl.facts dt, dl.facts dd { font-size: 8.3pt; line-height: 1.34; }
  .cols { gap: 0 26px; }
  .chips { margin-bottom: 8pt; }
  table, .panel, .sign { break-inside: avoid; }
  thead { display: table-header-group; }
  a { color: #000; text-decoration: none; }
}
"""


@cache
def logo():
    """
    The logo of the package, as svg markup to inline.

    Inlined rather than linked, because a report that fetches an image is no longer one file, and the
    colours are CSS variables, so the same markup follows the light and the dark theme. A package without
    the file still writes a report; it simply has no logo on it.
    """
    try:
        markup = pathlib.Path(data_path(LOGO_FILE)).read_text(encoding='utf-8')
    except (BrucellaIDError, OSError):
        return ''
    # The file carries its own width and height; the stylesheet sizes it instead
    return re.sub(r'\s+(width|height)="[^"]*"', '', markup, count=2)


def escape(value):
    """A value as HTML text. Everything that reaches the report goes through this."""
    return _escape('' if value is None else str(value), quote=True)


def anchor(name):
    """A sample name as an HTML id: sample names come from file names and may hold anything."""
    return 'sample-' + re.sub(r'[^A-Za-z0-9_-]+', '-', str(name)).strip('-').lower()


def pill(text, kind=''):
    """A small coloured label: a status, or what decided a species."""
    return '<span class="pill {}">{}</span>'.format(escape(kind), escape(text)) if kind else \
        '<span class="pill plain">{}</span>'.format(escape(text))


def heading(number, title, anchor_id):
    """A numbered section heading, so that the contents list and the page agree."""
    return '<h2 id="{}"{}><span class="n">{}</span> {}</h2>'.format(
        escape(anchor_id), ' class="first"' if number == 1 else '', number, escape(title))


class Raw(str):
    """A cell value that is markup already, such as a pill, and must not be escaped again."""


def cell(value, nowrap=False):
    text = value if isinstance(value, Raw) else escape(value)
    return '<td class="nw">{}</td>'.format(text) if nowrap else '<td>{}</td>'.format(text)


def grouped_table(headers, sections, caption):
    """
    One table whose rows are divided by a labelled row, instead of one table per group.

    Four groups of one sample each should not be four tables with four header rows; the groups are
    separated by a row that names them.

    :param sections: [(label, count label, rows)], a row being values or (values, row class)
    """
    out = ['<table>', '<caption>{}</caption>'.format(escape(caption)),
           '<thead><tr>{}</tr></thead>'.format(''.join('<th>{}</th>'.format(escape(h)) for h in headers)),
           '<tbody>']
    for label, count, rows in sections:
        out.append('<tr class="grouprow"><td colspan="{}">{}<span class="count">{}</span></td></tr>'.format(
            len(headers), escape(label), escape(count)))
        for row in rows:
            values, row_class = row if isinstance(row, tuple) else (row, '')
            cells = ''.join(cell(value) for value in values)
            out.append('<tr class="{}">{}</tr>'.format(escape(row_class), cells) if row_class
                       else '<tr>{}</tr>'.format(cells))
    out.append('</tbody></table>')
    return '\n'.join(out)


def table(headers, rows, caption='', classes=(), nowrap=()):
    """
    A table from already escaped-on-write values.

    :param headers: column titles
    :param rows: list of rows of values, or (values, row class) pairs when a row is highlighted
    :param nowrap: indexes of the columns that must stay on one line, such as a date
    """
    out = ['<table class="{}">'.format(' '.join(classes)) if classes else '<table>']
    if caption:
        out.append('<caption>{}</caption>'.format(escape(caption)))
    out.append('<thead><tr>{}</tr></thead>'.format(
        ''.join('<th>{}</th>'.format(escape(h)) for h in headers)))
    out.append('<tbody>')
    for row in rows:
        values, row_class = row if isinstance(row, tuple) else (row, '')
        cells = ''.join(cell(v, i in nowrap) for i, v in enumerate(values))
        out.append('<tr class="{}">{}</tr>'.format(escape(row_class), cells) if row_class
                   else '<tr>{}</tr>'.format(cells))
    out.append('</tbody></table>')
    return '\n'.join(out)


def facts(pairs):
    """A definition list of (term, value) pairs, skipping the ones with nothing to say."""
    out = ['<dl class="facts">']
    for term, value in pairs:
        if value in (None, '', ()):
            continue
        out.append('<dt>{}</dt><dd class="wrap">{}</dd>'.format(escape(term), escape(value)))
    out.append('</dl>')
    return '\n'.join(out)


VERDICTS = {
    AGREE: 'The sequence type names the same species as the closest reference genome, which named it on its '
           'own: the two kinds of evidence agree.',
    ARBITRATED: 'A reference of a second species was near enough to the closest one that the distance did '
                'not really choose. The sequence type chose between the species it left within reach.',
    CONFLICT: 'The sequence type belongs to a species that was not among the ones the distance left within '
              'reach. The species reported is the one MASH measured, on a thousand min-hashes against nine '
              'loci, and this sample needs a second look.',
    UNRESOLVED: 'A reference of a second species was just as close, and the sequence type could not choose '
                'between them. The species reported is the closest reference genome\'s.',
}


def summary_section(results, run):
    """The first thing a reviewer reads: what the run called, grouped, and a box to sign."""
    counts = counted(results)
    chips = ['<span class="chip"><b>{}</b> sample(s)</span>'.format(len(results))]
    for kind, label in ((OK, 'ok'), (WARNING, 'warning'), (FAILED, 'failed')):
        if counts[kind]:
            chips.append('<span class="chip {}"><b>{}</b> {}</span>'.format(kind, counts[kind], label))
    out = [heading(1, 'Summary', 'summary'),
           '<div class="chips">{}</div>'.format(''.join(chips)),
           facts([('Run started', run.get('started', '')),
                  ('Run finished', run.get('finished', '')),
                  ('Operator', run.get('operator', '')),
                  ('Command', run.get('command', ''))])]
    sections = []
    for (species, st), members in groups(results):
        label = '{} \u00b7 ST {}'.format(species, st) if species else 'Failed'
        rows = []
        for result in members:
            state = status(result)
            rows.append(([result.sample, result.species or UNKNOWN,
                          result.call.source if result.call else '',
                          result.call.verdict if result.call and result.call.verdict else '',
                          result.sequence_type, result.mlst.profile if result.mlst else '',
                          result.closest_reference, Raw(pill(state, state))],
                         state if state != OK else ''))
        sections.append((label, '{} sample(s)'.format(len(members)), rows))
    out.append(grouped_table(['Sample', 'Species', 'Species from', 'Decision', 'ST', 'Allele profile',
                              'Closest reference', 'Status'], sections,
                             'Samples called the same way are grouped, the largest group first and the '
                             'failed samples last.'))
    out.append('<div class="panel"><h4>Reviewed by</h4><div class="sign">'
               '<div class="line"></div><div class="line"></div><div class="line"></div>'
               '<div class="what">Name</div><div class="what">Date</div><div class="what">Signature</div>'
               '</div></div>')
    return out


def warnings_section(results):
    """Every warning and error of the run, one row each."""
    rows = [([sample, Raw(pill(state, state)), text], state if state != OK else '')
            for sample, state, text in warnings(results)]
    out = [heading(2, 'Notes, warnings and errors', 'warnings')]
    if not rows:
        out.append('<p>None: every sample was identified and typed without a note.</p>')
        return out
    out.append(table(['Sample', 'Status', 'What it says'], rows,
                     caption='Errors first, then the notes of each sample, in the order they were raised. '
                             'A note on a sample whose status is "ok" is a remark, not a problem.'))
    return out


def identification_facts(result):
    """What the closest reference was, and how close."""
    out = ['<h4>Identification: the MASH distance</h4>']
    if result.identification is None:
        out.append('<p>MASH did not run for this sample.</p>')
        return out
    best = result.identification.best
    out.append(facts([
        ('Closest reference', result.closest_reference),
        ('What it is', result.reference.description if result.reference else 'not in the reference table'),
        ('NCBI accessions', result.reference.accessions.replace(';', ', ') if result.reference
         and result.reference.accessions else ''),
        ('MASH distance', '{:.6f}'.format(best.distance)),
        ('Estimated identity', '{:.3f}%'.format(best.identity)),
        ('Shared min-hashes', '{}/{} (at least {} required)'.format(
            best.shared_hashes, best.total_hashes, result.identification.min_hashes)),
        ('K-mer copies required', kmer_copies_note(result.identification)),
    ]))
    return out


def identification_table(result, references):
    """The closest reference genomes, and how many were left out of the list."""
    if result.identification is None:
        return []
    best = result.identification.best
    out = []
    all_hits = close_hits(result, references)
    shown, hidden = top_hits(all_hits)
    rows = []
    for hit, reference, close in shown:
        row_class = 'best' if hit is best else ('close' if close else '')
        rows.append(([hit.reference, reference.species if reference else '(not in the reference table)',
                      reference.strain if reference else '',
                      '{:.6f}'.format(hit.distance), '{:.3f}'.format(hit.identity),
                      '{}/{}'.format(hit.shared_hashes, hit.total_hashes), '{:.3g}'.format(hit.p_value)],
                     row_class))
    caption = hits_caption(all_hits, hidden)
    out.append(table(['Reference genome', 'Species', 'Strain', 'Distance', 'Identity %', 'Shared hashes',
                      'p-value'], rows, caption=caption))
    return out


def decision_block(result):
    """How the species was decided: what each kind of evidence said, and which one the report followed."""
    call = result.call
    out = ['<h4>The species call</h4>']
    if call is None:
        out.append('<p>No species was called for this sample.</p>')
        return out
    entry = call.entry
    out.append(facts([
        ('Species reported', call.species),
        ('Decided by', call.source),
        ('Species within reach of the closest', ', '.join(call.candidates) if len(call.candidates) > 1
         else (call.candidates[0] if call.candidates else '')),
        ('Sequence type', call.st),
        ('What that ST is, in the isolate database', entry.evidence if entry else ''),
        ('Biovars of its isolates', entry.biovars.replace(';', ', ') if entry and entry.biovars else ''),
        ('Why the ST was not used', call.st_reason),
    ]))
    if call.verdict:
        out.append('<p>{} {}</p>'.format(pill(call.verdict, 'decided'), escape(VERDICTS[call.verdict])))
    return out


def mlst_facts(result):
    """The sequence type and what it was called from."""
    out = ['<h4>MLST: the sequence type</h4>']
    if result.mlst is None:
        out.append('<p>The MLST was not run for this sample.</p>')
        return out
    mlst = result.mlst
    out.append(facts([
        ('Sequence type', mlst.st),
        ('Allele profile', mlst.profile),
        ('Every locus called', 'yes' if mlst.complete else 'no'),
        ('Reads in the input', '{:,}'.format(mlst.total_reads) if mlst.total_reads is not None else ''),
        ('Reads matching an allele', '{:,}'.format(mlst.assigned_reads)),
        ('About the ST', mlst.note),
    ]))
    return out


def mlst_table(result):
    """The nine allele calls, each with the allele it beat and by how much."""
    if result.mlst is None:
        return []
    mlst, out = result.mlst, []
    rows = []
    for call in mlst.calls:
        ahead, ratio = margin(call)
        runner_up = '{}_{} ({})'.format(call.locus, call.runner_up, call.runner_up_reads) \
            if call.runner_up else ''
        rows.append([call.locus, call.allele or UNKNOWN, call.reads or '', runner_up,
                     '' if ahead is None else ahead,
                     '' if ratio is None else '{:.2f}'.format(ratio), call.note])
    out.append(table(['Locus', 'Allele', 'Reads', 'Runner-up (its reads)', 'Ahead by', 'Ratio', 'Note'],
                     rows,
                     caption='The margin each call had over the next allele of its locus; a ratio near 1.00 '
                             'is the one to check.'))
    return out


def inputs_block(files):
    """The files that went in, as they were: path, size, date and checksum."""
    out = ['<h4>Input files</h4>']
    if not files:
        out.append('<p>The input files could not be read when the report was written.</p>')
        return out
    links = any(item.target for item in files)
    headers = ['File', 'Size', 'Modified', 'MD5'] + (['Symbolic link to'] if links else [])
    rows = []
    for item in files:
        row = [item.name, item.size_text, item.modified, item.md5 or 'not computed (--no-md5)']
        rows.append(row + [item.target] if links else row)
    out.append(table(headers, rows, nowrap=(1, 2)))
    out.append('<p class="mono wrap" style="color:var(--muted)">{}</p>'.format(
        escape(' '.join(item.path for item in files))))
    return out


def samples_section(results, references, inputs):
    """One block per sample, in the order of the summary."""
    out = [heading(3, 'Samples', 'samples')]
    for _, members in groups(results):
        for result in members:
            state = status(result)
            out.append('<section class="sample {}" id="{}">'.format(state, anchor(result.sample)))
            out.append('<h3><span class="name">{}</span> {} {}</h3>'.format(
                escape(result.sample), pill(state, state),
                pill('{} ST {}'.format(result.species or UNKNOWN, result.sequence_type), 'decided')))
            if result.error:
                out.append('<div class="panel"><b>This sample failed.</b><p class="wrap">{}</p></div>'.format(
                    escape(result.error)))
            # The fact lists go side by side and the tables across the page: that is what keeps a
            # sample on one printed page, which is how a reviewer reads them
            out.append('<div class="cols"><div>')
            out.extend(identification_facts(result))
            out.append('</div><div>')
            out.extend(decision_block(result))
            out.append('</div></div>')
            out.extend(identification_table(result, references))
            out.append('<div class="cols"><div>')
            out.extend(mlst_facts(result))
            out.append('</div><div>')
            out.extend(inputs_block(inputs.get(result.sample, [])))
            out.append('</div></div>')
            out.extend(mlst_table(result))
            if result.notes:
                out.append('<h4>Notes</h4><ul>')
                out.extend('<li class="wrap">{}</li>'.format(escape(note)) for note in result.notes if note)
                out.append('</ul>')
            out.append('</section>')
    return out


METHODS = [
    ('The references listed for each sample',
     'A sample block lists the five closest reference genomes, and any other near enough to the closest '
     'that the distance does not separate them, which is the evidence the species call rests on. '
     'brucellaid_mash.tsv holds the distance to every reference of the sketch.'),
    ('The species, from a MASH distance',
     'The reads of a sample are reduced to the 1000 smallest hashes of their 21-mers, ignoring the k-mers '
     'seen fewer times than the copy threshold derived from the coverage, and that sketch is compared with '
     'the sketch of each reference genome. The distance is an '
     'estimate of the fraction of positions that differ, so 0.000024 is about 99.998% identity. The closest '
     'reference names the species, provided the reads share at least the required number of its 1000 '
     'min-hashes; below that the sample is reported as ND, because the closest genome of the set is not an '
     'answer. The reference data table above says which set this run measured against.'),
    ('Why there is no biovar in this report',
     'BrucellaID does not report a biovar, because it cannot measure one. Three things stand in the way. '
     'Hardly any reference carries a biovar at all: 7 of the 432 genomes of the default set do, because '
     'the NCBI records of the rest carry none, so 93% of the biovar-labelled B. abortus isolates of the '
     'PubMLST isolate database are of a biovar no reference has. The distance cannot separate the biovars '
     'that are there: the two B. suis biovar 2 genomes are 0.00163 apart, further than the biovar 3 genome '
     'is from one of them (0.00150). And a biovar is a phenotype that genome clustering does not track -- '
     'the best any rule from '
     'sequence type to biovar could manage over the isolate database is 66% for B. abortus and 63% for '
     'B. melitensis, because one ST spans several biovars. The biovars recorded for the ST are shown above '
     'for what they are worth. For a biovar, use classical biotyping, or brucetyper, which reads the '
     'Bruce-ladder and Suis-ladder assays in silico.'),
    ('The sequence type, from the alleles',
     'The same reads are given to Seal with the alleles of the nine loci as references, counting each read '
     'for every allele it shares the highest number of 25-mers with and for no other. The allele the sample '
     'carries is therefore counted for every read of its locus, and a neighbouring allele only for the reads '
     'that fall where the two agree, so the allele with the most reads wins. No mismatch is allowed in a '
     'k-mer, because the alleles of a locus differ by as little as one base. Two alleles with the same count '
     'are reported as a tie rather than separated by a coin toss, and a locus below the minimum read count '
     'is not called.'),
    ('What the ST values mean',
     'A number is a profile that is in the table. "new" means all nine alleles were called but the '
     'combination is not in the table, which is either a genuinely new ST or a novel allele called as its '
     'nearest neighbour. "ND" means at least one locus had no call; the note says which.'),
    ('How the species is decided when two are equally close',
     'Some Brucella species are closer to each other than the variation within either of them: in the '
     'sketches that ship with BrucellaID, the closest pair of different species is always a B. canis and a '
     'B. suis: 0.00024 apart in the 16 genomes of vSNP, 0.000264 in the 432 of the default set. Where a '
     'reference of a second species is within the close-match margin of the closest '
     'one, the distance has not really chosen, and the sequence type chooses among those species, because '
     'the two share no ST. An ST is only used when enough of its isolates agree on a species, and an ST '
     'that contradicts a species the distance picked cleanly is reported as a conflict and not obeyed.'),
    ('What this report does not establish',
     'That the sample is what the label on the tube says, that the library is not a mixture of two '
     'isolates, or that the reference data is current. It establishes what the software was given, what it '
     'measured, and what it concluded from those measurements.'),
]

REFERENCES_CITED = [
    'Whatmore AM, Perrett LL, MacMillan AP. Characterisation of the genetic diversity of Brucella by '
    'multilocus sequencing. BMC Microbiology 2007;7:34. doi:10.1186/1471-2180-7-34',
    'Ondov BD, Treangen TJ, Melsted P, et al. Mash: fast genome and metagenome distance estimation using '
    'MinHash. Genome Biology 2016;17:132. doi:10.1186/s13059-016-0997-x',
    'Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the '
    'PubMLST.org website and their applications. Wellcome Open Research 2018;3:124. '
    'doi:10.12688/wellcomeopenres.14826.1',
    'Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/',
]


def methods_section():
    out = [heading(4, 'Definitions and methods', 'methods')]
    for title, text in METHODS:
        out.append('<h3>{}</h3><p>{}</p>'.format(escape(title), escape(text)))
    return out


def run_section(run, data):
    """The run itself: who, where, with what, against which reference data."""
    out = [heading(5, 'Run information', 'run')]
    software = run.get('software', {})
    out.append(facts([('Operator', run.get('operator', '')),
                      ('User account', run.get('user', '')),
                      ('Computer', run.get('host', '')),
                      ('Operating system', run.get('system', '')),
                      ('Working directory', run.get('working_directory', '')),
                      ('Started', run.get('started', '')),
                      ('Finished', run.get('finished', '')),
                      ('Command', run.get('command', ''))]))
    out.append('<h3>Software</h3>')
    out.append(table(['Program', 'Version'], [[name, value] for name, value in software.items()]))
    parameters = [[name, '' if value is None else value] for name, value in run.get('parameters', {}).items()]
    out.append('<h3>Parameters</h3>')
    out.append(table(['Parameter', 'Value'], parameters))
    out.append('<h3>Reference data</h3>')
    if data:
        rows = [[item['name'], item['detail'], 'the package\'s' if item['packaged'] else 'given with an option',
                 item['md5'] or 'not computed (--no-md5)', item['path']] for item in data]
        out.append(table(['Data', 'What it holds', 'Source', 'MD5', 'Path'], rows,
                         caption='The checksum is what lets a second run, months later, be shown to have '
                                 'used the same reference data as this one.'))
    else:
        out.append('<p>The reference data could not be read when the report was written.</p>')
    out.append('<h3>References</h3><ol>')
    out.extend('<li class="wrap">{}</li>'.format(escape(line)) for line in REFERENCES_CITED)
    out.append('</ol>')
    return out


def build(results, references, run, inputs=None, data=None, note=''):
    """
    The whole report, as one string of HTML.

    :param results: list of pipeline.SampleResult, in the order they were analysed
    :param references: {reference: references.Reference}
    :param run: the run dictionary of report.run_info()
    :param inputs: {sample: [qa.FileFacts]}, from qa.describe_inputs()
    :param data: [reference data facts], from qa.data_facts()
    :param note: a line to put at the top of the report, for a report that is an example or a draft
    """
    inputs = inputs or {}
    counts = counted(results)
    subtitle = '{} sample(s): {} ok, {} with a warning, {} failed. BrucellaID {}, run on {}.'.format(
        len(results), counts[OK], counts[WARNING], counts[FAILED], __version__, run.get('started', ''))
    out = ['<!doctype html>', '<html lang="en">', '<head>', '<meta charset="utf-8">',
           '<meta name="viewport" content="width=device-width, initial-scale=1">',
           '<title>{}</title>'.format(escape(TITLE)),
           '<meta name="generator" content="BrucellaID {}">'.format(escape(__version__)),
           '<style>{}</style>'.format(CSS), '</head>', '<body>', '<main>',
           '<div class="logo">{}</div>'.format(logo()) if logo() else '',
           '<h1>{}</h1>'.format(escape(HEADING if logo() else TITLE)),
           '<p class="sub">{}</p>'.format(escape(subtitle)),
           '<div class="panel" style="border-color:var(--accent);background:var(--accent-soft)">{}</div>'.format(
               escape(note)) if note else '',
           '<nav class="toc"><b>Contents</b><ol>'
           '<li><a href="#summary">Summary</a></li>'
           '<li><a href="#warnings">Notes, warnings and errors</a></li>'
           '<li><a href="#samples">Samples</a></li>'
           '<li><a href="#methods">Definitions and methods</a></li>'
           '<li><a href="#run">Run information</a></li></ol></nav>']
    out += summary_section(results, run)
    out += warnings_section(results)
    out += samples_section(results, references, inputs)
    out += methods_section()
    out += run_section(run, data)
    out.append('<footer>BrucellaID {} &middot; {} &middot; {}@{} &middot; this report is self-contained: it '
               'needs no network to render.</footer>'.format(
                   escape(__version__), escape(run.get('started', '')), escape(run.get('user', '')),
                   escape(run.get('host', ''))))
    out += ['</main>', '</body>', '</html>']
    return '\n'.join(out) + '\n'


def write(results, references, run, path, inputs=None, data=None, note=''):
    """Write the report and return its path."""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(build(results, references, run, inputs=inputs, data=data, note=note))
    return path
