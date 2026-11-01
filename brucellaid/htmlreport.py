"""
The quality assurance report: one self-contained HTML file holding everything needed to check a run.

Self-contained is the whole point. No stylesheet, no script, no font and no image is fetched from anywhere,
so the file can be attached to an email, dropped into a laboratory information system or archived for seven
years and still render exactly as it did on the day it was written. It prints to clean pages as it stands.

What it holds, in the order a reviewer reads it: a summary of the run with the samples grouped by what they
were called and a box to sign; every warning and error, one row each; then, per sample, the identification
with *every* distance MASH measured rather than only the winning one, how the species was decided, the nine
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
from .mlst import UNKNOWN
from .qa import FAILED, OK, WARNING, close_hits, counted, groups, margin, status, warnings
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
  --bg: #ffffff; --panel: #f7f8fa; --ink: #16191d; --muted: #5c6670; --line: #dfe3e8;
  --accent: #1a6ab3; --accent-soft: #e8f1fa;
  --ok: #1a7f44; --ok-soft: #e7f5ec; --warn: #9a5b00; --warn-soft: #fdf2e2;
  --bad: #ad2525; --bad-soft: #fdecec; --mono: ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace;
  /* The logo of the package writes its four colours as these variables, so one file serves both themes */
  --logo-ink: #2f2a3d; --logo-accent: #7d4a9e; --logo-reference: #9aa5ad; --logo-muted: #5d6a72;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #15181c; --panel: #1d2127; --ink: #e8eaed; --muted: #9aa3ad; --line: #2e343c;
    --accent: #6fb3f2; --accent-soft: #19283a;
    --ok: #67c48d; --ok-soft: #14271c; --warn: #e0a95c; --warn-soft: #2a2114;
    --bad: #f08b8b; --bad-soft: #2b1818;
    --logo-ink: #f0edf5; --logo-accent: #b185d8; --logo-reference: #6d787f;
    --logo-muted: #a9b6bd;
  }
}
:root[data-theme="dark"] {
  --bg: #15181c; --panel: #1d2127; --ink: #e8eaed; --muted: #9aa3ad; --line: #2e343c;
  --accent: #6fb3f2; --accent-soft: #19283a;
  --ok: #67c48d; --ok-soft: #14271c; --warn: #e0a95c; --warn-soft: #2a2114;
  --bad: #f08b8b; --bad-soft: #2b1818;
  --logo-ink: #f0edf5; --logo-accent: #b185d8; --logo-reference: #6d787f;
  --logo-muted: #a9b6bd;
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 0 16px 64px; background: var(--bg); color: var(--ink);
  font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
}
main { max-width: 1100px; margin: 0 auto; }
.logo { margin: 26px 0 2px; }
.logo svg { width: 360px; max-width: 100%; height: auto; display: block; }
h1 { font-size: 1.6rem; margin: 10px 0 4px; letter-spacing: -0.01em; }
h2 {
  font-size: 1.2rem; margin: 40px 0 12px; padding-bottom: 6px; border-bottom: 2px solid var(--line);
}
h3 { font-size: 1.02rem; margin: 22px 0 8px; }
h4 { font-size: 0.92rem; margin: 18px 0 6px; color: var(--muted); text-transform: uppercase;
     letter-spacing: 0.06em; }
p, ul, ol { margin: 8px 0; }
a { color: var(--accent); }
code, .mono { font-family: var(--mono); font-size: 0.88em; }
.sub { color: var(--muted); margin: 0 0 18px; }
.chips { display: flex; flex-wrap: wrap; gap: 8px; margin: 14px 0 6px; }
.chip {
  display: inline-flex; align-items: baseline; gap: 6px; padding: 4px 10px; border-radius: 999px;
  background: var(--panel); border: 1px solid var(--line); font-size: 0.85rem;
}
.chip b { font-size: 0.95rem; }
.chip.ok { background: var(--ok-soft); border-color: var(--ok); color: var(--ok); }
.chip.warning { background: var(--warn-soft); border-color: var(--warn); color: var(--warn); }
.chip.failed { background: var(--bad-soft); border-color: var(--bad); color: var(--bad); }
.tag {
  display: inline-block; padding: 1px 7px; border-radius: 4px; font-size: 0.78rem; font-weight: 600;
  border: 1px solid var(--line); background: var(--panel); white-space: nowrap;
}
.tag.ok { background: var(--ok-soft); border-color: var(--ok); color: var(--ok); }
.tag.warning { background: var(--warn-soft); border-color: var(--warn); color: var(--warn); }
.tag.failed { background: var(--bad-soft); border-color: var(--bad); color: var(--bad); }
.tag.decided { background: var(--accent-soft); border-color: var(--accent); color: var(--accent); }
table { border-collapse: collapse; width: 100%; margin: 10px 0 4px; font-size: 0.88rem; }
caption { text-align: left; color: var(--muted); font-size: 0.85rem; padding-bottom: 6px; }
th, td { text-align: left; padding: 5px 9px; border-bottom: 1px solid var(--line); vertical-align: top; }
th { background: var(--panel); font-weight: 600; white-space: nowrap; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.best td { background: var(--accent-soft); }
tr.close td { background: var(--panel); }
tbody tr:last-child td { border-bottom: 1px solid var(--line); }
.group { margin: 0 0 4px; padding: 6px 10px; background: var(--accent-soft); border-left: 3px solid
         var(--accent); font-size: 0.88rem; }
.panel { background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 12px 14px;
         margin: 12px 0; }
.sample { border: 1px solid var(--line); border-radius: 10px; padding: 4px 16px 16px; margin: 18px 0;
          background: var(--bg); }
.sample > h3 { display: flex; flex-wrap: wrap; gap: 10px; align-items: baseline; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; margin: 8px 0; }
dl.facts dt { color: var(--muted); font-size: 0.85rem; }
dl.facts dd { margin: 0; }
.sign { display: grid; grid-template-columns: repeat(3, 1fr); gap: 18px; margin-top: 12px; }
.sign div { border-top: 1px solid var(--ink); padding-top: 4px; font-size: 0.82rem; color: var(--muted); }
.sign .line { height: 42px; border: 0; }
nav.toc ol { margin: 4px 0; padding-left: 22px; }
nav.toc { font-size: 0.9rem; }
footer { margin-top: 40px; padding-top: 10px; border-top: 1px solid var(--line); color: var(--muted);
         font-size: 0.8rem; }
.wrap { overflow-wrap: anywhere; }
@media (max-width: 640px) {
  dl.facts { grid-template-columns: 1fr; }
  dl.facts dd { margin-bottom: 6px; }
  .sign { grid-template-columns: 1fr; }
  table { font-size: 0.8rem; }
}
@media print {
  @page { size: A4; margin: 14mm 12mm 16mm; }
  body { padding: 0; font-size: 10pt; background: #fff; color: #000; }
  :root { --bg: #fff; --panel: #f2f2f2; --ink: #000; --muted: #444; --line: #bbb; --accent: #0b4f8a;
          --accent-soft: #eaf1f8; }
  nav.toc, .noprint { display: none; }
  .logo svg { width: 230px; }
  h2 { break-before: page; }
  h2.first { break-before: avoid; }
  /* A sample block is often taller than a page, so it may not be kept whole: keeping its tables whole
     is what matters, and a heading must not be left alone at the foot of a page. */
  table, .panel { break-inside: avoid; }
  h2, h3, h4 { break-after: avoid; }
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


def tag(text, kind=''):
    """A coloured label: a status, or a verdict."""
    return '<span class="tag {}">{}</span>'.format(escape(kind), escape(text)) if kind else \
        '<span class="tag">{}</span>'.format(escape(text))


def cell(value, numeric=False):
    return '<td class="num">{}</td>'.format(escape(value)) if numeric else '<td>{}</td>'.format(escape(value))


def table(headers, rows, caption='', numeric=(), classes=()):
    """
    A table from already escaped-on-write values.

    :param headers: column titles
    :param rows: list of rows of values, or (values, row class) pairs when a row is highlighted
    :param numeric: indexes of the columns to align right
    """
    out = ['<table class="{}">'.format(' '.join(classes)) if classes else '<table>']
    if caption:
        out.append('<caption>{}</caption>'.format(escape(caption)))
    head = ''.join('<th class="num">{}</th>'.format(escape(h)) if i in numeric else
                   '<th>{}</th>'.format(escape(h)) for i, h in enumerate(headers))
    out.append('<thead><tr>{}</tr></thead>'.format(head))
    out.append('<tbody>')
    for row in rows:
        values, row_class = row if isinstance(row, tuple) else (row, '')
        cells = ''.join(cell(v, i in numeric) for i, v in enumerate(values))
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
    out = ['<h2 class="first" id="summary">Summary</h2>',
           '<div class="chips">{}</div>'.format(''.join(chips)),
           facts([('Run started', run.get('started', '')),
                  ('Run finished', run.get('finished', '')),
                  ('Operator', run.get('operator', '')),
                  ('Command', run.get('command', ''))])]
    for (species, st), members in groups(results):
        if species:
            out.append('<p class="group"><b>{}</b>{} &mdash; ST {} &mdash; {} sample(s)</p>'.format(
                escape(species), '', escape(st), len(members)))
        else:
            out.append('<p class="group" style="border-color:var(--bad)"><b>Failed</b> &mdash; '
                       '{} sample(s)</p>'.format(len(members)))
        rows = []
        for result in members:
            state = status(result)
            rows.append(([result.sample, result.species or UNKNOWN, result.biovar,
                          result.call.source if result.call else '',
                          result.call.verdict if result.call and result.call.verdict else '',
                          result.sequence_type, result.mlst.profile if result.mlst else '',
                          result.closest_reference, state], state if state != OK else ''))
        out.append(table(['Sample', 'Species', 'Biovar', 'Species from', 'Decision', 'ST', 'Allele profile',
                          'Closest reference', 'Status'], rows))
    out.append('<div class="panel"><h4>Reviewed by</h4><div class="sign">'
               '<div class="line"></div><div class="line"></div><div class="line"></div>'
               '<div>Name</div><div>Date</div><div>Signature</div></div></div>')
    return out


def warnings_section(results):
    """Every warning and error of the run, one row each."""
    rows = [([sample, state, text], state if state != OK else '')
            for sample, state, text in warnings(results)]
    out = ['<h2 id="warnings">Notes, warnings and errors</h2>']
    if not rows:
        out.append('<p>None: every sample was identified and typed without a note.</p>')
        return out
    out.append(table(['Sample', 'Status', 'What it says'], rows,
                     caption='Errors first, then the notes of each sample, in the order they were raised. '
                             'A note on a sample whose status is "ok" is a remark, not a problem.'))
    return out


def identification_block(result, references):
    """Every distance MASH measured for one sample, not only the winning one."""
    out = ['<h4>Identification: the MASH distances</h4>']
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
        ('p-value', '{:.3g}'.format(best.p_value)),
    ]))
    rows = []
    for hit, reference, close in close_hits(result, references):
        row_class = 'best' if hit is best else ('close' if close else '')
        rows.append(([hit.reference, reference.species if reference else '(not in the reference table)',
                      reference.biovar if reference else '', reference.strain if reference else '',
                      '{:.6f}'.format(hit.distance), '{:.3f}'.format(hit.identity),
                      '{}/{}'.format(hit.shared_hashes, hit.total_hashes), '{:.3g}'.format(hit.p_value)],
                     row_class))
    out.append(table(['Reference genome', 'Species', 'Biovar', 'Strain', 'Distance', 'Identity %',
                      'Shared hashes', 'p-value'], rows, numeric=(4, 5, 6, 7),
                     caption='Closest first. The closest is highlighted; a shaded row is near enough to the '
                             'closest that the distance alone does not separate the two.'))
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
        ('Biovar', '{} (of {}, the closest reference of that species)'.format(call.biovar, call.reference.strain)
         if call.biovar and call.reference else call.biovar),
        ('Decided by', call.source),
        ('Species within reach of the closest', ', '.join(call.candidates) if len(call.candidates) > 1
         else (call.candidates[0] if call.candidates else '')),
        ('Sequence type', call.st),
        ('What that ST is, in the isolate database', entry.evidence if entry else ''),
        ('Biovars recorded for it', entry.biovars.replace(';', ', ') if entry and entry.biovars else ''),
        ('Why the ST was not used', call.st_reason),
    ]))
    if call.verdict:
        out.append('<p>{} {}</p>'.format(tag(call.verdict, 'decided'), escape(VERDICTS[call.verdict])))
    return out


def mlst_block(result):
    """The nine allele calls, each with the allele it beat and by how much."""
    out = ['<h4>MLST: the nine allele calls</h4>']
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
    rows = []
    for call in mlst.calls:
        ahead, ratio = margin(call)
        rows.append([call.locus, call.allele or UNKNOWN, call.reads or '',
                     '{}_{}'.format(call.locus, call.runner_up) if call.runner_up else '',
                     call.runner_up_reads or '', '' if ahead is None else ahead,
                     '' if ratio is None else '{:.2f}'.format(ratio), call.note])
    out.append(table(['Locus', 'Allele', 'Reads', 'Runner-up', 'Its reads', 'Ahead by', 'Ratio', 'Note'],
                     rows, numeric=(2, 4, 5, 6),
                     caption='"Ahead by" and "Ratio" are the margin the call had over the next allele of the '
                             'locus. A ratio near 1.00 is the call to check first.'))
    return out


def inputs_block(files):
    """The files that went in, as they were: path, size, date and checksum."""
    out = ['<h4>Input files</h4>']
    if not files:
        out.append('<p>The input files could not be read when the report was written.</p>')
        return out
    rows = []
    for item in files:
        rows.append([item.name, item.size_text, item.modified, item.md5 or 'not computed (--no-md5)',
                     item.target])
    out.append(table(['File', 'Size', 'Modified', 'MD5', 'Symbolic link to'], rows, numeric=(1,)))
    out.append('<p class="mono wrap" style="color:var(--muted)">{}</p>'.format(
        escape(' '.join(item.path for item in files))))
    return out


def samples_section(results, references, inputs):
    """One block per sample, in the order of the summary."""
    out = ['<h2 id="samples">Samples</h2>']
    for _, members in groups(results):
        for result in members:
            state = status(result)
            out.append('<section class="sample" id="{}">'.format(anchor(result.sample)))
            out.append('<h3>{} {} {}</h3>'.format(
                escape(result.sample), tag(state, state),
                tag('{} ST {}'.format(result.species or UNKNOWN, result.sequence_type), 'decided')))
            if result.error:
                out.append('<div class="panel"><b>This sample failed.</b><p class="wrap">{}</p></div>'.format(
                    escape(result.error)))
            out.extend(identification_block(result, references))
            out.extend(decision_block(result))
            out.extend(mlst_block(result))
            out.extend(inputs_block(inputs.get(result.sample, [])))
            if result.notes:
                out.append('<h4>Notes</h4><ul>')
                out.extend('<li class="wrap">{}</li>'.format(escape(note)) for note in result.notes if note)
                out.append('</ul>')
            out.append('</section>')
    return out


METHODS = [
    ('The species, from a MASH distance',
     'The reads of a sample are reduced to the 1000 smallest hashes of their 21-mers, ignoring the k-mers '
     'seen once, and that sketch is compared with the sketch of each reference genome. The distance is an '
     'estimate of the fraction of positions that differ, so 0.000024 is about 99.998% identity. The closest '
     'reference names the species, provided the reads share at least the required number of its 1000 '
     'min-hashes; below that the sample is reported as ND, because the closest of sixteen genomes is not an '
     'answer.'),
    ('Why the biovar is not a biovar call',
     'The biovar in this report is the biovar of the closest reference genome of the species reported. The '
     'Brucella biovar references differ by very little, so the distance table in each sample block is the '
     'thing to read before quoting a biovar. For a biovar, use classical biotyping or an in silico '
     'Bruce-ladder.'),
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
     'sketch that ships with BrucellaID, B. canis and B. suis bv. 4 are 0.00024 apart, the closest pair of '
     'the sixteen. Where a reference of a second species is within the close-match margin of the closest '
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
    out = ['<h2 id="methods">Definitions and methods</h2>']
    for heading, text in METHODS:
        out.append('<h3>{}</h3><p>{}</p>'.format(escape(heading), escape(text)))
    return out


def run_section(run, data):
    """The run itself: who, where, with what, against which reference data."""
    out = ['<h2 id="run">Run information</h2>']
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
