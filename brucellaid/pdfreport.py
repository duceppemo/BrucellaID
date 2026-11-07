"""
The same quality assurance report as a PDF: a file that can be signed, filed and printed unchanged.

The HTML report is the full one and needs nothing installed. This is the version for the places that want a
PDF -- an accreditation file, an attachment to a case record, a paper copy with a signature on it -- and it
holds the same evidence laid out for a page: the summary with the samples grouped and a box to sign, the
warnings, a block per sample with the closest reference genomes and every allele call with its margin, what
the numbers mean, and the run with the checksum of each reference data file.

It is laid out to look like the HTML report rather than merely to contain the same words: the same palette,
the same numbered sections, the same grouped summary, the same two columns of facts above each table, and
one sample to a page. A laboratory that reads one on screen and files the other should not have to work out
twice where to look.

reportlab is an optional dependency ("pip install brucellaid[pdf]"), the way xlsxwriter is for --xlsx, so
that BrucellaID itself still installs with nothing but the standard library.
"""

import re

from . import __version__
from .errors import BrucellaIDError
from .htmlreport import HEADING, METHODS, REFERENCES_CITED, VERDICTS
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

TITLE = 'BrucellaID quality assurance report'
LOGO_FILE = 'logo.png'
LOGO_WIDTH = 62  # Millimetres; the logo is 607 by 122 user units, so about 12.5 mm tall
LOGO_ASPECT = 122 / 607

# The palette of the HTML report, so that the two are recognisably one document.
INK = '#14171a'
SOFT_INK = '#3d454d'
MUTED = '#6b747d'
LINE = '#c2c8ce'
HAIR = '#e2e6ea'
PANEL = '#f4f5f7'
ACCENT = '#15639f'
ACCENT_SOFT = '#eaf2fa'
STATUS_COLOURS = {OK: '#16723c', WARNING: '#8a5000', FAILED: '#a51f1f'}


def _modules():
    """
    reportlab, or a message telling the operator how to get it.

    Imported through a function so that the module can be imported, and its layout tested, without
    reportlab installed.
    """
    try:
        from reportlab.lib.enums import TA_LEFT
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.pdfbase.pdfmetrics import stringWidth
        from reportlab.platypus import (
            Image,
            KeepTogether,
            PageBreak,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
        )
    except ImportError:
        raise BrucellaIDError(
            '--pdf needs the reportlab package: install it with "pip install reportlab" (or "pip install '
            'brucellaid[pdf]"). The HTML report holds the same information and needs nothing '
            'installed.') from None
    return {'TA_LEFT': TA_LEFT, 'A4': A4, 'ParagraphStyle': ParagraphStyle,
            'getSampleStyleSheet': getSampleStyleSheet, 'mm': mm, 'KeepTogether': KeepTogether,
            'PageBreak': PageBreak, 'Paragraph': Paragraph, 'SimpleDocTemplate': SimpleDocTemplate,
            'Spacer': Spacer, 'Table': Table, 'TableStyle': TableStyle, 'Image': Image,
            'stringWidth': stringWidth}


def escape(value):
    """A value as reportlab paragraph text: its mini-XML needs the three XML entities escaped."""
    text = '' if value is None else str(value)
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def fitted_widths(natural, total):
    """
    The width of each column: what its longest line asks for when that all fits the frame, and otherwise
    the widest columns capped at the largest value that does, so that only the columns with the most to
    say wrap. A table is never stretched to the frame, and a narrow column always keeps its natural
    width: spreading the spare space over every column alike is what pushes a value away from the next.
    """
    if sum(natural) <= total:
        return list(natural)
    remaining, count = total, len(natural)
    for width in sorted(natural):
        if width * count < remaining:   # Fits whole under any cap that spends the rest: keep it
            remaining -= width
            count -= 1
    return [min(width, remaining / count) for width in natural]


class Rich(str):
    """A cell value that is paragraph markup already, such as a coloured status, and must not be escaped."""


def coloured(text, colour, bold=True):
    """Text in one of the status colours, for a cell that has to stand out of a quiet table."""
    markup = '<font color="{}">{}</font>'.format(colour, escape(text))
    return Rich('<b>{}</b>'.format(markup) if bold else markup)


class Report:
    """Builds the flowables of the report. One instance per file written."""

    def __init__(self, rl):
        self.rl = rl
        self.story = []
        base = rl['getSampleStyleSheet']()
        mm = rl['mm']
        style = rl['ParagraphStyle']
        self.s = {
            'title': style('t', parent=base['Title'], fontName='Helvetica-Bold', fontSize=15.5, leading=18,
                           alignment=rl['TA_LEFT'], textColor=INK, spaceAfter=1),
            'sub': style('sub', parent=base['Normal'], fontSize=8.4, leading=11, textColor=MUTED,
                         spaceAfter=7),
            'h2': style('h2', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=12, leading=15,
                        textColor=INK, spaceBefore=2, spaceAfter=5),
            'h3': style('h3', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=9.6, leading=12,
                        textColor=INK, spaceBefore=9, spaceAfter=2),
            'h4': style('h4', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=6.8, leading=9,
                        textColor=MUTED, spaceBefore=7, spaceAfter=1),
            'body': style('body', parent=base['Normal'], fontSize=8.4, leading=11, textColor=SOFT_INK,
                          spaceAfter=3),
            'name': style('name', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=11.5,
                          leading=14, textColor=INK, spaceBefore=4, spaceAfter=3),
            'th': style('th', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=6.4, leading=8,
                        textColor=MUTED),
            'cell': style('cell', parent=base['Normal'], fontSize=7.3, leading=9, textColor=SOFT_INK),
            'cellink': style('cellink', parent=base['Normal'], fontSize=7.3, leading=9, textColor=INK),
            'term': style('term', parent=base['Normal'], fontSize=7.3, leading=9.4, textColor=MUTED),
            'value': style('value', parent=base['Normal'], fontSize=7.3, leading=9.4, textColor=INK),
            'group': style('group', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=7.6,
                           leading=9.5, textColor=INK),
            'mono': style('mono', parent=base['Normal'], fontName='Courier', fontSize=6.4, leading=8,
                          textColor=MUTED),
        }
        self.width = rl['A4'][0] - 24 * mm  # The frame of the document, less its margins

    # -- pieces -------------------------------------------------------------------------------------
    def para(self, text, style='body'):
        # reportlab has no text-transform, so the small headings are upper-cased here, as the stylesheet
        # of the HTML report does with CSS
        self.story.append(self.rl['Paragraph'](escape(text.upper() if style == 'h4' else text),
                                               self.s[style]))

    def rich(self, text, style='body'):
        """A paragraph whose markup is already written (a <b> or a <font color>)."""
        self.story.append(self.rl['Paragraph'](text, self.s[style]))

    def space(self, points=4):
        self.story.append(self.rl['Spacer'](1, points))

    def section(self, number, title):
        """A numbered section heading, the number in the accent colour as the HTML report has it."""
        self.rich('<font color="{}">{}</font>   {}'.format(ACCENT, number, escape(title)), 'h2')
        self.rule()

    def rule(self, colour=LINE, thickness=0.7, below=5):
        """A hairline across the frame, under a heading."""
        rl = self.rl
        line = rl['Table']([['']], colWidths=[self.width], rowHeights=[0.1])
        line.setStyle(rl['TableStyle']([('LINEBELOW', (0, 0), (-1, -1), thickness, colour),
                                        ('TOPPADDING', (0, 0), (-1, -1), 0),
                                        ('BOTTOMPADDING', (0, 0), (-1, -1), 0)]))
        self.story.append(line)
        self.space(below)

    def capture(self, build):
        """The flowables one block produces, instead of adding them to the story: for the two columns."""
        saved, self.story = self.story, []
        build()
        captured, self.story = self.story, saved
        return captured

    def columns(self, left, right, gap=9):
        """Two blocks side by side, which is what keeps a sample on one page."""
        rl = self.rl
        half = self.width / 2.0
        table = rl['Table']([[left, right]], colWidths=[half, half])
        table.setStyle(rl['TableStyle']([
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LEFTPADDING', (0, 0), (0, 0), 0), ('RIGHTPADDING', (0, 0), (0, 0), gap),
            ('LEFTPADDING', (1, 0), (1, 0), gap), ('RIGHTPADDING', (1, 0), (1, 0), 0),
            ('TOPPADDING', (0, 0), (-1, -1), 0), ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
        ]))
        self.story.append(table)

    def facts(self, pairs, width=None):
        """A two column list of (term, value): the PDF form of the HTML definition lists."""
        rows = [[self.rl['Paragraph'](escape(term), self.s['term']),
                 self.rl['Paragraph'](escape(value), self.s['value'])]
                for term, value in pairs if value not in (None, '', ())]
        if not rows:
            return
        total = width or self.width
        table = self.rl['Table'](rows, colWidths=[total * 0.42, total * 0.58])
        table.setStyle(self.rl['TableStyle']([
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 1.2), ('BOTTOMPADDING', (0, 0), (-1, -1), 1.2),
            ('LEFTPADDING', (0, 0), (-1, -1), 0), ('RIGHTPADDING', (0, 0), (-1, -1), 2),
        ]))
        self.story.append(table)
        self.space(2)

    def cells(self, values, header=False):
        """One row of a table as paragraphs, leaving markup that is already written alone."""
        style = self.s['th'] if header else self.s['cell']
        return [value if isinstance(value, self.rl['Paragraph'])
                else self.rl['Paragraph'](value if isinstance(value, Rich) else escape(value), style)
                for value in values]

    def cell_text(self, value):
        """The text a cell will show, for measuring it: markup stripped and the XML entities restored."""
        if isinstance(value, self.rl['Paragraph']):
            return value.getPlainText()
        text = str(value)
        if isinstance(value, Rich):
            text = re.sub(r'<[^>]+>', '', text)
            text = text.replace('&lt;', '<').replace('&gt;', '>').replace('&amp;', '&')
        return text

    def natural_widths(self, headers, rows):
        """What each column asks for: its longest value on one line, plus padding and room for bold."""
        measure = self.rl['stringWidth']
        widths = [measure(header.upper(), self.s['th'].fontName, self.s['th'].fontSize)
                  for header in headers]
        for row in rows:
            for column, value in enumerate(row):
                widths[column] = max(widths[column], measure(self.cell_text(value),
                                                             self.s['cell'].fontName,
                                                             self.s['cell'].fontSize))
        return [width + 8 for width in widths]

    def table(self, headers, rows, highlight=None, keep=True, width=None, groups_at=()):
        """
        A quiet table: no grid, a rule above and below the table and under the header, and a hairline
        between the rows. Each column is
        as wide as its content, capped by fitted_widths() when the whole does not fit the frame.

        :param rows: list of rows of values; a Rich value is written as markup
        :param highlight: {row index: colour} for the rows a reader should land on
        :param groups_at: {row index: label} rows that name a group and span the table
        """
        rl = self.rl
        data = [self.cells([h.upper() for h in headers], header=True)]
        for index, row in enumerate(rows):
            if index in (groups_at or {}):
                data.append([rl['Paragraph'](escape(groups_at[index]), self.s['group'])]
                            + [''] * (len(headers) - 1))
            data.append(self.cells(row))
        columns = fitted_widths(self.natural_widths(headers, rows), width or self.width)
        table = rl['Table'](data, colWidths=columns, repeatRows=1, hAlign='LEFT')
        commands = [
            ('LINEABOVE', (0, 0), (-1, 0), 0.7, LINE),
            ('LINEBELOW', (0, 0), (-1, 0), 0.7, LINE),
            ('LINEBELOW', (0, 1), (-1, -2), 0.25, HAIR),
            ('LINEBELOW', (0, -1), (-1, -1), 0.7, LINE),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 2.4), ('BOTTOMPADDING', (0, 0), (-1, -1), 2.4),
            ('LEFTPADDING', (0, 0), (-1, -1), 3), ('RIGHTPADDING', (0, 0), (-1, -1), 3),
        ]
        # The group rows, and the highlighted rows, in the data's own numbering
        offset = 1
        for index in range(len(rows)):
            if index in (groups_at or {}):
                commands += [('SPAN', (0, offset), (-1, offset)),
                             ('BACKGROUND', (0, offset), (-1, offset), PANEL),
                             ('LINEBELOW', (0, offset), (-1, offset), 0.25, LINE)]
                offset += 1
            if index in (highlight or {}):
                commands.append(('BACKGROUND', (0, offset), (-1, offset), highlight[index]))
            offset += 1
        table.setStyle(rl['TableStyle'](commands))
        self.story.append(rl['KeepTogether'](table) if keep else table)
        self.space(3)

    def sign_off(self):
        """A box for the reviewer: name, date and signature, with room to write."""
        rl = self.rl
        # The panel adds 9 points of padding on each side, so the lines inside it are that much
        # narrower; without this they run past the edge of the box they sit in
        inner = self.width - 20
        third = inner / 3.0
        rows = [['', '', ''],
                [rl['Paragraph']('Name', self.s['term']), rl['Paragraph']('Date', self.s['term']),
                 rl['Paragraph']('Signature', self.s['term'])]]
        table = rl['Table'](rows, colWidths=[third] * 3, rowHeights=[13 * rl['mm'], None])
        table.setStyle(rl['TableStyle']([
            ('LINEBELOW', (0, 0), (-1, 0), 0.6, SOFT_INK),
            ('VALIGN', (0, 0), (-1, -1), 'BOTTOM'),
            ('LEFTPADDING', (0, 0), (-1, -1), 0), ('RIGHTPADDING', (0, 0), (-1, -1), 10),
        ]))
        heading = rl['Paragraph']('REVIEWED BY', self.s['h4'])
        panel = rl['Table']([[heading], [table]], colWidths=[inner + 20])
        panel.setStyle(rl['TableStyle']([
            ('BACKGROUND', (0, 0), (-1, -1), PANEL),
            ('BOX', (0, 0), (-1, -1), 0.5, LINE),
            ('LEFTPADDING', (0, 0), (-1, -1), 9), ('RIGHTPADDING', (0, 0), (-1, -1), 9),
            ('TOPPADDING', (0, 0), (-1, 0), 7), ('BOTTOMPADDING', (0, -1), (-1, -1), 7),
        ]))
        self.story.append(rl['KeepTogether'](panel))

    # -- sections -----------------------------------------------------------------------------------
    def masthead(self, results, run, note=''):
        """The logo, the title and what the run came to."""
        try:
            path = data_path(LOGO_FILE)
        except BrucellaIDError:
            path = None  # A package without the image still writes the report
        if path:
            mm = self.rl['mm']
            self.story.append(self.rl['Image'](path, width=LOGO_WIDTH * mm,
                                               height=LOGO_WIDTH * LOGO_ASPECT * mm, hAlign='LEFT'))
            self.space(5)
        counts = counted(results)
        self.para(HEADING, 'title')
        self.para('{} sample(s): {} ok, {} with a warning, {} failed. BrucellaID {}, run on {}.'.format(
            len(results), counts[OK], counts[WARNING], counts[FAILED], __version__,
            run.get('started', '')), 'sub')
        if note:
            self.rich('<b>{}</b>'.format(escape(note)))
            self.space(4)

    def summary(self, results, run):
        """One table, the samples grouped by what they were called, and a box to sign."""
        self.section(1, 'Summary')
        self.facts([('Run started', run.get('started', '')), ('Run finished', run.get('finished', '')),
                    ('Operator', run.get('operator', '')), ('Command', run.get('command', ''))])
        self.space(4)
        rows, labels, highlight = [], {}, {}
        for (species, st), members in groups(results):
            labels[len(rows)] = '{} · ST {}   ({} sample(s))'.format(species, st, len(members)) \
                if species else 'Failed   ({} sample(s))'.format(len(members))
            for result in members:
                state = status(result)
                if state != OK:
                    highlight[len(rows)] = ACCENT_SOFT
                rows.append([result.sample, result.species or UNKNOWN,
                             result.call.source if result.call else '',
                             result.call.verdict if result.call and result.call.verdict else '',
                             result.sequence_type, result.mlst.profile if result.mlst else '',
                             result.closest_reference, coloured(state, STATUS_COLOURS[state])])
        self.table(['Sample', 'Species', 'From', 'Decision', 'ST', 'Allele profile', 'Closest reference',
                    'Status'], rows, groups_at=labels,
                   highlight=highlight, keep=False)
        self.space(6)
        self.sign_off()

    def warnings(self, results):
        self.story.append(self.rl['PageBreak']())
        self.section(2, 'Notes, warnings and errors')
        rows = [[sample, coloured(state, STATUS_COLOURS[state]), text]
                for sample, state, text in warnings(results)]
        if not rows:
            self.para('None: every sample was identified and typed without a note.')
            return
        self.table(['Sample', 'Status', 'What it says'], rows, keep=False)

    # -- one sample ---------------------------------------------------------------------------------
    def identification_facts(self, result, width):
        self.para('Identification: the MASH distance', 'h4')
        if result.identification is None:
            self.para('MASH did not run for this sample.')
            return
        best = result.identification.best
        self.facts([
            ('Closest reference', result.closest_reference),
            ('What it is', result.reference.description if result.reference
             else 'not in the reference table'),
            ('NCBI accessions', result.reference.accessions.replace(';', ', ')
             if result.reference and result.reference.accessions else ''),
            ('MASH distance', '{:.6f}'.format(best.distance)),
            ('Estimated identity', '{:.3f}%'.format(best.identity)),
            ('Shared min-hashes', '{}/{} (at least {} required)'.format(
                best.shared_hashes, best.total_hashes, result.identification.min_hashes)),
            ('K-mer copies required', kmer_copies_note(result.identification)),
        ], width=width)

    def species_call(self, result, width):
        self.para('The species call', 'h4')
        call = result.call
        if call is None:
            self.para('No species was called for this sample.')
            return
        self.facts([
            ('Species reported', call.species),
            ('Decided by', call.source),
            ('Species within reach', ', '.join(call.candidates)),
            ('Sequence type', call.st),
            ('What that ST is', call.entry.evidence if call.entry else ''),
            ('Biovars of its isolates', call.entry.biovars.replace(';', ', ')
             if call.entry and call.entry.biovars else ''),
            ('Why the ST was not used', call.st_reason),
        ], width=width)
        if call.verdict:
            self.rich('<font color="{}"><b>{}</b></font>  {}'.format(
                ACCENT, escape(call.verdict), escape(VERDICTS[call.verdict])), 'cell')

    def identification_table(self, result, references):
        if result.identification is None:
            return
        best = result.identification.best
        all_hits = close_hits(result, references)
        shown, hidden = top_hits(all_hits)
        rows, highlight = [], {}
        for index, (hit, reference, close) in enumerate(shown):
            if hit is best or close:
                highlight[index] = ACCENT_SOFT if hit is best else PANEL
            rows.append([hit.reference, reference.species if reference else '(not in the table)',
                         reference.strain if reference else '',
                         '{:.6f}'.format(hit.distance), '{:.3f}'.format(hit.identity),
                         '{}/{}'.format(hit.shared_hashes, hit.total_hashes)])
        caption = hits_caption(all_hits, hidden)
        self.para(caption, 'cell')
        self.table(['Reference genome', 'Species', 'Strain', 'Distance', 'Identity %', 'Shared hashes'],
                   rows, highlight=highlight, keep=False)

    def mlst_facts(self, result, width):
        self.para('MLST: the sequence type', 'h4')
        if result.mlst is None:
            self.para('The MLST was not run for this sample.')
            return
        mlst = result.mlst
        self.facts([('Sequence type', mlst.st), ('Allele profile', mlst.profile),
                    ('Every locus called', 'yes' if mlst.complete else 'no'),
                    ('Reads in the input', '{:,}'.format(mlst.total_reads)
                     if mlst.total_reads is not None else ''),
                    ('Reads matching an allele', '{:,}'.format(mlst.assigned_reads)),
                    ('About the ST', mlst.note)], width=width)

    def mlst_table(self, result):
        if result.mlst is None:
            return
        rows = []
        for call in result.mlst.calls:
            ahead, ratio = margin(call)
            runner_up = '{}_{} ({})'.format(call.locus, call.runner_up, call.runner_up_reads) \
                if call.runner_up else ''
            rows.append([call.locus, call.allele or UNKNOWN, call.reads or '', runner_up,
                         '' if ahead is None else ahead,
                         '' if ratio is None else '{:.2f}'.format(ratio), call.note])
        self.para('The margin each call had over the next allele of its locus; a ratio near 1.00 is the one '
                  'to check.', 'cell')
        self.table(['Locus', 'Allele', 'Reads', 'Runner-up (its reads)', 'Ahead', 'Ratio', 'Note'],
                   rows, keep=False)

    def inputs(self, files, width):
        self.para('Input files', 'h4')
        if not files:
            self.para('The input files could not be read when the report was written.')
            return
        # Half the frame cannot hold a checksum column as well, so the checksum and the path go on the
        # lines below, where a long value has the width to be read: file, then what it really is.
        rows = [[item.name, item.size_text, item.modified] for item in files]
        self.table(['File', 'Size', 'Modified'], rows, keep=False, width=width)
        for item in files:
            self.rich(escape(item.path) + (' -&gt; {}'.format(escape(item.target)) if item.target else ''),
                      'mono')
            self.rich('MD5 {}'.format(escape(item.md5 or 'not computed (--no-md5)')), 'mono')

    def sample(self, result, references, inputs, first=False):
        # The first sample follows the section heading on the same page; the rest each start one
        if not first:
            self.story.append(self.rl['PageBreak']())
        state = status(result)
        self.rich('{}  <font size="7.5" color="{}"><b>{}</b></font>  <font size="7.5" color="{}"><b>{} '
                  'ST {}</b></font>'.format(escape(result.sample), STATUS_COLOURS[state], state, ACCENT,
                                            escape(result.species or UNKNOWN),
                                            escape(result.sequence_type)), 'name')
        self.rule(colour=HAIR, thickness=0.5, below=3)
        if result.error:
            self.rich('<b>This sample failed.</b> {}'.format(escape(result.error)))

        half = self.width / 2.0 - 9
        self.columns(self.capture(lambda: self.identification_facts(result, half)),
                     self.capture(lambda: self.species_call(result, half)))
        self.space(5)
        self.identification_table(result, references)
        self.space(3)
        self.columns(self.capture(lambda: self.mlst_facts(result, half)),
                     self.capture(lambda: self.inputs(inputs.get(result.sample, []), half)))
        self.space(5)
        self.mlst_table(result)
        if result.notes:
            self.para('Notes', 'h4')
            for note in result.notes:
                if note:
                    self.rich('- {}'.format(escape(note)), 'cell')

    # -- the rest -----------------------------------------------------------------------------------
    def methods(self):
        self.story.append(self.rl['PageBreak']())
        self.section(4, 'Definitions and methods')
        for title, text in METHODS:
            self.para(title, 'h3')
            self.para(text)

    def run_information(self, run, data):
        self.story.append(self.rl['PageBreak']())
        self.section(5, 'Run information')
        self.facts([('Operator', run.get('operator', '')), ('User account', run.get('user', '')),
                    ('Computer', run.get('host', '')), ('Operating system', run.get('system', '')),
                    ('Working directory', run.get('working_directory', '')),
                    ('Started', run.get('started', '')), ('Finished', run.get('finished', '')),
                    ('Command', run.get('command', ''))])
        self.para('Software', 'h3')
        self.table(['Program', 'Version'],
                   [[name, value] for name, value in run.get('software', {}).items()])
        self.para('Parameters', 'h3')
        self.table(['Parameter', 'Value'],
                   [[name, '' if value is None else value] for name, value in run.get('parameters', {}).items()],
                   keep=False)
        self.para('Reference data', 'h3')
        if data:
            self.table(['Data', 'What it holds', 'Source', 'MD5'],
                       [[item['name'], item['detail'],
                         'the package\'s' if item['packaged'] else 'given with an option',
                         item['md5'] or 'not computed'] for item in data],
                       keep=False)
            for item in data:
                self.rich(escape(item['path']), 'mono')
        else:
            self.para('The reference data could not be read when the report was written.')
        self.para('References', 'h3')
        for index, line in enumerate(REFERENCES_CITED, start=1):
            self.rich('{}. {}'.format(index, escape(line)), 'cell')


def numbered_canvas(run, margin, baseline):
    """
    A canvas that stamps "Page 2 of 9" on every page.

    The total is only known once the document is finished, so each page is held back as it is completed and
    the footer is drawn on all of them at the end. A page of a quality assurance report that leaves the
    stapler on its own has to say which report it came from and whether anything is missing.
    """
    from reportlab.pdfgen.canvas import Canvas

    left = 'BrucellaID {}  |  {}  |  {}@{}'.format(__version__, run.get('started', ''), run.get('user', ''),
                                                   run.get('host', ''))[:150]

    class NumberedCanvas(Canvas):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._pages = []

        def showPage(self):
            self._pages.append(self.__dict__.copy())
            self._startPage()

        def save(self):
            total = len(self._pages)
            for state in self._pages:
                self.__dict__.update(state)
                self.saveState()
                self.setStrokeColor(HAIR)
                self.setLineWidth(0.5)
                self.line(margin, baseline + 7, self._pagesize[0] - margin, baseline + 7)
                self.setFont('Helvetica', 6.3)
                self.setFillColor(MUTED)
                self.drawString(margin, baseline, left)
                self.drawRightString(self._pagesize[0] - margin, baseline,
                                     'Page {} of {}'.format(self._pageNumber, total))
                self.restoreState()
                super().showPage()
            super().save()

    return NumberedCanvas


def write(results, references, run, path, inputs=None, data=None, note=''):
    """
    Write the PDF report and return its path.

    :raises BrucellaIDError: when reportlab is not installed
    """
    rl = _modules()
    inputs = inputs or {}
    mm = rl['mm']
    document = rl['SimpleDocTemplate'](
        str(path), pagesize=rl['A4'], leftMargin=12 * mm, rightMargin=12 * mm, topMargin=13 * mm,
        bottomMargin=15 * mm, title=TITLE, author=run.get('operator', ''), subject='BrucellaID run',
        creator='BrucellaID {}'.format(__version__))
    report = Report(rl)
    report.masthead(results, run, note=note)
    report.summary(results, run)
    report.warnings(results)
    report.story.append(rl['PageBreak']())
    report.section(3, 'Samples')
    ordered = [result for _, members in groups(results) for result in members]
    for index, result in enumerate(ordered):
        report.sample(result, references, inputs, first=index == 0)
    report.methods()
    report.run_information(run, data)
    document.build(report.story, canvasmaker=numbered_canvas(run, 12 * mm, 9 * mm))
    return path
