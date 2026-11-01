"""
The same quality assurance report as a PDF: a file that can be signed, filed and printed unchanged.

The HTML report is the full one and needs nothing installed. This is the version for the places that want a
PDF -- an accreditation file, an attachment to a case record, a paper copy with a signature on it -- and it
holds the same evidence laid out for a page: the summary with the samples grouped and a box to sign, the
warnings, a block per sample with every distance MASH measured and every allele call with its margin, what
the numbers mean, and the run with the checksum of each reference data file.

reportlab is an optional dependency ("pip install brucellaid[pdf]"), the way xlsxwriter is for --xlsx, so
that BrucellaID itself still installs with nothing but the standard library.
"""

from . import __version__
from .errors import BrucellaIDError
from .htmlreport import METHODS, REFERENCES_CITED, VERDICTS
from .mlst import UNKNOWN
from .qa import FAILED, OK, WARNING, close_hits, counted, groups, margin, status, warnings

TITLE = 'BrucellaID quality assurance report'
# Grey for the rules and the header rows, blue for the row a reader should land on first, and the three
# status colours of the HTML report, so that the two look like the same document.
INK = '#16191d'
MUTED = '#5c6670'
LINE = '#c9ced5'
PANEL = '#eef0f3'
ACCENT = '#1a6ab3'
ACCENT_SOFT = '#e8f1fa'
STATUS_COLOURS = {OK: '#1a7f44', WARNING: '#9a5b00', FAILED: '#ad2525'}


def _modules():
    """
    reportlab, or a message telling the operator how to get it.

    Imported through a function so that the module can be imported, and its layout tested, without
    reportlab installed.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_LEFT
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.platypus import (
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
    return {'colors': colors, 'TA_LEFT': TA_LEFT, 'A4': A4, 'ParagraphStyle': ParagraphStyle,
            'getSampleStyleSheet': getSampleStyleSheet, 'mm': mm, 'KeepTogether': KeepTogether,
            'PageBreak': PageBreak, 'Paragraph': Paragraph, 'SimpleDocTemplate': SimpleDocTemplate,
            'Spacer': Spacer, 'Table': Table, 'TableStyle': TableStyle}


def escape(value):
    """A value as reportlab paragraph text: its mini-XML needs the three XML entities escaped."""
    text = '' if value is None else str(value)
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


class Report:
    """Builds the flowables of the report. One instance per file written."""

    def __init__(self, rl):
        self.rl = rl
        self.story = []
        base = rl['getSampleStyleSheet']()
        mm = rl['mm']
        style = rl['ParagraphStyle']
        self.s = {
            'title': style('t', parent=base['Title'], fontName='Helvetica-Bold', fontSize=16, leading=19,
                           alignment=rl['TA_LEFT'], textColor=INK, spaceAfter=2),
            'sub': style('sub', parent=base['Normal'], fontSize=8.5, leading=11, textColor=MUTED,
                         spaceAfter=8),
            'h2': style('h2', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=12.5, leading=15,
                        textColor=INK, spaceBefore=10, spaceAfter=4),
            'h3': style('h3', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=10, leading=12.5,
                        textColor=INK, spaceBefore=8, spaceAfter=2),
            'h4': style('h4', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=8, leading=10,
                        textColor=MUTED, spaceBefore=7, spaceAfter=1),
            'body': style('body', parent=base['Normal'], fontSize=8.5, leading=11, textColor=INK,
                          spaceAfter=3),
            'cell': style('cell', parent=base['Normal'], fontSize=7.2, leading=8.8, textColor=INK),
            'cellb': style('cellb', parent=base['Normal'], fontName='Helvetica-Bold', fontSize=7.2,
                           leading=8.8, textColor=INK),
            'mono': style('mono', parent=base['Normal'], fontName='Courier', fontSize=6.8, leading=8.4,
                          textColor=MUTED),
        }
        self.width = rl['A4'][0] - 24 * mm  # The frame of the document, less its margins

    # -- pieces -------------------------------------------------------------------------------------
    def para(self, text, style='body'):
        self.story.append(self.rl['Paragraph'](escape(text), self.s[style]))

    def rich(self, text, style='body'):
        """A paragraph whose markup is already written (a <b> or a <font color>)."""
        self.story.append(self.rl['Paragraph'](text, self.s[style]))

    def space(self, points=4):
        self.story.append(self.rl['Spacer'](1, points))

    def facts(self, pairs, widths=(46, 54)):
        """A two column table of (term, value): the PDF form of the HTML definition lists."""
        rows = [[self.rl['Paragraph'](escape(term), self.s['cell']),
                 self.rl['Paragraph'](escape(value), self.s['cell'])]
                for term, value in pairs if value not in (None, '', ())]
        if not rows:
            return
        total = self.width
        table = self.rl['Table'](rows, colWidths=[total * widths[0] / 100.0, total * widths[1] / 100.0])
        table.setStyle(self.rl['TableStyle']([
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TEXTCOLOR', (0, 0), (0, -1), MUTED),
            ('LINEBELOW', (0, 0), (-1, -2), 0.25, LINE),
            ('TOPPADDING', (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('LEFTPADDING', (0, 0), (-1, -1), 0),
        ]))
        self.story.append(table)
        self.space(3)

    def table(self, headers, rows, widths=None, numeric=(), highlight=None, keep=True):
        """
        A table of values, the header row shaded and repeated on every page it continues onto.

        :param rows: list of rows of values
        :param widths: column widths as percentages of the frame
        :param highlight: {row index: colour} for the rows a reader should land on
        """
        rl = self.rl
        data = [[rl['Paragraph'](escape(h), self.s['cellb']) for h in headers]]
        for row in rows:
            data.append([rl['Paragraph'](escape(v), self.s['cell']) for v in row])
        total = self.width
        columns = [total * w / 100.0 for w in widths] if widths else [total / len(headers)] * len(headers)
        table = rl['Table'](data, colWidths=columns, repeatRows=1)
        commands = [
            ('BACKGROUND', (0, 0), (-1, 0), PANEL),
            ('GRID', (0, 0), (-1, -1), 0.25, LINE),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 2.5),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2.5),
        ]
        for column in numeric:
            commands.append(('ALIGN', (column, 1), (column, -1), 'RIGHT'))
        for index, colour in (highlight or {}).items():
            commands.append(('BACKGROUND', (0, index + 1), (-1, index + 1), colour))
        table.setStyle(rl['TableStyle'](commands))
        self.story.append(rl['KeepTogether'](table) if keep else table)
        self.space(3)

    def sign_off(self):
        """A box for the reviewer: name, date and signature, with room to write."""
        rl = self.rl
        third = self.width / 3.0
        rows = [['', '', ''],
                [rl['Paragraph']('Name', self.s['cell']), rl['Paragraph']('Date', self.s['cell']),
                 rl['Paragraph']('Signature', self.s['cell'])]]
        table = rl['Table'](rows, colWidths=[third] * 3, rowHeights=[13 * rl['mm'], None])
        table.setStyle(rl['TableStyle']([
            ('LINEBELOW', (0, 0), (-1, 0), 0.5, INK),
            ('TEXTCOLOR', (0, 1), (-1, 1), MUTED),
            ('VALIGN', (0, 0), (-1, -1), 'BOTTOM'),
            ('LEFTPADDING', (0, 0), (-1, -1), 0),
        ]))
        self.para('Reviewed by', 'h4')
        self.story.append(rl['KeepTogether'](table))

    # -- sections -----------------------------------------------------------------------------------
    def summary(self, results, run):
        counts = counted(results)
        self.para(TITLE, 'title')
        self.para('{} sample(s): {} ok, {} with a warning, {} failed. BrucellaID {}, run on {}.'.format(
            len(results), counts[OK], counts[WARNING], counts[FAILED], __version__,
            run.get('started', '')), 'sub')
        self.para('Summary', 'h2')
        self.facts([('Run started', run.get('started', '')), ('Run finished', run.get('finished', '')),
                    ('Operator', run.get('operator', '')), ('Command', run.get('command', ''))])
        for (species, st), members in groups(results):
            label = '<b>{}</b> &mdash; ST {} &mdash; {} sample(s)'.format(
                escape(species), escape(st), len(members)) if species else \
                '<font color="{}"><b>Failed</b></font> &mdash; {} sample(s)'.format(FAILED and
                                                                                    STATUS_COLOURS[FAILED],
                                                                                    len(members))
            self.rich(label, 'h4')
            rows = []
            for result in members:
                rows.append([result.sample, result.species or UNKNOWN, result.biovar,
                             result.call.source if result.call else '',
                             result.call.verdict if result.call and result.call.verdict else '',
                             result.sequence_type, result.mlst.profile if result.mlst else '',
                             status(result)])
            self.table(['Sample', 'Species', 'Biovar', 'From', 'Decision', 'ST', 'Allele profile', 'Status'],
                       rows, widths=(15, 16, 8, 10, 12, 5, 21, 13),
                       highlight={i: ACCENT_SOFT for i, r in enumerate(members) if status(r) != OK})
        self.sign_off()

    def warnings(self, results):
        self.story.append(self.rl['PageBreak']())
        self.para('Notes, warnings and errors', 'h2')
        rows = [[sample, state, text] for sample, state, text in warnings(results)]
        if not rows:
            self.para('None: every sample was identified and typed without a note.')
            return
        self.table(['Sample', 'Status', 'What it says'], rows, widths=(16, 10, 74), keep=False)

    def sample(self, result, references, inputs):
        self.story.append(self.rl['PageBreak']())
        state = status(result)
        self.rich('{} &nbsp;<font size="8" color="{}">{}</font> &nbsp;<font size="8" color="{}">{} '
                  'ST {}</font>'.format(escape(result.sample), STATUS_COLOURS[state], state, ACCENT,
                                        escape(result.species or UNKNOWN), escape(result.sequence_type)),
                  'h2')
        if result.error:
            self.rich('<b>This sample failed.</b> {}'.format(escape(result.error)))

        self.para('Identification: the MASH distances', 'h4')
        if result.identification is None:
            self.para('MASH did not run for this sample.')
        else:
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
            ])
            rows, highlight = [], {}
            for index, (hit, reference, close) in enumerate(close_hits(result, references)):
                if hit is best or close:
                    highlight[index] = ACCENT_SOFT if hit is best else PANEL
                rows.append([hit.reference, reference.species if reference else '(not in the table)',
                             reference.biovar if reference else '',
                             '{:.6f}'.format(hit.distance), '{:.3f}'.format(hit.identity),
                             '{}/{}'.format(hit.shared_hashes, hit.total_hashes)])
            self.table(['Reference genome', 'Species', 'Biovar', 'Distance', 'Identity %', 'Shared hashes'],
                       rows, widths=(31, 23, 8, 13, 12, 13), numeric=(3, 4, 5), highlight=highlight,
                       keep=False)

        self.para('The species call', 'h4')
        call = result.call
        if call is None:
            self.para('No species was called for this sample.')
        else:
            self.facts([
                ('Species reported', call.species),
                ('Biovar', '{} (of {}, the closest reference of that species)'.format(
                    call.biovar, call.reference.strain) if call.biovar and call.reference else call.biovar),
                ('Decided by', call.source),
                ('Species within reach', ', '.join(call.candidates)),
                ('Sequence type', call.st),
                ('What that ST is', call.entry.evidence if call.entry else ''),
                ('Why the ST was not used', call.st_reason),
            ])
            if call.verdict:
                self.rich('<b>{}</b>: {}'.format(escape(call.verdict), escape(VERDICTS[call.verdict])))

        self.para('MLST: the nine allele calls', 'h4')
        if result.mlst is None:
            self.para('The MLST was not run for this sample.')
        else:
            mlst = result.mlst
            self.facts([('Sequence type', mlst.st), ('Allele profile', mlst.profile),
                        ('Every locus called', 'yes' if mlst.complete else 'no'),
                        ('Reads in the input', '{:,}'.format(mlst.total_reads)
                         if mlst.total_reads is not None else ''),
                        ('Reads matching an allele', '{:,}'.format(mlst.assigned_reads)),
                        ('About the ST', mlst.note)])
            rows = []
            for locus_call in mlst.calls:
                ahead, ratio = margin(locus_call)
                rows.append([locus_call.locus, locus_call.allele or UNKNOWN, locus_call.reads or '',
                             '{}_{}'.format(locus_call.locus, locus_call.runner_up)
                             if locus_call.runner_up else '', locus_call.runner_up_reads or '',
                             '' if ahead is None else ahead,
                             '' if ratio is None else '{:.2f}'.format(ratio), locus_call.note])
            self.table(['Locus', 'Allele', 'Reads', 'Runner-up', 'Its reads', 'Ahead', 'Ratio', 'Note'],
                       rows, widths=(11, 8, 8, 14, 9, 8, 8, 34), numeric=(2, 4, 5, 6), keep=False)

        self.para('Input files', 'h4')
        files = inputs.get(result.sample, [])
        if not files:
            self.para('The input files could not be read when the report was written.')
        else:
            self.table(['File', 'Size', 'Modified', 'MD5'],
                       [[f.name, f.size_text, f.modified, f.md5 or 'not computed (--no-md5)'] for f in files],
                       widths=(32, 10, 26, 32), numeric=(1,), keep=False)
            for item in files:
                self.rich(escape(item.path) + (' -> {}'.format(escape(item.target))
                                                               if item.target else ''), 'mono')

        if result.notes:
            self.para('Notes', 'h4')
            for note in result.notes:
                if note:
                    self.rich('- {}'.format(escape(note)), 'cell')

    def methods(self):
        self.story.append(self.rl['PageBreak']())
        self.para('Definitions and methods', 'h2')
        for heading, text in METHODS:
            self.para(heading, 'h3')
            self.para(text)

    def run_information(self, run, data):
        self.story.append(self.rl['PageBreak']())
        self.para('Run information', 'h2')
        self.facts([('Operator', run.get('operator', '')), ('User account', run.get('user', '')),
                    ('Computer', run.get('host', '')), ('Operating system', run.get('system', '')),
                    ('Working directory', run.get('working_directory', '')),
                    ('Started', run.get('started', '')), ('Finished', run.get('finished', '')),
                    ('Command', run.get('command', ''))])
        self.para('Software', 'h3')
        self.table(['Program', 'Version'], [[name, value] for name, value in run.get('software', {}).items()],
                   widths=(28, 72))
        self.para('Parameters', 'h3')
        self.table(['Parameter', 'Value'],
                   [[name, '' if value is None else value] for name, value in run.get('parameters', {}).items()],
                   widths=(34, 66), keep=False)
        self.para('Reference data', 'h3')
        if data:
            self.table(['Data', 'What it holds', 'Source', 'MD5'],
                       [[item['name'], item['detail'],
                         'the package\'s' if item['packaged'] else 'given with an option',
                         item['md5'] or 'not computed'] for item in data],
                       widths=(18, 36, 18, 28), keep=False)
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
                self.setFont('Helvetica', 6.5)
                self.setFillColor(MUTED)
                self.drawString(margin, baseline, left)
                self.drawRightString(self._pagesize[0] - margin, baseline,
                                     'Page {} of {}'.format(self._pageNumber, total))
                self.restoreState()
                super().showPage()
            super().save()

    return NumberedCanvas


def write(results, references, run, path, inputs=None, data=None):
    """
    Write the PDF report and return its path.

    :raises BrucellaIDError: when reportlab is not installed
    """
    rl = _modules()
    inputs = inputs or {}
    mm = rl['mm']
    document = rl['SimpleDocTemplate'](
        str(path), pagesize=rl['A4'], leftMargin=12 * mm, rightMargin=12 * mm, topMargin=14 * mm,
        bottomMargin=16 * mm, title=TITLE, author=run.get('operator', ''), subject='BrucellaID run',
        creator='BrucellaID {}'.format(__version__))
    report = Report(rl)
    report.summary(results, run)
    report.warnings(results)
    for _, members in groups(results):
        for result in members:
            report.sample(result, references, inputs)
    report.methods()
    report.run_information(run, data)
    document.build(report.story, canvasmaker=numbered_canvas(run, 12 * mm, 10 * mm))
    return path
