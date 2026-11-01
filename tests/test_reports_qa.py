"""The quality assurance reports: the self-contained HTML, the PDF and the MultiQC file."""

import json
import sys
from html.parser import HTMLParser

import pytest

from brucellaid import htmlreport, qa, report
from brucellaid.errors import BrucellaIDError
from brucellaid.mash import Identification, parse_dist
from brucellaid.mlst import UNKNOWN, AlleleCall, MlstResult
from brucellaid.pipeline import SampleResult, species_notes
from brucellaid.references import COLUMNS, parse_references
from brucellaid.species import call_species, load_st_species

REFERENCES = parse_references([
    '\t'.join(COLUMNS) + '\n',
    'NC_006932-NC_006933.fasta\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1;NC_006933.1\n',
    'NC_010103-NC_010104.fasta\tcanis\tBrucella canis\t\tATCC 23365\tNC_010103.1\n',
    'B-REF-BS4-40.fasta\tsuis4\tBrucella suis\t4\t40\t\n',
])
CANIS_HITS = ('NC_010103-NC_010104.fasta\tq\t0.000239\t0\t990/1000\n'
              'B-REF-BS4-40.fasta\tq\t0.000420\t0\t980/1000\n'
              'NC_006932-NC_006933.fasta\tq\t0.0029\t0\t888/1000\n')
RUN = {'command': 'brucellaid -i reads -o out', 'started': '2026-10-31T09:15:00-04:00',
       'finished': '2026-10-31T09:15:37-04:00', 'operator': 'A Reviewer', 'user': 'lab', 'host': 'bench1',
       'system': 'Linux-6.8', 'working_directory': '/data/run42',
       'software': {'brucellaid': '1.0.0', 'MASH': '2.3'},
       'parameters': {'Threads': 8, 'MLST': True, 'Minimum reads per allele': None}}
DATA = [{'name': 'MASH sketch', 'path': '/pkg/references.msh', 'detail': '16 reference genome(s)',
         'md5': 'abc123', 'packaged': True},
        {'name': 'ST to species', 'path': '/mine/st.tsv', 'detail': '118 ST(s)', 'md5': '', 'packaged': False}]


def good(name='ab1', hits='NC_006932-NC_006933.fasta\tq\t0.0012345\t0\t950/1000\n', st='1', files=()):
    found = Identification(parse_dist(hits))
    reference, notes = species_notes(found, REFERENCES)
    # st=None stands for --no-mlst: there is no MlstResult at all, not one with an empty ST
    mlst = MlstResult(st=st, calls=[AlleleCall('gap', '2', 40, '3', 12),
                                    AlleleCall('aroA', None, 0, '1', 2, note='no read matched')],
                      assigned_reads=900, total_reads=1000) if st else None
    call = call_species(found, REFERENCES, mlst=mlst, table=load_st_species())
    return SampleResult(name, list(files), identification=found, reference=reference, mlst=mlst, call=call,
                        notes=list(notes))


def failed(name='bad'):
    return SampleResult(name, ['bad_R1.fastq.gz'], error='Seal failed (exit code 1):\tno\noutput')


def bare(name='nothing'):
    """A sample MASH never ran on: the report must not fall over on it."""
    return SampleResult(name, [])


class Structure(HTMLParser):
    """Checks that every tag is closed, that no id repeats and that every anchor points somewhere."""

    VOID = {'meta', 'br', 'hr', 'img', 'input', 'link', 'source', 'col'}

    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.stack, self.problems, self.ids, self.refs = [], [], [], []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if attributes.get('id'):
            self.ids.append(attributes['id'])
        if attributes.get('href', '').startswith('#'):
            self.refs.append(attributes['href'][1:])
        if tag not in self.VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in self.VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.problems.append('</{}> does not close <{}>'.format(tag, self.stack[-1:] or ['nothing']))
        if self.stack:
            self.stack.pop()


def parse(markup):
    checker = Structure()
    checker.feed(markup)
    assert checker.problems == []
    assert checker.stack == []
    assert len(set(checker.ids)) == len(checker.ids)
    assert set(checker.refs) <= set(checker.ids)
    return checker


# -- the HTML report ------------------------------------------------------------------------------------

def test_html_report_is_well_formed_and_holds_every_section(tmp_path):
    results = [good(), good('canis', hits=CANIS_HITS, st='20'), failed(), bare()]
    inputs = qa.describe_inputs(results)
    markup = htmlreport.build(results, REFERENCES, RUN, inputs=inputs, data=DATA)
    parse(markup)
    for section in ('id="summary"', 'id="warnings"', 'id="samples"', 'id="methods"', 'id="run"'):
        assert section in markup
    assert '<!doctype html>' in markup
    assert 'Reviewed by' in markup and 'Signature' in markup


def test_html_report_carries_the_logo():
    markup = htmlreport.build([good()], REFERENCES, RUN)
    assert '<div class="logo"><svg' in markup
    assert 'var(--logo-accent' in markup  # Themed, so one drawing serves light and dark
    assert '<h1>Quality assurance report</h1>' in markup  # The wordmark is in the logo
    assert '<title>BrucellaID quality assurance report</title>' in markup
    parse(markup)


def test_html_report_without_the_logo(monkeypatch):
    """The report is still written, with the name as a heading instead of the drawing."""
    htmlreport.logo.cache_clear()
    monkeypatch.setattr('brucellaid.htmlreport.data_path',
                        lambda name: (_ for _ in ()).throw(BrucellaIDError('no logo')))
    markup = htmlreport.build([good()], REFERENCES, RUN)
    htmlreport.logo.cache_clear()
    assert '<div class="logo">' not in markup
    assert '<h1>BrucellaID quality assurance report</h1>' in markup


def test_html_report_needs_no_network():
    """The point of the report: it renders the same in ten years, on a machine with no network."""
    markup = htmlreport.build([good()], REFERENCES, RUN)
    assert '<script' not in markup
    assert 'src=' not in markup
    assert '<link' not in markup
    # The only URLs are the ones printed inside citations, not resources the page fetches
    assert 'href="http' not in markup


def test_html_report_holds_every_distance_not_only_the_winner():
    markup = htmlreport.build([good('canis', hits=CANIS_HITS, st='20')], REFERENCES, RUN)
    for reference in ('NC_010103-NC_010104.fasta', 'B-REF-BS4-40.fasta', 'NC_006932-NC_006933.fasta'):
        assert reference in markup
    assert 'class="best"' in markup  # The closest reference
    assert 'class="close"' in markup  # The one near enough that the distance does not separate them
    assert '0.000420' in markup


def test_html_report_explains_how_the_species_was_decided():
    markup = htmlreport.build([good('canis', hits=CANIS_HITS, st='20')], REFERENCES, RUN)
    assert 'ST 20 is Brucella canis in 70 of 70 isolate(s)' in markup
    assert 'arbitrated' in markup
    assert htmlreport.VERDICTS['arbitrated'].split('.')[0] in markup


def test_html_report_says_why_an_st_was_not_used():
    """ST 17 is B. suis in 18 of 20 isolates: below the bar, and the report has to say so."""
    markup = htmlreport.build([good(st='17')], REFERENCES, RUN)
    assert 'Why the ST was not used' in markup
    assert 'below the 95% agreement required' in markup


def test_html_report_shows_the_allele_margins():
    markup = htmlreport.build([good()], REFERENCES, RUN)
    assert 'Ahead by' in markup and 'Ratio' in markup
    assert '3.33' in markup  # 40 reads against 12 for the runner-up
    assert 'no read matched' in markup


def test_html_report_of_a_failed_sample_and_one_with_no_identification():
    markup = htmlreport.build([failed(), bare()], REFERENCES, RUN)
    parse(markup)
    assert 'This sample failed.' in markup
    assert 'MASH did not run for this sample.' in markup
    assert 'No species was called for this sample.' in markup
    assert 'The MLST was not run for this sample.' in markup
    assert 'The input files could not be read when the report was written.' in markup
    # The error is shown twice: once in the warnings table, once on the sample itself. HTML collapses the
    # tab and the line break Seal put in it, so it reads as one line either way.
    assert markup.count('Seal failed (exit code 1)') == 2


def test_html_report_skips_an_empty_note():
    result = good()
    result.notes = ['', 'nearly as close: something']
    markup = htmlreport.build([result], REFERENCES, RUN)
    assert '<li class="wrap"></li>' not in markup
    assert '<li class="wrap">nearly as close: something</li>' in markup


def test_html_report_without_warnings_says_so():
    markup = htmlreport.build([good()], REFERENCES, RUN)
    assert 'None: every sample was identified and typed without a note.' in markup


def test_html_report_without_reference_data_says_so():
    markup = htmlreport.build([good()], REFERENCES, RUN, data=[])
    assert 'The reference data could not be read when the report was written.' in markup


def test_html_report_lists_the_reference_data_with_its_checksums():
    markup = htmlreport.build([good()], REFERENCES, RUN, data=DATA)
    assert 'abc123' in markup
    assert "the package&#x27;s" in markup
    assert 'given with an option' in markup
    assert 'not computed (--no-md5)' in markup


def test_html_report_escapes_what_a_sample_name_may_hold():
    result = good('<script>alert(1)</script>')
    markup = htmlreport.build([result], REFERENCES, RUN)
    assert '<script>' not in markup
    assert '&lt;script&gt;' in markup


def test_html_anchor_of_an_awkward_sample_name():
    assert htmlreport.anchor('A/b c.1') == 'sample-a-b-c-1'
    assert htmlreport.anchor('***') == 'sample-'


def test_html_report_input_files_with_and_without_md5(tmp_path):
    reads = tmp_path / 'ab1_R1.fastq.gz'
    reads.write_bytes(b'reads')
    results = [good(files=[str(reads)])]
    with_md5 = htmlreport.build(results, REFERENCES, RUN, inputs=qa.describe_inputs(results))
    assert qa.file_md5(reads) in with_md5
    without = htmlreport.build(results, REFERENCES, RUN, inputs=qa.describe_inputs(results, md5=False))
    assert 'not computed (--no-md5)' in without


def test_html_report_is_written_to_a_file(tmp_path):
    path = htmlreport.write([good()], REFERENCES, RUN, tmp_path / 'report.html')
    assert path.read_text(encoding='utf-8').startswith('<!doctype html>')


def test_html_table_helper_aligns_numbers_and_takes_a_caption():
    markup = htmlreport.table(['A', 'N'], [['x', 1], (['y', 2], 'best')], caption='Why', numeric=(1,))
    assert '<caption>Why</caption>' in markup
    assert '<th class="num">N</th>' in markup
    assert '<td class="num">1</td>' in markup
    assert '<tr class="best">' in markup


def test_html_tag_helper():
    assert htmlreport.tag('ok', 'ok') == '<span class="tag ok">ok</span>'
    assert htmlreport.tag('plain') == '<span class="tag">plain</span>'


# -- MultiQC --------------------------------------------------------------------------------------------

def test_multiqc_file(tmp_path):
    path = tmp_path / 'brucellaid_mqc.json'
    report.write_multiqc([good(), good('canis', hits=CANIS_HITS, st='20'), failed()], path)
    content = json.loads(path.read_text())
    assert content['id'] == 'brucellaid'
    assert content['plot_type'] == 'table'
    assert set(content['data']) == {'ab1', 'canis', 'bad'}
    assert content['data']['ab1'] == {'Species': 'Brucella abortus', 'Biovar': '1',
                                      'SpeciesSource': 'MASH+ST (agree)', 'STSpecies': 'Brucella abortus',
                                      'SequenceType': '1', 'SharedHashes': '950/1000',
                                      'ClosestReference': 'NC_006932-NC_006933.fasta', 'Status': 'ok'}
    # Nothing is left blank: MultiQC renders an empty cell as a hole in the table
    assert all(value for row in content['data'].values() for value in row.values())
    assert content['data']['bad']['Status'] == 'failed'
    assert content['data']['canis']['SpeciesSource'] == 'MASH+ST (arbitrated)'
    assert set(content['headers']) == set(report.MQC_COLUMNS)


def test_multiqc_keeps_sequence_types_as_text(tmp_path):
    """JSON, not TSV, so that MultiQC does not average a column of sequence types."""
    path = tmp_path / 'x_mqc.json'
    report.write_multiqc([good()], path)
    assert '"SequenceType": "1"' in path.read_text()


# -- the PDF report -------------------------------------------------------------------------------------

def test_pdf_report_says_how_to_install_reportlab(monkeypatch, tmp_path):
    """Without reportlab, --pdf has to say what to install, not raise ImportError at the operator."""
    from brucellaid import pdfreport
    # Every reportlab module has to go: with a submodule left in sys.modules, the import succeeds without
    # ever looking at the package this test hides.
    for name in [name for name in sys.modules if name == 'reportlab' or name.startswith('reportlab.')]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.setitem(sys.modules, 'reportlab', None)
    with pytest.raises(BrucellaIDError, match='pip install brucellaid\\[pdf\\]'):
        pdfreport.write([good()], REFERENCES, RUN, tmp_path / 'r.pdf')


def test_pdf_report_carries_the_logo(tmp_path):
    """The wordmark is the logo image, so the page must really hold an image."""
    pytest.importorskip('reportlab')
    from brucellaid import pdfreport
    path = pdfreport.write([good()], REFERENCES, RUN, tmp_path / 'r.pdf')
    reader, _ = read_pdf(path)
    images = list(reader.pages[0].images)
    assert len(images) == 1
    assert images[0].image.width > 500  # The logo, not a stray pixel


def test_pdf_report_without_the_logo(monkeypatch, tmp_path):
    """A package whose data folder lost the image still writes the report."""
    pytest.importorskip('reportlab')
    from brucellaid import pdfreport
    monkeypatch.setattr('brucellaid.pdfreport.data_path',
                        lambda name: (_ for _ in ()).throw(BrucellaIDError('no logo')))
    path = pdfreport.write([good()], REFERENCES, RUN, tmp_path / 'r.pdf')
    reader, text = read_pdf(path)
    assert list(reader.pages[0].images) == []
    assert 'Quality assurance report' in text


def test_pdf_escape():
    from brucellaid import pdfreport
    assert pdfreport.escape('a & b <i> c') == 'a &amp; b &lt;i&gt; c'
    assert pdfreport.escape(None) == ''


def read_pdf(path):
    """
    (reader, text) of a PDF, with the whitespace of the text normalised.

    Text extraction breaks lines wherever the PDF happens to start a new text object, so a sentence of the
    report comes back with newlines in the middle of it. Collapsing the whitespace is what makes it possible
    to assert that a sentence is there at all.
    """
    pypdf = pytest.importorskip('pypdf', reason='pypdf reads the report back')
    reader = pypdf.PdfReader(str(path))
    text = ' '.join((page.extract_text() or '') for page in reader.pages)
    return reader, ' '.join(text.split())


def test_pdf_report(tmp_path):
    pytest.importorskip('reportlab')
    reads = tmp_path / 'ab1_R1.fastq.gz'
    reads.write_bytes(b'reads')
    results = [good(files=[str(reads)]), good('canis', hits=CANIS_HITS, st='20'), failed(), bare()]
    inputs = qa.describe_inputs(results)
    from brucellaid import pdfreport
    path = pdfreport.write(results, REFERENCES, RUN, tmp_path / 'report.pdf', inputs=inputs, data=DATA)
    reader, text = read_pdf(path)
    assert len(reader.pages) > 4
    assert reader.metadata['/Title'] == pdfreport.TITLE
    assert reader.metadata['/Author'] == 'A Reviewer'
    for expected in ('Quality assurance report', 'Summary', 'Reviewed by', 'Signature',
                     'Notes, warnings and errors', 'Definitions and methods', 'Run information',
                     'ST 20 is Brucella canis', 'B-REF-BS4-40.fasta', 'Seal failed',
                     'MASH did not run for this sample.', 'abc123', 'doi:10.1186/1471-2180-7-34'):
        assert expected in text, expected
    # Every page says which report it belongs to and whether a page is missing
    assert 'Page 1 of {}'.format(len(reader.pages)) in text


def test_pdf_report_without_warnings_or_reference_data(tmp_path):
    pytest.importorskip('reportlab')
    from brucellaid import pdfreport
    path = pdfreport.write([good(st=None)], REFERENCES, RUN, tmp_path / 'r.pdf', data=[])
    _, text = read_pdf(path)
    assert 'None: every sample was identified and typed without a note.' in text
    assert 'The reference data could not be read when the report was written.' in text
    assert 'The MLST was not run for this sample.' in text


def test_pdf_report_without_any_run_information(tmp_path):
    """A report built from an empty run dictionary: every fact table is empty and is left out."""
    pytest.importorskip('reportlab')
    from brucellaid import pdfreport
    result = good()
    result.notes = ['', 'nearly as close: something']  # An empty note is skipped, not printed as a bullet
    path = pdfreport.write([result], REFERENCES, {}, tmp_path / 'r.pdf')
    _, text = read_pdf(path)
    assert 'Quality assurance report' in text
    assert '- nearly as close: something' in text


def test_pdf_report_of_a_symbolic_link_and_no_md5(tmp_path):
    pytest.importorskip('reportlab')
    from brucellaid import pdfreport
    real = tmp_path / 'real_R1.fastq.gz'
    real.write_bytes(b'reads')
    link = tmp_path / 'ab1_R1.fastq.gz'
    link.symlink_to(real)
    results = [good(files=[str(link)])]
    path = pdfreport.write(results, REFERENCES, RUN, tmp_path / 'r.pdf',
                           inputs=qa.describe_inputs(results, md5=False))
    _, text = read_pdf(path)
    assert 'not computed (--no-md5)' in text
    assert 'real_R1.fastq.gz' in text


# -- write_reports --------------------------------------------------------------------------------------

def test_write_reports_writes_the_qa_reports(tmp_path):
    from datetime import datetime
    results = [good()]
    rows, paths = report.write_reports(results, REFERENCES, tmp_path, datetime(2026, 10, 31), run=RUN)
    names = [path.name for path in paths]
    assert htmlreport.REPORT_FILE in names
    assert report.MQC_FILE in names
    assert report.PDF_FILE not in names  # --pdf was not asked for
    assert rows[0][0] == 'Sample'


def test_write_reports_can_skip_the_html_report(tmp_path):
    from datetime import datetime
    _, paths = report.write_reports([good()], REFERENCES, tmp_path, datetime(2026, 10, 31), run=RUN,
                                    html=False, md5=False)
    assert htmlreport.REPORT_FILE not in [path.name for path in paths]
    assert not (tmp_path / htmlreport.REPORT_FILE).exists()


def test_write_reports_writes_the_pdf_when_asked(tmp_path):
    pytest.importorskip('reportlab')
    from datetime import datetime
    _, paths = report.write_reports([good()], REFERENCES, tmp_path, datetime(2026, 10, 31), run=RUN, pdf=True)
    assert report.PDF_FILE in [path.name for path in paths]
    assert (tmp_path / report.PDF_FILE).stat().st_size > 1000


def test_json_report_holds_the_input_files_and_the_status(tmp_path):
    reads = tmp_path / 'ab1_R1.fastq.gz'
    reads.write_bytes(b'reads')
    results = [good(files=[str(reads)])]
    content = report.json_report(results, REFERENCES, run=RUN, inputs=qa.describe_inputs(results))
    sample = content['samples'][0]
    assert sample['status'] == 'ok'
    assert sample['input_files'][0]['md5'] == qa.file_md5(reads)
    assert sample['input_files'][0]['symlink_to'] is None
    assert sample['species_call']['st_not_used'] is None


def test_json_report_says_why_an_st_was_not_used():
    content = report.json_report([good(st='17')], REFERENCES, run=RUN)
    assert 'below the 95% agreement required' in content['samples'][0]['species_call']['st_not_used']
    assert content['samples'][0]['input_files'] == []


def test_run_info_holds_the_environment_and_the_reference_data():
    from datetime import datetime
    info = report.run_info('brucellaid -r1 a.fastq -o out', {'Threads': 4}, datetime(2026, 10, 31),
                           operator='A Reviewer', data=DATA)
    assert info['operator'] == 'A Reviewer'
    assert info['user'] == qa.user_name()
    assert info['reference_data'] == DATA
    assert info['software']['brucellaid']
    assert UNKNOWN not in info  # Nothing odd leaked into the dictionary
