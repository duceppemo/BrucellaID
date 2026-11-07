#!/usr/bin/env python3
"""
Make the example quality assurance reports of the wiki, with the machine they ran on left out.

    python3 scripts/make_example_report.py READS_FOLDER            # writes docs/example_report.html and .pdf

A quality assurance report names the computer it was produced on, the account that produced it and the full
path of every file it touched. That is the point of it, and it is exactly what must not be committed to a
public repository.

So the example is not a report with those values painted over afterwards: the run is real, and the reports
are rendered from it with those values replaced before anything is written. That is the only way the PDF
can be redacted at all -- there is no sensible way to edit text out of a finished PDF -- and it means the
HTML and the PDF are redacted identically, by construction rather than by two sets of rules.

What is left out: the user account, the host name, the operating system and the working directory, each
replaced by [removed], and the folder of every input and reference data file, leaving the file name. What
stays: every measurement, every checksum, the parameters, the versions, the timestamps and the operator
name, which is a byline rather than a machine.

The reads are the four public read sets of the tutorial. Fetch them first:
https://github.com/duceppemo/BrucellaID/wiki/Tutorial
"""
import argparse
import dataclasses
import pathlib
import shlex
import sys
from datetime import datetime

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from brucellaid import htmlreport, pdfreport, qa  # noqa: E402
from brucellaid.mlst import PROFILES_FILE, alleles_fasta, load_profiles  # noqa: E402
from brucellaid.pipeline import analyse_samples  # noqa: E402
from brucellaid.references import (  # noqa: E402
    EXTENDED_REFERENCES_FILE,
    data_path,
    load_references,
    reference_sketch,
)
from brucellaid.report import run_info  # noqa: E402
from brucellaid.samples import find_samples  # noqa: E402
from brucellaid.species import ST_SPECIES_FILE, load_st_species  # noqa: E402

REMOVED = '[removed]'
HIDDEN = ('user', 'host', 'system', 'working_directory')
NOTE = ('This is an example report: a real run of the four public read sets of the BrucellaID tutorial, '
        'produced by scripts/make_example_report.py, which leaves out the computer it ran on. The user '
        'account, the host name, the operating system and the working directory read [removed], and the '
        'input and reference files are named without the folder they were in. Everything else is as '
        'BrucellaID wrote it. A report of your own names them all, which is what makes it worth signing.')


def redact_run(run):
    """The run information with the machine taken out, and the reference data files named without a folder."""
    run = dict(run)
    for key in HIDDEN:
        if run.get(key):
            run[key] = REMOVED
    run['reference_data'] = [{**item, 'path': '.../brucellaid/data/{}'.format(pathlib.PurePath(
        item['path']).name)} for item in run.get('reference_data', [])]
    return run


def redact_inputs(inputs):
    """The input files named without the folder they were in, and with no symbolic link target."""
    return {sample: [dataclasses.replace(facts, path=pathlib.PurePath(facts.path).name, target='')
                     for facts in files]
            for sample, files in inputs.items()}


def secrets(run, inputs, reads):
    """
    Every string that must not survive into a report, taken from the run before it was redacted.

    Checked against the real values rather than against the redacted ones -- the first version of this
    script compared the output with the values it had already replaced, which can only ever pass.
    """
    values = {run.get(key) for key in HIDDEN}
    values.add(str(pathlib.Path.cwd()))
    values.add(str(pathlib.Path(reads).resolve()))
    for files in inputs.values():
        for facts in files:
            values.add(str(pathlib.PurePath(facts.path).parent))
    for item in run.get('reference_data', []):
        values.add(str(pathlib.PurePath(item['path']).parent))
    return {value for value in values if value and value not in (REMOVED, '.', '/')}


def leaks(text, values):
    """The values that are still in a report, so that the script can refuse to write it."""
    return sorted(value for value in values if value in text)


def pdf_text(path):
    """The text of a PDF, with its whitespace collapsed, for checking what is in it."""
    try:
        from pypdf import PdfReader
    except ImportError:
        raise SystemExit('Checking the PDF needs pypdf: pip install pypdf. The report is not written '
                         'unless it can be checked.') from None
    reader = PdfReader(str(path))
    return ' '.join(' '.join((page.extract_text() or '') for page in reader.pages).split())


def reference_data(sketch, references, loci, st_table):
    """The reference data of the run, as the command line builds it."""
    alleles = alleles_fasta(None)
    return qa.data_facts([
        ('MASH sketch', sketch, '{} reference genome(s), k = 21, 1000 hashes'.format(len(references))),
        ('Reference table', data_path(EXTENDED_REFERENCES_FILE),
         '{} genome(s) with their species, biovar, strain and accessions'.format(len(references))),
        ('MLST alleles', alleles, '{} allele(s) of the {} loci'.format(qa.fasta_count(alleles), len(loci))),
        ('MLST profiles', data_path(PROFILES_FILE),
         '{} sequence type(s)'.format(len(load_profiles(None)[1]))),
        ('ST to species', data_path(ST_SPECIES_FILE),
         '{} sequence type(s) with a species'.format(len(st_table))),
    ])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('reads', help='folder holding the fastq files of the tutorial')
    parser.add_argument('-o', '--out', default='docs', metavar='FOLDER',
                        help='where to write example_report.html and example_report.pdf (default: %(default)s)')
    parser.add_argument('--operator', default='M-O Duceppe', help='the name on the report')
    parser.add_argument('-t', '--threads', type=int, default=8)
    parser.add_argument('-p', '--parallel', type=int, default=2)
    args = parser.parse_args(argv)

    samples = find_samples(args.reads)
    if not samples:
        raise SystemExit('No fastq file in {}.'.format(args.reads))
    print('{} sample(s): {}'.format(len(samples), ', '.join(s.name for s in samples)), file=sys.stderr)

    references = load_references()
    st_table = load_st_species()
    loci = load_profiles(None)[0]
    started = datetime.now().astimezone()
    results = analyse_samples(samples, jobs=args.parallel, threads=args.threads, references=references,
                              st_table=st_table)

    # The command is recorded with the reads folder named, not its path: it is a path like any other
    command = shlex.join(['python3', 'scripts/make_example_report.py',
                          pathlib.PurePath(args.reads.rstrip('/')).name])
    run = run_info(command, {'Threads': args.threads, 'Parallel samples': args.parallel, 'MLST': True,
                             'Input checksums': 'MD5'}, started, operator=args.operator,
                   data=reference_data(reference_sketch(None), references, loci, st_table))
    real_inputs = qa.describe_inputs(results)
    must_go = secrets(run, real_inputs, args.reads)
    run = redact_run(run)
    inputs = redact_inputs(real_inputs)

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    markup = htmlreport.build(results, references, run, inputs=inputs, data=run['reference_data'], note=NOTE)
    left = leaks(markup, must_go)
    if left:
        raise SystemExit('Still in the HTML report: {}'.format(', '.join(left)))
    (out / 'example_report.html').write_text(markup, encoding='utf-8')

    # The PDF is written to a temporary name and only kept once it has been read back and found clean
    draft = out / 'example_report.pdf.part'
    pdfreport.write(results, references, run, draft, inputs=inputs, data=run['reference_data'], note=NOTE)
    left = leaks(pdf_text(draft), must_go)
    if left:
        draft.unlink()
        raise SystemExit('Still in the PDF report: {}'.format(', '.join(left)))
    draft.replace(out / 'example_report.pdf')

    print('checked against {} value(s) that had to go'.format(len(must_go)), file=sys.stderr)
    for name in ('example_report.html', 'example_report.pdf'):
        print('{}: {:,} bytes'.format(out / name, (out / name).stat().st_size), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
