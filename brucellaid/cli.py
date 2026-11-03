"""Command line entry point."""

import logging
import os
import re
import shlex
import sys
import time
from argparse import ArgumentParser, ArgumentTypeError, RawDescriptionHelpFormatter
from datetime import datetime
from pathlib import Path

from . import __version__, qa
from .errors import BrucellaIDError
from .mash import DEFAULT_MIN_HASHES, KMER_COPIES_AUTO, check_mash, sketch_references
from .mlst import DEFAULT_MIN_READS, PROFILES_FILE, alleles_fasta, load_profiles
from .pipeline import analyse_samples
from .references import (
    EXTENDED_REFERENCES_FILE,
    EXTENDED_SKETCH_FILE,
    REFERENCES_FILE,
    compare_with_sketch,
    data_path,
    load_references,
    reference_sketch,
)
from .report import print_table, run_info, write_reports
from .samples import find_samples, single_sample
from .seal import DEFAULT_KMER_SIZE, DEFAULT_MISMATCHES, auto_memory, check_seal
from .species import DEFAULT_MIN_ISOLATES, DEFAULT_MIN_PURITY, ST_SPECIES_FILE, load_st_species

log = logging.getLogger('brucellaid')

# Expected errors are reported without a traceback
USER_ERRORS = (BrucellaIDError, OSError)

EPILOG = '''examples:
  brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
  brucellaid -r1 sample.fastq.gz -o results/
  brucellaid -i folder_of_fastq/ -o results/
  brucellaid -i folder_of_fastq/ -o results/ -p 4 -t 16 --xlsx

documentation: https://github.com/duceppemo/BrucellaID/wiki'''


def setup_logging(verbose=False):
    """Configure the package logger only, so importing brucellaid never alters the root logger."""
    for handler in list(log.handlers):
        log.removeHandler(handler)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-7s %(message)s', datefmt='%H:%M:%S'))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)


def positive_int(value):
    try:
        value = int(value)
    except ValueError:
        raise ArgumentTypeError('"{}" is not an integer'.format(value)) from None
    if value < 1:
        raise ArgumentTypeError('must be >= 1')
    return value


def java_memory(value):
    if not value[:-1].isdigit() or value[-1].lower() not in 'mg':
        raise ArgumentTypeError('"{}" is not a Java memory size such as 500m or 2g'.format(value))
    return value


def available_cpus():
    try:
        return len(os.sched_getaffinity(0))  # Respects cgroup/taskset limits on Linux (e.g. SLURM jobs)
    except AttributeError:
        return os.cpu_count() or 1


def format_elapsed(seconds):
    """A readable duration: "1h2m3s"."""
    minutes, seconds = divmod(int(round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
    return ''.join('{}{}'.format(value, name) for name, value in periods if value) or '0s'


def build_parser():
    max_cpu = available_cpus()
    parser = ArgumentParser(prog='brucellaid', formatter_class=RawDescriptionHelpFormatter, epilog=EPILOG,
                            description='Identify the Brucella species of sequencing reads by their MASH distance\n'
                                        'to 16 reference genomes, and type them with the nine-locus MLST scheme.\n'
                                        'Where two species are too close for MASH to choose, the sequence type\n'
                                        'chooses between them.')
    inputs = parser.add_argument_group('input (one sample with -r1, or a folder of samples with -i)')
    source = inputs.add_mutually_exclusive_group(required=True)
    source.add_argument('-r1', '--r1', metavar='FILE',
                        help='Reads: single-end, or R1 of paired-end. fastq, gzipped or not.')
    source.add_argument('-i', '--input', metavar='FOLDER',
                        help='Folder of fastq files, searched recursively. Each single-end fastq file or pair of '
                             'R1/R2 fastq files is a sample.')
    inputs.add_argument('-r2', '--r2', metavar='FILE',
                        help='R2 reads, for paired-end fastq files.')
    inputs.add_argument('-s', '--sample', metavar='NAME',
                        help='With -r1: sample name used in the reports. Default: the input file name, without '
                             'extension and read suffix (_R1, _1, ...).')

    output = parser.add_argument_group('output')
    output.add_argument('-o', '--output', metavar='FOLDER', required=True,
                        help='Folder to hold the reports. Created if needed.')
    output.add_argument('--xlsx', action='store_true',
                        help='Also write the legacy vSNP-style report, stat_alignment_summary_<date>.xlsx. Needs '
                             'the xlsxwriter package.')
    output.add_argument('--no-html', action='store_true',
                        help='Do not write brucellaid_report.html, the self-contained quality assurance report.')
    output.add_argument('--pdf', action='store_true',
                        help='Also write the quality assurance report as a PDF, brucellaid_report.pdf. Needs '
                             'the reportlab package ("pip install brucellaid[pdf]").')
    output.add_argument('--no-md5', action='store_true',
                        help='Do not checksum the input files for the reports. MD5 takes a second or two per '
                             'gigabyte.')
    output.add_argument('--operator', metavar='NAME',
                        help='Name to put on the quality assurance report. Default: the user account.')

    calling = parser.add_argument_group('identification and typing')
    calling.add_argument('--min-hashes', metavar='N', type=int, default=DEFAULT_MIN_HASHES,
                         help='Minimum number of the reference\'s 1000 min-hashes the reads must share for a '
                              'species to be reported. Default: {}.'.format(DEFAULT_MIN_HASHES))
    calling.add_argument('--min-kmer-copies', metavar='N', type=int, default=KMER_COPIES_AUTO,
                         help='Number of copies of a k-mer the reads must hold for MASH to sketch it, which '
                              'is how the sequencing errors are kept out of the identification. The right '
                              'value rises with the coverage: the default, 0, derives it from the number of '
                              'bases in the files (one more copy per 20x, from 2 to 20).')
    calling.add_argument('--min-reads', metavar='N', type=int, default=DEFAULT_MIN_READS,
                         help='Minimum number of reads assigned to an MLST allele for the allele to be called. '
                              'Default: {}.'.format(DEFAULT_MIN_READS))
    calling.add_argument('--no-mlst', action='store_true',
                         help='Identify the species only, and skip the MLST typing.')
    calling.add_argument('--no-st-species', action='store_true',
                         help='Do not let the sequence type take part in the species call: report the species of '
                              'the closest reference genome, whether or not another species is just as close.')
    calling.add_argument('--min-st-purity', metavar='FRACTION', type=float, default=DEFAULT_MIN_PURITY,
                         help='Share of the isolates of a sequence type that must carry the same species for '
                              'that ST to take part in the species call. Default: {}.'.format(DEFAULT_MIN_PURITY))
    calling.add_argument('--min-st-isolates', metavar='N', type=positive_int, default=DEFAULT_MIN_ISOLATES,
                         help='Isolates a sequence type must rest on to take part in the species call. '
                              'Default: {}.'.format(DEFAULT_MIN_ISOLATES))
    calling.add_argument('-k', '--kmer-size', metavar='N', type=positive_int, default=DEFAULT_KMER_SIZE,
                         help='K-mer size used to assign the reads to the MLST alleles. '
                              'Default: {}.'.format(DEFAULT_KMER_SIZE))
    calling.add_argument('-m', '--mismatch', metavar='N', type=int, default=DEFAULT_MISMATCHES,
                         help='Mismatches allowed between an allele k-mer and a read (maximum Hamming distance). '
                              'Default: {}, as the alleles of a locus differ by a few bases.'
                              .format(DEFAULT_MISMATCHES))

    data = parser.add_argument_group('reference data (the package ships all of it)')
    data.add_argument('--reference-set', choices=('default', 'extended'), default='default', dest='reference_set',
                      help='Which packaged reference set to identify against. "default" is the 16 genomes of '
                           'vSNP, 8 species; "extended" is 428 curated genomes covering the whole Brucella '
                           'clade, 29 of the 30 named species plus unnamed lineages. Default: default.')
    data.add_argument('--sketch', metavar='FILE',
                      help='MASH sketch of the reference genomes, instead of the one of the package.')
    data.add_argument('--references', metavar='FILE', dest='references',
                      help='Table of the reference genomes and their species, instead of the one of the package.')
    data.add_argument('--alleles', metavar='FILE',
                      help='Fasta file of the MLST alleles, instead of the one of the package.')
    data.add_argument('--profiles', metavar='FILE',
                      help='Table of the MLST profiles, instead of the one of the package.')
    data.add_argument('--st-species', metavar='FILE',
                      help='Table of the species of each sequence type, instead of the one of the package. '
                           'Built by scripts/make_st_species.py.')

    other = parser.add_argument_group('performance and other options')
    other.add_argument('-t', '--threads', metavar='N', type=positive_int, default=max_cpu,
                       help='Number of threads. Default: all available ({}).'.format(max_cpu))
    other.add_argument('-p', '-j', '--parallel', '--jobs', metavar='N', type=positive_int, default=1, dest='parallel',
                       help='With -i: number of samples analysed at the same time, sharing the threads. Default: 1.')
    other.add_argument('--memory', metavar='SIZE', type=java_memory,
                       help='Memory for Seal (Java heap size), per job, e.g. 2g. Default: from the size of the '
                            'allele file, within 85%% of the machine.')
    other.add_argument('--verbose', action='store_true', help='Show debug messages, including the commands run.')
    other.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    return parser


def check_arguments(parser, args):
    if args.input and (args.r2 or args.sample):
        parser.error('-r2 and --sample are only used with -r1')
    if args.r1 and args.parallel > 1:
        parser.error('--parallel is only used with -i')
    if args.min_kmer_copies < 0:
        parser.error('--min-kmer-copies cannot be negative (0 derives it from the coverage)')
    if args.min_hashes < 1:
        parser.error('--min-hashes must be >= 1')
    if args.min_reads < 1:
        parser.error('--min-reads must be >= 1')
    if args.mismatch < 0:
        parser.error('--mismatch cannot be negative')
    if not 0.0 < args.min_st_purity <= 1.0:
        parser.error('--min-st-purity must be a fraction between 0 and 1, not {}'.format(args.min_st_purity))
    if args.st_species and args.no_st_species:
        parser.error('--st-species gives a table to use and --no-st-species uses none: pick one')
    if args.mismatch >= args.kmer_size:
        parser.error('--mismatch ({}) must be smaller than --kmer-size ({})'.format(args.mismatch, args.kmer_size))
    if args.sample is not None and (not args.sample.strip() or '/' in args.sample or args.sample in ('.', '..')):
        parser.error('--sample must be a name, not a path: "{}"'.format(args.sample))
    if args.sample is not None and re.search(r'[\t\r\n]', args.sample):
        parser.error('--sample cannot contain tabs or line breaks (they would break the table)')


def file_bases(path):
    """Rough number of bases in a fasta file: its size, which is close enough to size the Seal heap."""
    try:
        return Path(path).stat().st_size
    except OSError:
        return 0


def reference_data(args, sketch, references, loci, st_table, mlst=True, md5=True):
    """
    The reference data files the run measured against, each with what it holds, for the reports.

    A quality assurance report that names "the MLST profile table" proves nothing. These are the files, with
    their checksums, so that two runs can be shown to have used the same data, or not.
    """
    alleles = alleles_fasta(args.alleles)
    files = [('MASH sketch', sketch, '{} reference genome(s), k = 21, 1000 hashes'.format(
                  len(sketch_references(sketch)))),
             ('Reference table', args.references or data_path(REFERENCES_FILE),
              '{} genome(s) with their species, biovar, strain and accessions'.format(len(references)))]
    if mlst:
        files.append(('MLST alleles', alleles, '{} allele(s) of the {} loci'.format(
            qa.fasta_count(alleles), len(loci))))
        files.append(('MLST profiles', args.profiles or data_path(PROFILES_FILE),
                      '{} sequence type(s)'.format(len(load_profiles(args.profiles)[1]))))
    if st_table is not None:
        files.append(('ST to species', args.st_species or data_path(ST_SPECIES_FILE),
                      '{} sequence type(s) with a species'.format(len(st_table))))
    return qa.data_facts(files, md5=md5)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    check_arguments(parser, args)
    setup_logging(args.verbose)
    command = shlex.join([parser.prog, *(sys.argv[1:] if argv is None else argv)])
    log.debug('brucellaid %s: %s', __version__, command)
    started = datetime.now().astimezone()
    start = time.monotonic()
    output = Path(args.output).expanduser()
    mlst = not args.no_mlst
    try:
        check_mash()
        if mlst:
            check_seal()
        if args.reference_set == 'extended':
            if args.sketch or args.references:
                raise BrucellaIDError('--reference-set extended names the packaged extended sketch and reference '
                                      'table: it cannot be combined with --sketch or --references.')
            args.sketch = data_path(EXTENDED_SKETCH_FILE)
            args.references = data_path(EXTENDED_REFERENCES_FILE)
        sketch = reference_sketch(args.sketch)
        references = load_references(args.references)
        alleles = alleles_fasta(args.alleles)
        loci = load_profiles(args.profiles)[0] if mlst else ()
        st_table = None if args.no_st_species else load_st_species(args.st_species)
        data = reference_data(args, sketch, references, loci, st_table, mlst=mlst, md5=not args.no_md5)
        log.info('%d reference genome(s) in %s', len(sketch_references(sketch)), sketch)
        for message in compare_with_sketch(sketch_references(sketch), references):
            log.warning('%s', message)
        if mlst:
            log.info('%d MLST locus/loci: %s', len(loci), ', '.join(loci))
        if st_table:
            log.info('%d sequence type(s) with a species, to decide the species MASH cannot', len(st_table))
        if args.input:
            samples = find_samples(args.input, exclude=[output])
            log.info('%d sample(s) found in %s', len(samples), args.input)
        else:
            samples = [single_sample(args.r1, args.r2, name=args.sample)]
        jobs = min(args.parallel, len(samples))
        memory = args.memory or auto_memory(file_bases(alleles), jobs=jobs)
        results = analyse_samples(samples, jobs=jobs, threads=args.threads, sketch=sketch, references=references,
                                  min_hashes=args.min_hashes, min_copies=args.min_kmer_copies,
                                  min_reads=args.min_reads, kmer_size=args.kmer_size,
                                  mismatch=args.mismatch, memory=memory, alleles=args.alleles,
                                  profiles=args.profiles, mlst=mlst, st_table=st_table,
                                  use_st_species=not args.no_st_species, min_st_purity=args.min_st_purity,
                                  min_st_isolates=args.min_st_isolates)
        run = run_info(command, {'Input': ' '.join(os.path.abspath(f) for f in (args.input, args.r1, args.r2) if f),
                                 'Output folder': str(output.resolve()),
                                 'Reference sketch': os.path.abspath(sketch),
                                 'MLST alleles': os.path.abspath(alleles) if mlst else None,
                                 'Minimum shared hashes': args.min_hashes,
                                 'Minimum k-mer copies': (args.min_kmer_copies or
                                                          'from the coverage of each sample'),
                                 'Minimum reads per allele': args.min_reads if mlst else None,
                                 'ST to species table': None if st_table is None else len(st_table),
                                 'Minimum ST purity': None if st_table is None else args.min_st_purity,
                                 'Minimum ST isolates': None if st_table is None else args.min_st_isolates,
                                 'K-mer size': args.kmer_size, 'Mismatches': args.mismatch,
                                 'MLST': mlst, 'Threads': args.threads, 'Parallel samples': jobs,
                                 'Seal memory': memory,
                                 'Input checksums': 'none' if args.no_md5 else 'MD5',
                                 'HTML report': not args.no_html, 'PDF report': args.pdf},
                       started, operator=args.operator, data=data)
        rows, reports = write_reports(results, references, output, started, command=command, run=run,
                                      xlsx=args.xlsx, html=not args.no_html, pdf=args.pdf,
                                      md5=not args.no_md5)
    except USER_ERRORS as e:
        log.error('%s', e)
        log.debug('Traceback:', exc_info=True)  # Shown with --verbose
        sys.exit(1)
    except KeyboardInterrupt:
        log.error('Interrupted')
        sys.exit(130)

    # The results go to stdout, the messages to stderr: "brucellaid ... > results.tsv" works
    print_table(rows)
    for report in reports:
        log.info('Report saved in %s', report)
    log.info('Done in %s', format_elapsed(time.monotonic() - start))
    failed = [r.sample for r in results if r.error]
    if failed:
        log.error('%d of %d sample(s) failed: %s', len(failed), len(results), ', '.join(failed))
        sys.exit(1)
