from accessoryFunctions.accessoryFunctions import filer, run_subprocess, write_to_logfile, make_path, relative_symlink
from glob import glob
import xlsxwriter
import os


class BrucellaIDMethods(object):

    @staticmethod
    def file_list(path):
        """
        Create a list of all FASTQ files present in the supplied path. Accepts .fastq.gz, .fastq, .fq, and .fq.gz
        extensions only
        :param path: type STR: absolute path of folder containing FASTQ files
        :return fastq_files: sorted list of all FASTQ files present in the supplied path
        """
        # Use glob to find the acceptable extensions of FASTQ files in the supplied path
        fastq_files = glob(os.path.join(path, '*.fastq'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fastq.gz'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fq'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fq.gz'))
        # Sort the list of fastq files
        fastq_files = sorted(fastq_files)
        # Ensure that there are actually files present in the path
        assert fastq_files, 'Cannot find FASTQ files in the supplied path: {path}'.format(path=path)

        return fastq_files

    @staticmethod
    def strain_list(fastq_files):
        """
        Use the filer method to parse a list of absolute paths of FASTQ files to yield the name of the strain.
        e.g. /path/to/03-1057_S10_L001_R1_001.fastq.gz will return /path/to/03-1057:
        /path/to/03-1057_S10_L001_R1_001.fastq.gz
        :param fastq_files: type LIST: List of absolute of paired and/or unpaired FASTQ files
        :return strain_dict: dictionary of the absolute paths to the base strain name: FASTQ files
        """
        # As filer returns a set of the names, transform this set into a sorted list
        strain_dict = filer(filelist=fastq_files, returndict=True)
        return strain_dict

    @staticmethod
    def strain_namer(strain_folders):
        """
        Extract the base strain name from a list of absolute paths with the strain name (this list is usually created
        using the strain_list method above
        e.g. /path/to/03-1057 will yield 03-1057
        :param strain_folders: type iterable: List or dictionary of absolute paths and base strain names
        :return: strain_names: list of strain base names
        """
        # Initialise a dictionary to store the strain names
        strain_name_dict = dict()
        for strain_folder in strain_folders:
            # Extract the base name from the absolute path plus base name
            strain_name = os.path.basename(strain_folder)
            strain_name_dict[strain_name] = strain_folder
        return strain_name_dict

    @staticmethod
    def file_link(strain_folder_dict, strain_name_dict):
        """
        Create folders for each strain. Create relative symlinks to the original FASTQ files from within the folder
        :param strain_folder_dict: type DICT: Dictionary of strain folder path: FASTQ files
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :return: strain_fastq_dict: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        """
        #
        strain_fastq_dict = dict()
        for strain_name, strain_folder in strain_name_dict.items():
            # Create the strain folder path if required
            make_path(strain_folder)
            # Use the strain_folder value from the strain_name_dict as the key to extract the list of FASTQ files
            # associated with each strain
            for fastq_file in strain_folder_dict[strain_folder]:
                # Create relative symlinks between the original FASTQ files and the strain folder
                symlink_path = relative_symlink(src_file=fastq_file,
                                                output_dir=strain_folder,
                                                export_output=True)
                # Add the absolute path of the symlink to the dictionary
                try:
                    strain_fastq_dict[strain_name].append(symlink_path)
                except KeyError:
                    strain_fastq_dict[strain_name] = [symlink_path]
        return strain_fastq_dict

    @staticmethod
    def call_mash_sketch(strain_fastq_dict, strain_name_dict, logfile):
        """
        Run MASH sketch on the provided FASTQ files
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :param logfile: type STR: Absolute path to logfile basename
        :return: fastq_sketch_dict: Dictionary of strain name: absolute path to MASH sketch file
        """
        # Initialise a dictionary to store the absolute path of the sketch file
        fastq_sketch_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            # Extract the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute paths of the sketch file with and without the .msh extension (used for calling MASH)
            fastq_sketch_no_ext = os.path.join(strain_folder, '{sn}_sketch'.format(sn=strain_name))
            fastq_sketch = fastq_sketch_no_ext + '.msh'
            # Create the system call - cat together the FASTQ files, and pipe them into MASH
            # -p requests the number of desired threads, -m sets the minimum copies of each k-mer required to pass
            # noise filter for reads to 2 (ignores single copy kmers), - indicates that MASH should use stdin
            # -o is the absolute path to the sketch output file
            mash_sketch_command = 'cat {fastq} | mash sketch -m 2 - -o {output_file}' \
                .format(fastq=' '.join(fastq_files),
                        output_file=fastq_sketch_no_ext)
            # Only make the system call if the output sketch file doesn't already exist
            if not os.path.isfile(fastq_sketch):
                out, err = run_subprocess(command=mash_sketch_command)
                # Write the stdout, and stderr to the main logfile, as well as to the strain-specific logs
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the absolute path of the sketch file
            fastq_sketch_dict[strain_name] = fastq_sketch
        return fastq_sketch_dict

    @staticmethod
    def call_mash_dist(strain_fastq_dict, strain_name_dict, fastq_sketch_dict, ref_sketch_file, logfile):
        """
        Run a MASH dist of a pre-sketched set of FASTQ reads against the custom MASH sketch file of the reference
        genomes supported by vSNP
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :param fastq_sketch_dict: type DICT: Dictionary of strain name: absolute path to MASH sketch file
        :param ref_sketch_file: type STR: Absolute path to the custom sketch file of reference sequences
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_mash_outputs: Dictionary of strain name: absolute path of MASH dist output table
        """
        # Initialise the dictionary to store the absolute path of the MASH dist output file
        strain_mash_outputs = dict()
        for strain_name in strain_fastq_dict:
            # Extract the absolute path of the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Extract the absolute path of the strain-specific sketch file
            fastq_sketch_file = fastq_sketch_dict[strain_name]
            # Set the absolute path of the MASH dist output table
            out_tab = os.path.join(strain_folder, '{sn}_mash.tab'.format(sn=strain_name))
            # Create the system call: -p is the number of threads requested, -m sets the minimum copies of each k-mer
            # required to pass noise filter for reads to 2 (ignores single copy kmers). MASH outputs are piped to
            # the sort function, which sorts the data as follows: g: general numeric sort, K: keydef
            mash_dist_command = 'mash dist -m 2 {ref_sketch_file} {fastq_sketch} | sort -gk3 > {out}' \
                .format(ref_sketch_file=ref_sketch_file,
                        fastq_sketch=fastq_sketch_file,
                        out=out_tab)
            if not os.path.isfile(out_tab):
                out, err = run_subprocess(command=mash_dist_command)
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            strain_mash_outputs[strain_name] = out_tab
        return strain_mash_outputs

    @staticmethod
    def parse_mash_accession_species(mash_species_file):
        """
        Parse the reference genus accession: species code .csv file included in the mash dependencies path
        :param mash_species_file: type STR: Absolute path to file containing reference file name: species code
        :return: accession_species_dict: Dictionary of reference accession: species code
        """
        # Initialise a dictionary to store the species code
        accession_species_dict = dict()
        with open(mash_species_file, 'r') as species_file:
            for line in species_file:
                # Extract the accession and the species code pair from the line
                accession, species = line.rstrip().split(',')
                # Populate the dictionary with the accession: species pair
                accession_species_dict[accession] = species
        return accession_species_dict

    @staticmethod
    def mash_best_ref(mash_dist_dict, accession_species_dict):
        """
        Parse the MASH dist output table to determine the closest reference sequence, as well as the total
        number of matching hashes the strain and that reference genome share
        :param mash_dist_dict: type DICT: Dictionary of strain name: absolute path of MASH dist output table
        :param accession_species_dict: type DICT: Dictionary of reference accession: species code
        :return: strain_best_ref_dict: Dictionary of strain name: closest MASH-calculated reference genome
        :return: strain_ref_matches_dict: Dictionary of strain name: number of matching hashes between query and
        closest reference genome
        :return: strain_species_dict: Dictionary of strain name: species code
        """
        # Initialise dictionaries to store the strain-specific closest reference genome, number of matching hashes
        # between read sets and the reference genome, as well as the species code
        strain_best_ref_dict = dict()
        strain_ref_matches_dict = dict()
        strain_species_dict = dict()
        for strain_name, mash_dist_table in mash_dist_dict.items():
            with open(mash_dist_table, 'r') as mash_dist:
                # Extract all the data included on each line of the table outputs
                best_ref, query_id, mash_distance, p_value, matching_hashes = mash_dist.readline().rstrip().split('\t')
            # Split the total of matching hashes from the total number of hashes
            matching_hashes = int(matching_hashes.split('/')[0])
            # Populate the dictionaries appropriately
            if matching_hashes >= 150:
                strain_best_ref_dict[strain_name] = best_ref
                strain_ref_matches_dict[strain_name] = matching_hashes
                strain_species_dict[strain_name] = accession_species_dict[best_ref]
        return strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict

    @staticmethod
    def brucella_mlst(seqpath, mlst_db_path, outpath, logfile):
        """
        Run MLST analyses on the strains
        :param seqpath: type STR: Absolute path to folder containing FASTQ files
        :param mlst_db_path: type STR: Absolute path to folder containing Brucella MLST database files
        :param logfile: type STR: Absolute path to the logfile basename
        """
        # Set the system call to the MLST python package in sipprverse
        # -s path of input fastq
        # -t path of target files to process
        # -a Specify analysis type: mlst or rmlst
        #
        mlst_cmd = 'python -m MLSTsippr.mlst -s {} -t {} -a mlst {}' \
            .format(seqpath, mlst_db_path, outpath)
        # Run the system call - don't worry about checking if there are outputs, the package will parse these reports
        # rather than re-running the analyses
        out, err = run_subprocess(mlst_cmd)
        # Write STDOUT and STDERR to the logfile
        write_to_logfile(out=out,
                         err=err,
                         logfile=logfile)

    @staticmethod
    def parse_mlst_report(strain_name_dict, mlst_report):
        """
        Parse the MLST report to create a dictionary with the calculated sequence type and the number of exact matches
        per strain to the sequence type
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path to strain working directory
        :param mlst_report: type STR: Absolute path to MLST report
        :return: strain_mlst_dict: Dictionary of strain name: sequence type and number of matches to the sequence type
        """
        # Initialise a dictionary to store the parsed MLST results
        strain_mlst_dict = dict()
        # Open the MLST report file
        with open(mlst_report, 'r') as report:
            # Skip the header line
            next(report)
            for line in report:
                # Split the line on commas
                data = line.rstrip().split(',')
                # Extract the sequence type and the number of matches to the sequence type from the line. Populate the
                # dictionary with these values
                strain_mlst_dict[data[0]] = {
                    'sequence_type': data[2],
                    'matches': data[3],
                }
        # If the strain did not have MLST outputs, populate negative 'ND' values for the sequence type and number
        # of matches
        for strain in strain_name_dict:
            # Only populate negative values for strains absent from the outputs
            if strain not in strain_mlst_dict:
                strain_mlst_dict[strain] = {
                    'sequence_type': 'ND',
                    'matches': 'ND',
                }
        return strain_mlst_dict

    @staticmethod
    def create_vcf_report(start_time, strain_species_dict, strain_mlst_dict, report_path):
        """
        Create an Excel report of the vcf outputs
        :param start_time: type datetime.now(): Datetime object
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :param strain_mlst_dict: type DICT: Dictionary of strain name: sequence type and number of matches to the
        sequence type
        :param report_path: type STR: Absolute path to path in which reports are to be created
        :return: vcf_report: Absolute path to the Excel report
        """
        # Create a date string consistent with classic vSNP
        start = '{:%Y-%m-%d_%H-%M-%S}'.format(start_time)
        # Initialise a string to store the absolute path to the Excel report
        vcf_report = os.path.join(os.path.join(report_path, 'stat_alignment_summary_{start}.xlsx'
                                               .format(start=start)))
        # Set the list of the headers
        header_list = ['time_stamp', 'sample_name', 'species', 'mlst_type']
        # Create a workbook to store the report using xlsxwriter.
        workbook = xlsxwriter.Workbook(vcf_report)
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        # Set the column to zero
        col = 0
        # Write the header to the spreadsheet
        for header in header_list:
            worksheet.write(row, col, header, bold)
            col += 1
        for strain_name, info in strain_mlst_dict.items():
            # Increment the row and reset the column to zero in preparation of writing results
            row += 1
            col = 0
            # Allow strains that do not have a match to be added to the report
            try:
                strain_species = strain_species_dict[strain_name]
                ml_seq_type = strain_mlst_dict[strain_name]['sequence_type']
            except KeyError:
                strain_species = 'ND'
                ml_seq_type = 'ND'
            # Populate the list with the required data
            data_list = [start, strain_name, strain_species, ml_seq_type]
            # Write out the data to the spreadsheet
            for results in data_list:
                worksheet.write(row, col, results, courier)
                col += 1
        # Close the workbook
        workbook.close()
        return vcf_report
