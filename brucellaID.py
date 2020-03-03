#!/usr/bin/env python3

from argparse import ArgumentParser
from accessoryFunctions.accessoryFunctions import SetupLogging
import logging
import os
from datetime import datetime
from brucellaID_methods import BrucellaIDMethods


class BrucellaID(object):

    def main(self):
        self.fastq_manipulation()
        self.best_reference_calculation()
        self.typing()
        self.report()

    def fastq_manipulation(self):
        """
        Determine the number of strains to process. Create strain-specific working directories with relative symlinks
        to FASTQ files
        """
        logging.info('Locating FASTQ files, creating strain-specific working directories and symlinks to files')
        fastq_files = BrucellaIDMethods.file_list(path=self.seqpath)
        logging.debug('FASTQ files: \n{fastq_files}'.format(fastq_files='\n'.join(fastq_files)))
        strain_folder_dict = BrucellaIDMethods.strain_list(fastq_files=fastq_files)
        self.strain_name_dict = BrucellaIDMethods.strain_namer(strain_folders=strain_folder_dict)
        logging.debug(
            'Strain names: \n{strain_names}'.format(strain_names='\n'.join(sorted(self.strain_name_dict))))
        self.strain_fastq_dict = BrucellaIDMethods.file_link(strain_folder_dict=strain_folder_dict,
                                                             strain_name_dict=self.strain_name_dict)
        logging.debug(
            'Strain-specific symlinked FASTQ files: \n{symlinks}'.format(
                symlinks='\n'.join(['{strain_name}: {fastq_files}'.format(strain_name=sn, fastq_files=ff)
                                    for sn, ff in self.strain_fastq_dict.items()])))

    def best_reference_calculation(self):
        """
                Determine the best reference using MASH
                """
        logging.info('Running MASH analyses')
        fastq_sketch_dict = BrucellaIDMethods.call_mash_sketch(strain_fastq_dict=self.strain_fastq_dict,
                                                               strain_name_dict=self.strain_name_dict,
                                                               logfile=self.logfile)
        logging.debug(
            'Strain-specific MASH sketch files: \n{files}'.format(
                files='\n'.join(['{strain_name}: {sketch_file}'.format(strain_name=sn, sketch_file=sf)
                                 for sn, sf in fastq_sketch_dict.items()])))
        logging.info('Parsing MASH outputs to determine closest reference genomes')
        mash_dist_dict = BrucellaIDMethods.call_mash_dist(strain_fastq_dict=self.strain_fastq_dict,
                                                          strain_name_dict=self.strain_name_dict,
                                                          fastq_sketch_dict=fastq_sketch_dict,
                                                          ref_sketch_file=os.path.join(
                                                              self.dependency_path, 'mash', 'vsnp_reference.msh'),
                                                          logfile=self.logfile)
        logging.debug(
            'Strain-specific MASH output tables: \n{files}'.format(
                files='\n'.join(['{strain_name}: {table}'.format(strain_name=sn, table=tf)
                                 for sn, tf in mash_dist_dict.items()])))
        logging.info('Loading reference genome: species dictionary')
        accession_species_dict = BrucellaIDMethods.parse_mash_accession_species(mash_species_file=os.path.join(
            self.dependency_path, 'mash', 'species_accessions.csv'))
        logging.info('Determining closest reference genome and extracting corresponding species from MASH outputs')
        self.strain_best_ref_dict, self.strain_ref_matches_dict, self.strain_species_dict = \
            BrucellaIDMethods.mash_best_ref(mash_dist_dict=mash_dist_dict,
                                            accession_species_dict=accession_species_dict)
        logging.debug(
            'Strain-specific MASH-calculated best reference file: \n{files}'.format(
                files='\n'.join(['{strain_name}: {best_ref}'.format(strain_name=sn, best_ref=br)
                                 for sn, br in self.strain_best_ref_dict.items()])))
        logging.debug(
            'Number of matches to strain-specific MASH-calculated best reference file: \n{files}'.format(
                files='\n'.join(['{strain_name}: {num_matches}'.format(strain_name=sn, num_matches=nm)
                                 for sn, nm in self.strain_ref_matches_dict.items()])))
        logging.debug(
            'Species code for best reference file: \n{files}'.format(
                files='\n'.join(['{strain_name}: {species_code}'.format(strain_name=sn, species_code=sc)
                                 for sn, sc in self.strain_species_dict.items()])))

    def typing(self):
        logging.info('Performing MLST analyses')
        BrucellaIDMethods.brucella_mlst(seqpath=self.seqpath,
                                        mlst_db_path=self.dbpath,
                                        outpath=self.outpath,
                                        logfile=self.logfile)
        logging.info('Parsing MLST outputs')
        self.strain_mlst_dict = BrucellaIDMethods.parse_mlst_report(strain_name_dict=self.strain_name_dict,
                                                                    mlst_report=os.path.join(self.outpath,
                                                                                             'reports', 'mlst.csv'))
        logging.debug('MLST results: \n{files}'.format(
            files='\n'.join(['{strain_name}: {mlst_result}'.format(strain_name=sn, mlst_result=mr)
                             for sn, mr in self.strain_mlst_dict.items()])))

    def report(self):
        """
        Create the .xlsx report consistent with the legacy vSNP format
        """
        logging.info('Creating report')
        BrucellaIDMethods.create_vcf_report(
            start_time=self.start_time,
            strain_species_dict=self.strain_species_dict,
            strain_mlst_dict=self.strain_mlst_dict,
            report_path=self.report_path)

    def __init__(self, args):
        SetupLogging()

        self.seqpath = args.sequencepath
        self.dbpath = args.databasepath
        self.outpath = args.outputpath

        # Initialise class variables
        self.strain_name_dict = dict()
        self.strain_species_dict = dict()
        self.strain_mlst_dict = dict()
        self.strain_fastq_dict = dict()
        self.strain_best_ref_dict = dict()
        self.strain_ref_matches_dict = dict()

        # Extract the path of the folder containing this script
        self.script_path = os.path.abspath(os.path.dirname(__file__))
        # Use the script path to set the absolute path of the dependencies folder
        self.dependency_path = os.path.join(self.script_path, 'dependencies')
        assert os.path.isdir(self.dependency_path), 'Something went wrong with the install. Cannot locate the ' \
                                                    'dependencies folder in: {sp}'.format(sp=self.script_path)
        self.logfile = os.path.join(self.outpath, 'log')
        self.report_path = os.path.join(self.outpath, 'reports')
        self.start_time = datetime.now()

        self.main()


if __name__ == '__main__':
    parser = ArgumentParser(description='Perform in silico PCR using bbduk and SPAdes')

    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of folder containing .fasta/.fastq(.gz) files to process.')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='Path of folder containing the MLST database')
    parser.add_argument('-o', '--outputpath',
                        required=True,
                        help='Path of folder to hold output files')

    # Get the arguments into an object
    arguments = parser.parse_args()

    BrucellaID(arguments)
