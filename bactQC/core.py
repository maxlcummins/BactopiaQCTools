# bactQC/core.py

import os
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import json
import logging
import re

# Configure logging at the module level
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Genome:
    def __init__(self, sample_name, input_dir, taxid=None):
        self.sample_name = sample_name
        self.input_dir = input_dir
        self.taxid = str(taxid) if taxid is not None else None
        self.qc_data = {'sample': sample_name}
        self.qc_results = {'sample': sample_name}
        self.qc_requirements = {'sample': sample_name}
        
    def run(self):
        """
        Run all quality control checks for the sample.

        Populates self.qc_data with the results.
        """
        self.get_expected_genome_size()
        self.get_assembly_size()
        self.check_bracken()
        expected_genus = self.qc_data['genome_size']['organism_name'].split()[0]
        self.check_mlst(expected_genus)
        self.check_checkm()
        self.check_assembly_scan()
        self.check_fastp()

    def get_expected_genome_size(self):
        """
        Retrieve expected genome size information from the NCBI API or infer it based on bracken abundance data if taxid is not provided.

        Updates self.qc_data['genome_size'] with the results.
        """
        base_URL = "https://api.ncbi.nlm.nih.gov/genome/v0/expected_genome_size/expected_genome_size?species_taxid="

        if self.taxid is None:
            bracken_path = os.path.join(self.input_dir, self.sample_name, 'tools', 'bracken', f"{self.sample_name}.bracken.adjusted.abundances.txt")
            bracken_result = pd.read_csv(bracken_path, sep='\t')
            self.taxid = str(bracken_result.iloc[0]['taxonomy_id'])
            guessed_species_name = bracken_result.iloc[0]['name']
        else:
            self.taxid = str(self.taxid)

        url = f"{base_URL}{self.taxid}"
        response = requests.get(url)

        data = {}
        
        if response.content:
            try:
                root = ET.fromstring(response.content)
                data = {
                    'organism_name': root.findtext('organism_name'),
                    'species_taxid': root.findtext('species_taxid'),
                    'expected_ungapped_length': int(root.findtext('expected_ungapped_length')),
                    'minimum_ungapped_length': int(root.findtext('minimum_ungapped_length')),
                    'maximum_ungapped_length': int(root.findtext('maximum_ungapped_length'))
                }
            except ET.ParseError as e:
                print(f"Error parsing XML for taxid {self.taxid}: {e}")
        else:
            print(f"No content returned for taxid {self.taxid}")
        if data:
            self.qc_data['genome_size'] = data
        else:
            raise ValueError(f"Failed to retrieve genome size for taxid {self.taxid}")

    def get_assembly_size(self):
        """
        Retrieves the total contig length in base-pairs for a given sample from the assembler results.

        Updates self.qc_data['assembly_size'] with the results.
        """
        # Construct the file path to the assembler results
        file_path = os.path.join(self.input_dir, self.sample_name, 'main', 'assembler', f"{self.sample_name}.tsv")

        # Check if the file exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Assembly results file not found for sample '{self.sample_name}' at '{file_path}'")

        # Read the assembler results file
        assembly_scan_df = pd.read_csv(file_path, sep='\t')

        # Ensure the DataFrame has exactly one row
        if assembly_scan_df.shape[0] != 1:
            raise ValueError(f"Sample '{self.sample_name}': Expected one row in assembly scan data, but got {assembly_scan_df.shape[0]} rows")

        # Extract the total_contig_length from the first (and only) row
        total_length = assembly_scan_df.iloc[0]['total_contig_length']

        self.qc_data['assembly_size'] = {'total_length': total_length}
        
        
    def check_bracken(self, min_primary_abundance=0.80):
        """
        Check Bracken results for a given sample.

        Updates self.qc_data['bracken'] with the results.
        """
        # Get the file path
        file_path = os.path.join(self.input_dir, self.sample_name, 'tools', 'bracken', f"{self.sample_name}.bracken.tsv")
        
        # Check the path exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Bracken data not found at {file_path}")

        # Read the file
        bracken_result = pd.read_csv(file_path, sep='\t')

        # Check that bracken_result has only one row
        if bracken_result.shape[0] != 1:
            raise ValueError(f"Sample {self.sample_name}: Expected one row for bracken, but got {bracken_result.shape[0]} rows")

        # Add a column with a boolean for whether bracken_primary_species_abundance is greater than min_primary_abundance for first row
        bracken_result['passed_bracken_QC'] = bracken_result['bracken_primary_species_abundance'] > min_primary_abundance

        # Take the first row and convert it to a dictionary after dropping the 'sample' column
        bracken_result = bracken_result.iloc[0].drop('sample').to_dict()

        # Add an entry with the minimum primary abundance
        bracken_result['primary_abundance_requirement'] = min_primary_abundance

        # Check if the first word within bracken_primary_species differs from the first word within the bracken_secondary_species
        has_secondary_species  = bracken_result['bracken_secondary_species'] != 'No secondary abundance > 1%'

        # Check if our primary and secondary genera are Escherichia and Shigella (in either order)
        is_ecoli_or_shigella  = {'Escherichia', 'Shigella'} == {bracken_result['bracken_primary_species'], bracken_result['bracken_secondary_species']}

        # Determine if there is a conflict in the genera detected
        if has_secondary_species and is_ecoli_or_shigella:
            bracken_result['genus_conflict'] = bracken_result['bracken_primary_species'].split()[0] != bracken_result['bracken_secondary_species'].split()[0]
        else:
            bracken_result['genus_conflict'] = False

        self.qc_data['bracken'] = bracken_result
        
        self.qc_results['Passed bracken'] = bracken_result['passed_bracken_QC']
        
        self.qc_requirements['bracken'] = {'min_primary_abundance': min_primary_abundance}

    def check_mlst(self, expected_genus):
        """
        Check MLST results for a given sample.

        Updates self.qc_data['mlst'] with the results.
        """
        # Get the file path
        file_path = os.path.join(self.input_dir, self.sample_name, 'tools', 'mlst', f"{self.sample_name}.tsv")
        
        # Check the path exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"MLST data not found at {file_path}")

        # Read the file
        mlst_result = pd.read_csv(file_path, sep='\t', header=None)

        # Define the column names
        mlst_result.columns = ['sample', 'scheme', 'ST', 'allele1', 'allele2', 'allele3', 'allele4', 'allele5', 'allele6', 'allele7']

        # Remove '.fna.gz' or '.fna' from the sample column
        mlst_result['sample'] = mlst_result['sample'].str.replace(r'\.fna(\.gz)?$', '', regex=True)
        
        #Remove '.fasta.gz' or '.fasta' from the sample column
        mlst_result['sample'] = mlst_result['sample'].str.replace(r'\.fasta(\.gz)?$', '', regex=True)
        
        #Remove '.fna.gz' or '.fna' from the sample column
        mlst_result['sample'] = mlst_result['sample'].str.replace(r'\.fa(\.gz)?$', '', regex=True)

        # Check that mlst_result has only one row
        if mlst_result.shape[0] != 1:
            raise ValueError(f"Sample {self.sample_name}: Expected one row in MLST data, but got {mlst_result.shape[0]} rows")

        # Download the scheme species map if necessary
        cache_file = os.path.join(self.input_dir, 'scheme_species_map.tab')
        if not os.path.exists(cache_file):
            scheme_species_map = pd.read_csv('https://raw.githubusercontent.com/tseemann/mlst/master/db/scheme_species_map.tab', sep='\t')
            scheme_species_map.to_csv(cache_file, sep='\t', index=False)
        else:
            scheme_species_map = pd.read_csv(cache_file, sep='\t')

        # Filter the dataframe to only the expected genus
        scheme_species_map = scheme_species_map[scheme_species_map['GENUS'] == expected_genus]

        # Check if the scheme in mlst result partially matches with one or more schemes for the expected genus
        if not any(mlst_result['scheme'].str.contains('|'.join(scheme_species_map['#SCHEME']))):
            # Set mlst_result['passed_mlst'] to False
            mlst_result['passed_mlst'] = False
        else:
            # Set mlst_result['passed_mlst'] to True

            mlst_result['passed_mlst'] = True

        # Take the first row and convert it to a dictionary after dropping the 'sample' column
        mlst_result = mlst_result.iloc[0].drop('sample').to_dict()

        self.qc_data['mlst'] = mlst_result
                
        self.qc_results['Passed mlst'] = mlst_result['passed_mlst']
        
        self.qc_requirements['mlst']= {'expected_genus': expected_genus}

    def check_checkm(self, min_completeness=80, max_contamination=10):
        """
        Checks CheckM results for a given sample and evaluates quality metrics.

        Updates self.qc_data['checkm'] with the results.
        """
        # Search recursively for all files in the input directory called 'checkm.tsv'
        checkm_files = []
        for root, dirs, files in os.walk(self.input_dir):
            for file in files:
                if file == 'checkm.tsv':
                    checkm_files.append(os.path.join(root, file))

        # Check that we found at least one checkm.tsv file
        if len(checkm_files) == 0:
            raise ValueError(f"Bracken data not found at {file_path}")

        # Determine the most recent checkm.tsv file based on the timestamp name of directory two levels up
        checkm_files = sorted(checkm_files, key=lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))))

        # Select the most recent checkm.tsv file
        file_path = checkm_files[-1]

        # Read the file
        checkm_result = pd.read_csv(file_path, sep='\t')

        # Rename 'Bin Id' to 'sample'
        checkm_result.rename(columns={'Bin Id': 'sample'}, inplace=True)

        # Check if sample name exists in the checkm file
        if self.sample_name not in checkm_result['sample'].values:
            raise ValueError(f"Sample {self.sample_name} not found in {file_path}")

        # Filter the dataframe to only the sample of interest
        checkm_result = checkm_result[checkm_result['sample'] == self.sample_name]

        # Create a dictionary from the first row after dropping the 'sample' column
        checkm_result = checkm_result.iloc[0].drop('sample').to_dict()

        # Add an entry with the minimum completeness
        checkm_result['completeness_requirement'] = min_completeness

        # Add an entry with the maximum contamination
        checkm_result['contamination_requirement'] = max_contamination

        # Add an entry as a boolean for if the completeness is greater than min_completeness
        checkm_result['passed_completeness'] = checkm_result['Completeness'] > min_completeness

        # Add an entry as a boolean for if the contamination is less than max_contamination
        checkm_result['passed_contamination'] = checkm_result['Contamination'] < max_contamination

        # Add an entry as a boolean for if completion and contamination requirements are met
        checkm_result['passed_checkm_QC'] = checkm_result['passed_completeness'] and checkm_result['passed_contamination']

        self.qc_data['checkm'] = checkm_result
        
        self.qc_results['Passed checkm'] = checkm_result['passed_checkm_QC']
        
        self.qc_requirements['checkm'] = {'max_contamination': max_contamination,
                                          'min_completeness': min_completeness}

    def check_assembly_scan(self, maximum_contigs=500, minimum_N50=15000):
        """
        Check the quality of assembly scan results for a given sample.

        Updates self.qc_data['assembly_scan'] with the results.
        """
        # Ensure genome size data is available
        if 'genome_size' not in self.qc_data:
            self.get_expected_genome_size()

        data = self.qc_data['genome_size']

        # Get the file path
        file_path = os.path.join(self.input_dir, self.sample_name, 'main', 'assembler', f"{self.sample_name}.tsv")

        # Check the path exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Assembly scan data not found at {file_path}")

        # Read the file
        assembly_scan_df = pd.read_csv(file_path, sep='\t')

        if assembly_scan_df.shape[0] != 1:
            raise ValueError(f"Sample {self.sample_name}: Expected one row in assembly scan data, but got {assembly_scan_df.shape[0]} rows")

        # Take the first row and convert it to a dictionary after dropping the 'sample' column
        assembly_scan_results = assembly_scan_df.iloc[0].drop('sample').to_dict()

        # Add QC pass/fail flags
        assembly_scan_results['passed_contigs'] = assembly_scan_results['total_contig'] < maximum_contigs
        assembly_scan_results['passed_N50'] = assembly_scan_results['n50_contig_length'] > minimum_N50
        assembly_scan_results.update(data)

        # Calculate acceptable genome size range
        min_length = data['minimum_ungapped_length'] * 0.95
        max_length = data['maximum_ungapped_length'] * 1.05
        total_length = assembly_scan_results['total_contig_length']

        assembly_scan_results['passed_genome_size'] = (total_length > min_length) and (total_length < max_length)
        assembly_scan_results['passed_assembly_scan'] = all([
            assembly_scan_results['passed_contigs'],
            assembly_scan_results['passed_N50'],
            assembly_scan_results['passed_genome_size']
        ])

        self.qc_data['assembly_scan'] = assembly_scan_results
        
        self.qc_results['Passed assembly_scan'] = assembly_scan_results['passed_assembly_scan']
        
        self.qc_requirements['assembly_scan'] = {'maximum_contigs': maximum_contigs,
                                                 'minimum_N50': minimum_N50,
                                                 'minimum_ungapped_length': min_length,
                                                 'maximum_ungapped_length': max_length}

    def check_fastp(self, min_q30_bases=0.90, min_coverage=30):
        """
        Checks fastp quality control data for a given sample.

        Updates self.qc_data['fastp'] with the results.
        """

        # Get the file path
        file_path = os.path.join(self.input_dir, self.sample_name, 'main', 'qc', 'summary', f"{self.sample_name}.fastp.json")
        
        # Check the path exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Fastp data not found at {file_path}")
        
        # Read the JSON file
        with open(file_path, 'r') as file:
            data = json.load(file)

            # Create a dictionary with the name as a single-element list and set it as the index
            fastp_results = {}

            fastp_results['pre_filt_total_reads'] = int(data['summary']['before_filtering']['total_reads'])
            fastp_results['pre_filt_total_bases'] = int(data['summary']['before_filtering']['total_bases'])
            fastp_results['pre_filt_q30_rate'] = float(data['summary']['before_filtering']['q30_rate'])
            fastp_results['pre_filt_gc'] = float(data['summary']['before_filtering']['gc_content'])
            fastp_results['post_filt_total_reads'] = int(data['summary']['after_filtering']['total_reads'])
            fastp_results['post_filt_total_bases'] = int(data['summary']['after_filtering']['total_bases'])
            fastp_results['post_filt_q20_rate'] = float(data['summary']['after_filtering']['q20_rate'])
            fastp_results['post_filt_q30_rate'] = float(data['summary']['after_filtering']['q30_rate'])
            fastp_results['post_filt_gc'] = float(data['summary']['after_filtering']['gc_content'])

        # Ensure assembly size is available
        if 'assembly_size' not in self.qc_data:
            self.get_assembly_size()

        assembly_size = self.qc_data['assembly_size']

        # Compute coverage
        coverage = fastp_results['post_filt_total_bases'] / assembly_size['total_length']

        # Add the coverage to the dictionary
        fastp_results['coverage'] = int(coverage)

        # Add a boolean for whether the Q30 bases are greater than min_q30_bases
        fastp_results['passed_q30_bases'] = fastp_results['post_filt_q30_rate'] > min_q30_bases

        # Add a boolean for whether the coverage is greater than min_coverage
        fastp_results['passed_coverage'] = coverage > min_coverage

        # Add a boolean for whether the fastp QC passed
        fastp_results['passed_fastp_QC'] = fastp_results['passed_q30_bases'] and fastp_results['passed_coverage']      

        self.qc_data['fastp'] = fastp_results
        
        self.qc_results['Passed fastp'] = fastp_results['passed_fastp_QC']
        
        self.qc_requirements['fastp'] = {'min_q30_bases': min_q30_bases,
                                         'min_coverage': min_coverage}

    def get_qc_results(self):
        """
        Returns the quality control results.
        """
        # Create a DataFrame from the qc_results dictionary
        results_df = pd.DataFrame.from_dict(self.qc_results, orient='index')
        
        # Transpose the DataFrame
        results_df = results_df.T
        
        # Add columns for detected species
        results_df['Detected species (Bracken)'] = self.qc_data['bracken']['bracken_primary_species']
        results_df['Detected species (Mash)'] = self.qc_data['genome_size']['organism_name']

        # Change column order
        results_df = results_df[['sample', 'Detected species (Bracken)', 'Detected species (Mash)', 'Passed bracken', 'Passed mlst', 'Passed checkm', 'Passed assembly_scan', 'Passed fastp']]

        # Write the results to file
        results_df.to_csv(f"{self.sample_name}_qc_results.tsv", sep='\t')
        
        # Return pass or fail for each QC check
        return self.qc_results
    
    def get_qc_thresholds(self):
        """
        Returns the quality control results.
        """
        # Return requirements for each QC check
        return self.qc_requirements
        
        

# Define __all__ for explicit exports
__all__ = [
    'Genome',
]
