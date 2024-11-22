# bactQC/bactQC.py
import os
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import json

def get_expected_genome_size(sample_name, input_dir, taxid=None):
    """
    Retrieve expected genome size information from the NCBI API or infer it based on bracken abundance data if taxid is not provided.
    
    Parameters
    ----------
    sample_name : str
        The name of the sample.
    input_dir : str
        The directory path where input files are located.
    taxid : int or str, optional
        The taxonomy ID of the species. If not provided, the species is inferred from bracken data.
    
    Returns
    -------
    dict
        A dictionary containing expected genome size data, including:
            - 'sample': Sample name.
            - 'organism_name': Name of the organism.
            - 'species_taxid': Taxonomy ID of the species.
            - 'expected_ungapped_length': Expected genome size.
            - 'minimum_ungapped_length': Minimum expected genome size.
            - 'maximum_ungapped_length': Maximum expected genome size.
    """
    base_URL = "https://api.ncbi.nlm.nih.gov/genome/v0/expected_genome_size/expected_genome_size?species_taxid="
    data = {'sample': sample_name}
    
    if taxid is None:
        print(f"Expected species not provided, guessing species based on bracken data for {sample_name} - make sure you have generated it!")
        bracken_path = os.path.join(input_dir, sample_name, 'tools', 'bracken', f"{sample_name}.bracken.adjusted.abundances.txt")
        bracken_result = pd.read_csv(bracken_path, sep='\t')
        taxid = bracken_result.iloc[0]['taxonomy_id']
        guessed_species_name = bracken_result.iloc[0]['name']
        print(f"Getting expected, minimum and maximum genome size for {guessed_species_name} {taxid} from:")
        print(f"{base_URL}{taxid}")
    else:
        taxid = str(taxid)
    
    url = f"{base_URL}{taxid}"
    response = requests.get(url)
    
    if response.content:
        try:
            root = ET.fromstring(response.content)
            data.update({
                'organism_name': root.findtext('organism_name'),
                'species_taxid': root.findtext('species_taxid'),
                'expected_ungapped_length': int(root.findtext('expected_ungapped_length')),
                'minimum_ungapped_length': int(root.findtext('minimum_ungapped_length')),
                'maximum_ungapped_length': int(root.findtext('maximum_ungapped_length'))
            })
        except ET.ParseError:
            print(f"Error parsing XML for taxid {taxid}")
    else:
        print(f"No content returned for taxid {taxid}")
    
    return data

def get_assembly_size(sample_name, input_dir):
    """
    Retrieves the total contig length in base-pairs for a given sample from the assembler results.

    Parameters
    ----------
    sample_name : str
        The name of the sample.
    input_dir : str
        The directory path where input files are located.

    Returns
    -------
    dict
        A dictionary containing:
            - 'sample': Sample name.
            - 'total_length': Total length of contigs for the sample in base-pairs.
    """
    # Construct the file path to the assembler results
    file_path = os.path.join(input_dir, sample_name, 'main', 'assembler', f"{sample_name}.tsv")
    
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Assembly results file not found for sample '{sample_name}' at '{file_path}'")
    
    # Read the assembler results file
    assembly_scan_df = pd.read_csv(file_path, sep='\t')
    
    # Ensure the DataFrame has exactly one row
    if assembly_scan_df.shape[0] != 1:
        raise ValueError(f"Sample '{sample_name}': Expected one row in assembly scan data, but got {assembly_scan_df.shape[0]} rows")
    
    # Extract the total_contig_length from the first (and only) row
    total_length = assembly_scan_df.iloc[0]['total_contig_length']
    
    return {'sample': sample_name, 'total_length': total_length}

def check_bracken_results(sample_name, input_dir, min_primary_abundance=0.80):
    """
    Check Bracken results for a given sample.
    This function reads Bracken results from a TSV file, verifies the data, and determines
    whether the results pass quality control based on primary species abundance and genus consistency.
    Args:
        sample_name (str): The name of the sample.
        input_dir (str): The directory where the sample data is located.
        min_primary_abundance (float, optional): The minimum required abundance for the primary species. Defaults to 0.80.
    Returns:
        dict: A dictionary containing the Bracken results with additional QC flags, including whether the primary species abundance passes the threshold, the primary abundance requirement, and any genus conflicts.
    Raises:
        ValueError: If the Bracken results file does not contain exactly one row.
    """
    # Get the file path
    file_path = os.path.join(input_dir, sample_name, 'tools', 'bracken', f"{sample_name}.bracken.tsv")
    
    # Read the file
    bracken_result = pd.read_csv(file_path, sep='\t')
    
    # Check that bracken_result has only one row
    if bracken_result.shape[0] != 1:
        raise ValueError(f"Sample {sample_name}: Expected one row for bracken, but got {bracken_result.shape[0]} rows")
    
    # Add a column with a boolean for whether bracken_primary_species_abundance is greater than min_primary_abundance for first row
    bracken_result['passed bracken QC'] = bracken_result['bracken_primary_species_abundance'] > min_primary_abundance
     
    # Take the first row and convert it to a dictionary
    bracken_result = bracken_result.iloc[0].to_dict()
    
    # Add an entry with the minimum primary abundance
    bracken_result['primary_abundance_requirement'] = min_primary_abundance
    
    # Check if the first word within bracken_primary_species differs from the first word within the bracken_secondary_species
    only_one_dominant_primary_species = bracken_result['bracken_secondary_species'] != 'No secondary abundance > 1%'
    
    # Check if our genera are Escherichia and Shigella
    not_ecoli_and_shigella = {'Escherichia', 'Shigella'} == {bracken_result['bracken_primary_species'], bracken_result['bracken_secondary_species']}
    
    # Determine if there is a conflict in the genera detected
    if only_one_dominant_primary_species and not_ecoli_and_shigella:
        # Check if the first word within bracken_primary_species differs from the first word within the bracken_secondary_species
        bracken_result['Genus_conflict'] = bracken_result['bracken_primary_species'].split()[0] != bracken_result['bracken_secondary_species'].split()[0]
    else:
        bracken_result['Genus_conflict'] = False
    
    return bracken_result

def check_mlst_results(sample_name, input_dir, expected_genus):
    """
    Check MLST results for a given sample.
    This function reads the MLST (Multilocus Sequence Typing) result from a specified directory,
    verifies that the result matches the expected genus, and returns the processed MLST data.
    Parameters:
        sample_name (str): The name of the sample to check.
        input_dir (str): The directory where the input files are located.
        expected_genus (str): The expected genus to filter the scheme species map.
    Returns:
        dict: A dictionary containing the MLST result information, including a 'passed_mlst' key indicating whether the MLST check passed.
    Raises:
        ValueError: If the MLST data does not contain exactly one row.
    """
    # Get the file path
    file_path = os.path.join(input_dir, sample_name, 'tools', 'mlst', f"{sample_name}.tsv")
    
    # Read the file
    mlst_result = pd.read_csv(file_path, sep='\t', header=None)
    
    # Define the column names
    mlst_result.columns = ['sample', 'scheme', 'ST', 'allele1', 'allele2', 'allele3', 'allele4', 'allele5', 'allele6', 'allele7']
    
    # Remove '.fna.gz' or '.fna' from the sample column
    mlst_result['sample'] = mlst_result['sample'].str.replace('.fna.gz', '').str.replace('.fna', '')
    
    # Check that mlst_result has only one row
    if mlst_result.shape[0] != 1:
        raise ValueError(f"Sample {sample_name}: Expected one row in MLST data, but got {mlst_result.shape[0]} rows")
    
    # Download the scheme species map
    scheme_species_map = pd.read_csv('https://raw.githubusercontent.com/tseemann/mlst/refs/heads/master/db/scheme_species_map.tab', sep='\t')
    
    # Filter the dataframe to only the expected genus
    scheme_species_map = scheme_species_map[scheme_species_map['GENUS'] == expected_genus]
    
    # Check if the scheme in mlst result partially matches with one or more schemes for the expected genus
    if not any(mlst_result['scheme'].str.contains('|'.join(scheme_species_map['#SCHEME']))):
        # Set mlst_result['passed_mlst'] to False
        mlst_result['passed_mlst'] = False
    else:
        # Set mlst_result['passed_mlst'] to True
        mlst_result['passed_mlst'] = True 
    
    # Take the first row and convert it to a dictionary
    mlst_result = mlst_result.iloc[0].to_dict()
      
    return mlst_result

def check_checkm_results(sample_name, input_dir, min_completeness=0.80, max_contamination=10):
    """
    Checks CheckM results for a given sample and evaluates quality metrics.
    This function searches recursively for 'checkm.tsv' files within the specified
    input directory, selects the most recent file based on directory naming,
    and verifies the completeness and contamination levels for the specified sample.
    Parameters:
        sample_name (str): The name of the sample to check in the CheckM results.
        input_dir (str): The directory path to search for 'checkm.tsv' files.
        min_completeness (float, optional): The minimum required completeness threshold.
            Defaults to 0.80.
        max_contamination (float, optional): The maximum allowed contamination threshold.
            Defaults to 10.
    Returns:
        dict: A dictionary containing the CheckM results for the sample, including
            completeness and contamination values, the specified requirements, and
            boolean flags indicating whether the sample passes the completeness and
            contamination criteria.
    Raises:
        ValueError: If no 'checkm.tsv' files are found in the input directory.
        ValueError: If the specified sample name is not found in the CheckM results.
    """
    # Search recursively for all files in the input directory called 'checkm.tsv'
    checkm_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file == 'checkm.tsv':
                checkm_files.append(os.path.join(root, file))
                
    # Check that we found at least one checkm.tsv file
    if len(checkm_files) == 0:
        raise ValueError("No checkm.tsv files found")
    
    # Determine the most recent checkm.tsv file based on the timestamp name of directory two levels up
    checkm_files = sorted(checkm_files, key=lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))))   
    
    # Select the most recent checkm.tsv file
    file_path = checkm_files[-1]
  
    # Read the file
    checkm_result = pd.read_csv(file_path, sep='\t')
    
    # Rename 'Bin Id' to 'sample'
    checkm_result.rename(columns={'Bin Id': 'sample'}, inplace=True)
    
    # Check if sample name exists in the checkm file
    if sample_name not in checkm_result['sample'].values:
        raise ValueError(f"Sample {sample_name} not found in {file_path}")
    
    # Filter the dataframe to only the sample of interest
    checkm_result = checkm_result[checkm_result['sample'] == sample_name]
    
    # Create a dictionary from the first row
    checkm_result = checkm_result.iloc[0].to_dict()
    
    # Add an entry with the minimum completeness
    checkm_result['completeness_requirement'] = min_completeness
    
    # Add an entry with the maximum contamination
    checkm_result['contamination_requirement'] = max_contamination
    
    # Add an entry as a boolean for if the completeness is greater than min_completeness
    checkm_result['passed_completeness'] = checkm_result['Completeness'] > min_completeness
    
    # Add an entry as a boolean for if the contamination is less than max_contamination
    checkm_result['passed_contamination'] = checkm_result['Contamination'] < max_contamination
    
    # Add an entry as a boolean for if completion and contamination requirements are met
    checkm_result['passed_checkm QC'] = checkm_result['passed_completeness'] and checkm_result['passed_contamination']
    
    return checkm_result

def check_assembly_scan_results(sample_name, input_dir, maximum_contigs=500, minimum_N50=15000, taxid=None):
    """
    Check the quality of assembly scan results for a given sample.
    This function verifies assembly metrics such as the number of contigs, N50 value,
    and genome size against specified thresholds. It retrieves expected genome size
    information from the NCBI API or infers it based on bracken abundance data if taxid is not provided.
    
    Parameters
    ----------
    sample_name : str
        The name of the sample to be checked.
    input_dir : str
        The directory path where input files are located.
    maximum_contigs : int, optional
        The maximum allowed number of contigs (default is 500).
    minimum_N50 : int, optional
        The minimum required N50 contig length (default is 15000).
    taxid : int or str, optional
        The taxonomy ID of the species. If not provided, the species is inferred from bracken data.
    
    Returns
    -------
    dict
        A dictionary containing assembly scan results and QC pass/fail flags, including:
            - 'sample': Sample name.
            - 'organism_name': Name of the organism.
            - 'species_taxid': Taxonomy ID of the species.
            - 'expected_ungapped_length': Expected genome size.
            - 'minimum_ungapped_length': Minimum expected genome size.
            - 'maximum_ungapped_length': Maximum expected genome size.
            - 'passed_contigs': Boolean indicating if contig count is below maximum.
            - 'passed_N50': Boolean indicating if N50 is above minimum.
            - 'passed_genome_size': Boolean indicating if genome size is within acceptable range.
            - 'passed_assembly_scan': Boolean indicating if all QC checks are passed.
    """
    # Retrieve expected genome size data
    data = get_expected_genome_size(sample_name, input_dir, taxid)
    
    # Get the file path
    file_path = os.path.join(input_dir, sample_name, 'main', 'assembler', f"{sample_name}.tsv")
    
    # Read the file
    assembly_scan_df = pd.read_csv(file_path, sep='\t')
    
    if assembly_scan_df.shape[0] != 1:
        raise ValueError(f"Sample {sample_name}: Expected one row in assembly scan data, but got {assembly_scan_df.shape[0]} rows")
    
    # Take the first row and convert it to a dictionary
    assembly_scan_results = assembly_scan_df.iloc[0].to_dict()
    
    # Add QC pass/fail flags
    assembly_scan_results['passed_contigs'] = assembly_scan_results['total_contig'] < maximum_contigs
    assembly_scan_results['passed_N50'] = assembly_scan_results['n50_contig_length'] > minimum_N50
    assembly_scan_results.update(data)
    
    # Calculate acceptable genome size range
    min_length = data['minimum_ungapped_length'] * 0.8
    max_length = data['maximum_ungapped_length'] * 1.2
    total_length = assembly_scan_results['total_contig_length']
    
    assembly_scan_results['passed_genome_size'] = (total_length > min_length) and (total_length < max_length)
    assembly_scan_results['passed_assembly_scan'] = all([
        assembly_scan_results['passed_contigs'],
        assembly_scan_results['passed_N50'],
        assembly_scan_results['passed_genome_size']
    ])
    
    return assembly_scan_results

def check_fastp_data(sample_name, input_dir, min_q30_bases=0.80, min_coverage=30):
    """
    Checks fastp quality control data for a given sample.
    This function reads a fastp JSON summary file for the specified sample,
    extracts relevant metrics before and after filtering, calculates coverage
    based on the assembly size, and determines whether the sample passes
    predefined quality thresholds.
    Parameters:
        sample_name (str): The name of the sample.
        input_dir (str): The directory containing sample data.
        min_q30_bases (float, optional): Minimum required proportion of Q30 bases after filtering. Defaults to 0.80.
        min_coverage (int, optional): Minimum required coverage after filtering. Defaults to 30.
    Returns:
        dict: A dictionary containing extracted metrics and boolean flags indicating
              whether the sample passed Q30 base rate and coverage thresholds.
    """
    
    # Get the file path
    file_path = os.path.join(input_dir, sample_name, 'main', 'qc', 'summary', f"{sample_name}.fastp.json")
    
    # Read the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)
        
        # Create a dictionary with the name as a single-element list and set it as the index
        fastp_results = {'sample': sample_name}
    
        fastp_results['pre_filt_total_reads'] = int(data['summary']['before_filtering']['total_reads'])
        fastp_results['pre_filt_total_bases'] = int(data['summary']['before_filtering']['total_bases'])
        fastp_results['pre_filt_q30_rate'] = float(data['summary']['before_filtering']['q30_rate'])
        fastp_results['pre_filt_gc'] = float(data['summary']['before_filtering']['gc_content'])
        fastp_results['post_filt_total_reads'] = int(data['summary']['after_filtering']['total_reads'])
        fastp_results['post_filt_total_bases'] = int(data['summary']['after_filtering']['total_bases'])
        fastp_results['post_filt_q20_rate'] = float(data['summary']['after_filtering']['q20_rate'])
        fastp_results['post_filt_q30_rate'] = float(data['summary']['after_filtering']['q30_rate'])
        fastp_results['post_filt_gc'] = float(data['summary']['after_filtering']['gc_content'])
    
    # Get our assembly size
    assembly_size = get_assembly_size(sample_name, input_dir)
    
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
   
    return fastp_results

def run(sample_name, input_dir, taxid=None):
    """
    Run all quality control checks for a given sample.
    This function runs all quality control checks for a specified sample
    and returns the results as a dictionary.
    Parameters:
        sample_name (str): The name of the sample to check.
        input_dir (str): The directory path where input files are located.
        taxid (int or str, optional): The taxonomy ID of the species. Defaults to None.
    Returns:
        dict: A dictionary containing the results of all quality control checks.
    """
    qc_results = {}
    
    # Run all quality control checks
    qc_results['genome_size'] = get_expected_genome_size(sample_name, input_dir, taxid)
    qc_results['assembly_size'] = get_assembly_size(sample_name, input_dir)
    qc_results['bracken'] = check_bracken_results(sample_name, input_dir)
    qc_results['mlst'] = check_mlst_results(sample_name, input_dir, qc_results['genome_size']['organism_name'].split()[0])
    qc_results['checkm'] = check_checkm_results(sample_name, input_dir)
    qc_results['assembly_scan'] = check_assembly_scan_results(sample_name, input_dir, taxid=taxid)
    qc_results['fastp'] = check_fastp_data(sample_name, input_dir)
    
    return qc_results

# Define __all__ for explicit exports
__all__ = [
    'get_expected_genome_size',
    'get_assembly_size',
    'check_bracken_results',
    'check_mlst_results',
    'check_checkm_results',
    'check_assembly_scan_results',
    'check_fastp_data',
    'run',
]
