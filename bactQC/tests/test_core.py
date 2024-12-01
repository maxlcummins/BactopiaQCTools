# test_core.py

import pytest
from unittest.mock import patch, MagicMock, mock_open
import pandas as pd
import os
import json
from bactQC.core import Genome

@pytest.fixture
def genome():
    sample_name = 'sample1'
    input_dir = '/path/to/input'
    return Genome(sample_name, input_dir)

def test_initialization(genome):
    # Test that the Genome object initializes correctly
    assert genome.sample_names == ['sample1']
    assert genome.input_dir == '/path/to/input'
    assert genome.qc_data == {}
    assert genome.qc_results == {}
    assert genome.qc_requirements == {}

@patch('os.path.isfile')
@patch('pandas.read_csv')
@patch('bactQC.core.requests.get')
def test_get_expected_genome_size(mock_get, mock_read_csv, mock_isfile, genome):
    # Mock Bracken file existence
    mock_isfile.return_value = True
    # Mock Bracken DataFrame
    bracken_df = pd.DataFrame({
        'sample': ['sample1'],
        'taxonomy_id': [562],
        'name': ['Escherichia coli'],
        'fraction_total_reads': [0.85],
        'kraken_assigned_reads': [850],
        'added_reads': [100],
        'new_est_reads': [950],
        'fraction_total_reads2': [0.95]
    })
    mock_read_csv.return_value = bracken_df

    # Mock response content for NCBI API
    mock_response = MagicMock()
    mock_response.content = b'''
        <root>
            <organism_name>Escherichia coli</organism_name>
            <species_taxid>562</species_taxid>
            <expected_ungapped_length>5000000</expected_ungapped_length>
            <minimum_ungapped_length>4500000</minimum_ungapped_length>
            <maximum_ungapped_length>5500000</maximum_ungapped_length>
        </root>
    '''
    mock_get.return_value = mock_response

    # Call the method with sample_name
    genome.get_expected_genome_size('sample1')

    # Assert expected results
    expected_data = {
        'organism_name': 'Escherichia coli',
        'species_taxid': '562',
        'expected_ungapped_length': 5000000,
        'minimum_ungapped_length': 4500000,
        'maximum_ungapped_length': 5500000
    }
    assert genome.qc_data['sample1']['genome_size'] == expected_data

@patch('os.path.isfile')
@patch('pandas.read_csv')
def test_get_assembly_size(mock_read_csv, mock_isfile, genome):
    # Mock file existence
    mock_isfile.return_value = True
    # Mock DataFrame returned by pandas.read_csv
    mock_df = pd.DataFrame({
        'sample': ['sample1'],
        'total_contig_length': [5000000]
    })
    mock_read_csv.return_value = mock_df

    # Call the method with sample_name
    genome.get_assembly_size('sample1')

    # Assert expected results
    expected_data = {'total_length': 5000000}
    assert genome.qc_data['sample1']['assembly_size'] == expected_data

@patch('os.path.isfile')
@patch('pandas.read_csv')
def test_check_bracken(mock_read_csv, mock_isfile, genome):
    # Mock file existence
    mock_isfile.return_value = True
    # Mock DataFrame for Bracken results
    mock_df = pd.DataFrame({
        'sample': ['sample1'],
        'bracken_primary_species_abundance': [0.85],
        'bracken_primary_species': ['Escherichia coli'],
        'taxonomy_id': [562],
        'name': ['Escherichia coli'],
        'bracken_secondary_species': ['No secondary abundance > 1%'],
        'bracken_secondary_species_abundance': [0.0]
    })
    mock_read_csv.return_value = mock_df

    # Call the method with sample_name
    genome.check_bracken('sample1', min_primary_abundance=0.80)

    # Assert expected results
    bracken_result = genome.qc_data['sample1']['bracken']
    assert bracken_result['passed_bracken_QC']
    assert bracken_result['primary_abundance_requirement'] == 0.80
    assert not bracken_result['genus_conflict']

@patch('os.path.exists')
@patch('os.path.isfile')
@patch('pandas.read_csv')
def test_check_mlst(mock_read_csv, mock_isfile, mock_exists, genome):
    # Mock file existence
    mock_isfile.return_value = True
    mock_exists.return_value = True  # Simulate that cache_file exists

    # Mock DataFrame for MLST results
    mlst_df = pd.DataFrame([
        ['sample1.fasta', 'ecoli', 10, 'adk(1)', 'fumC(1)', 'gyrB(1)', 'icd(1)', 'mdh(1)', 'purA(1)', 'recA(1)']
    ])

    # Mock DataFrame for scheme_species_map.tab
    scheme_df = pd.DataFrame({
        '#SCHEME': ['ecoli'],
        'GENUS': ['Escherichia']
    })

    def read_csv_side_effect(*args, **kwargs):
        if 'scheme_species_map.tab' in args[0]:
            return scheme_df
        else:
            return mlst_df

    mock_read_csv.side_effect = read_csv_side_effect

    # Provide genome_size data to get expected_genus
    genome.qc_data['sample1'] = {
        'genome_size': {
            'organism_name': 'Escherichia coli'
        }
    }
    expected_genus = 'Escherichia'

    # Call the method with sample_name and expected_genus
    genome.check_mlst('sample1', expected_genus)

    # Assert expected results
    mlst_result = genome.qc_data['sample1']['mlst']
    assert mlst_result['passed_mlst']

def test_get_expected_genome_size_no_bracken(genome):
    # Mock Bracken file existence to be False
    with patch('os.path.isfile') as mock_isfile:
        mock_isfile.return_value = False
        # Expect a FileNotFoundError
        with pytest.raises(FileNotFoundError):
            genome.get_expected_genome_size('sample1')

@patch('os.walk')
@patch('pandas.read_csv')
def test_check_checkm(mock_read_csv, mock_walk, genome):
    # Mock os.walk to find checkm.tsv files
    mock_walk.return_value = [
        ('/path/to/input', ['dir1'], ['file1']),
        ('/path/to/input/dir1', [], ['checkm.tsv'])
    ]

    # Mock DataFrame for CheckM results
    mock_df = pd.DataFrame({
        'Bin Id': ['sample1'],
        'Completeness': [95.0],
        'Contamination': [2.0]
    })
    mock_read_csv.return_value = mock_df

    # Call the method with sample_name
    genome.check_checkm('sample1', min_completeness=80, max_contamination=10)

    # Assert expected results
    checkm_result = genome.qc_data['sample1']['checkm']
    assert checkm_result['passed_checkm_QC']
    assert checkm_result['passed_completeness']
    assert checkm_result['passed_contamination']

@patch('os.path.isfile')
@patch('pandas.read_csv')
def test_check_assembly_scan(mock_read_csv, mock_isfile, genome):
    # Mock file existence
    mock_isfile.return_value = True

    # Mock DataFrame for assembly scan results
    mock_df = pd.DataFrame({
        'sample': ['sample1'],
        'total_contig_length': [5000000],
        'total_contig': [100],
        'n50_contig_length': [20000]
    })
    mock_read_csv.return_value = mock_df

    # Provide genome size data
    genome.qc_data['sample1'] = {
        'genome_size': {
            'organism_name': 'Escherichia coli',
            'species_taxid': '562',
            'expected_ungapped_length': 5000000,
            'minimum_ungapped_length': 4500000,
            'maximum_ungapped_length': 5500000
        }
    }

    # Call the method with sample_name
    genome.check_assembly_scan('sample1', maximum_contigs=500, minimum_n50=15000)

    # Assert expected results
    assembly_scan_result = genome.qc_data['sample1']['assembly_scan']
    assert assembly_scan_result['passed_assembly_scan']
    assert assembly_scan_result['passed_contigs']
    assert assembly_scan_result['passed_N50']
    assert assembly_scan_result['passed_genome_size']

@patch('os.path.isfile')
@patch('builtins.open', new_callable=mock_open, read_data='{"summary": {"before_filtering": {"total_reads": 1000000, "total_bases": 160000000, "q30_rate": 0.95, "gc_content": 0.5}, "after_filtering": {"total_reads": 950000, "total_bases": 155000000, "q20_rate": 0.98, "q30_rate": 0.96, "gc_content": 0.5}}}')
def test_check_fastp(mock_file, mock_isfile, genome):
    # Mock file existence
    mock_isfile.return_value = True

    # Provide assembly size data
    genome.qc_data['sample1'] = {
        'assembly_size': {'total_length': 5000000}
    }

    # Call the method with sample_name
    genome.check_fastp('sample1', min_q30_bases=0.90, min_coverage=30)

    # Assert expected results
    fastp_result = genome.qc_data['sample1']['fastp']
    assert fastp_result['passed_fastp_QC']
    assert fastp_result['passed_q30_bases']
    assert fastp_result['passed_coverage']
    from pytest import approx
    assert fastp_result['coverage'] == approx(31.0, abs=0.1)

@patch('pandas.DataFrame.to_csv')
def test_get_qc_results(mock_to_csv, genome):
    # Provide qc_results data
    genome.qc_results['sample1'] = {
        'bracken': True,
        'mlst': True,
        'checkm': True,
        'assembly_scan': True,
        'fastp': True,
        'overall': True
    }

    # Provide qc_data with detected species
    genome.qc_data['sample1'] = {
        'bracken': {
            'bracken_primary_species': 'Escherichia coli'
        },
        'genome_size': {
            'organism_name': 'Escherichia coli'
        }
    }

    # Call the method
    results = genome.get_qc_results()

    # Assert expected results
    expected_results = pd.DataFrame([{
        'sample': 'sample1',
        'Detected species (Bracken)': 'Escherichia coli',
        'Detected species (Mash)': 'Escherichia coli',
        'bracken': True,
        'mlst': True,
        'checkm': True,
        'assembly_scan': True,
        'fastp': True,
        'overall': True
    }])

    # Since results is a DataFrame, compare it to expected_results
    pd.testing.assert_frame_equal(results.reset_index(drop=True), expected_results, check_like=True)

    # Check that to_csv was called
    mock_to_csv.assert_called_once()

@patch('pandas.DataFrame.to_csv')
def test_get_qc_thresholds(mock_to_csv, genome):
    # Provide qc_requirements data
    genome.qc_requirements['sample1'] = {
        'bracken': {'min_primary_abundance': 0.80},
        'mlst': {'expected_genus': 'Escherichia'},
        'checkm': {'max_contamination': 10, 'min_completeness': 80},
        'assembly_scan': {'maximum_contigs': 500, 'minimum_n50': 15000, 'minimum_ungapped_length': 4500000, 'maximum_ungapped_length': 5500000},
        'fastp': {'min_q30_bases': 0.90, 'min_coverage': 30}
    }

    # Call the method
    thresholds = genome.get_qc_thresholds()

    # Assert expected results
    expected_thresholds = pd.DataFrame([{
        'sample': 'sample1',
        'bracken_min_primary_abundance': 0.80,
        'mlst_expected_genus': 'Escherichia',
        'checkm_max_contamination': 10,
        'checkm_min_completeness': 80,
        'assembly_scan_maximum_contigs': 500,
        'assembly_scan_minimum_n50': 15000,
        'assembly_scan_minimum_ungapped_length': 4500000,
        'assembly_scan_maximum_ungapped_length': 5500000,
        'fastp_min_q30_bases': 0.90,
        'fastp_min_coverage': 30
    }])

    # Since thresholds is a DataFrame, compare it to expected_thresholds
    pd.testing.assert_frame_equal(thresholds.reset_index(drop=True), expected_thresholds, check_like=True)

    # Check that to_csv was called
    mock_to_csv.assert_called_once()
