import unittest
from unittest.mock import patch, MagicMock, mock_open
import pandas as pd
import os
import json
from bactQC.core import Genome

class TestGenome(unittest.TestCase):
    def setUp(self):
        # Initialize a Genome instance for testing
        self.sample_name = 'sample1'
        self.input_dir = '/path/to/input'
        self.genome = Genome(self.sample_name, self.input_dir)

    def test_initialization(self):
        # Test that the Genome object initializes correctly
        self.assertEqual(self.genome.sample_name, self.sample_name)
        self.assertEqual(self.genome.input_dir, self.input_dir)
        self.assertIsNone(self.genome.taxid)
        self.assertEqual(self.genome.qc_data, {'sample': self.sample_name})
        self.assertEqual(self.genome.qc_results, {'sample': self.sample_name})
        self.assertEqual(self.genome.qc_requirements, {'sample': self.sample_name})

    @patch('bactQC.core.requests.get')
    def test_get_expected_genome_size(self, mock_get):
        # Mock response content
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

        # Set taxid and call the method
        self.genome.taxid = '562'
        self.genome.get_expected_genome_size()

        # Assert expected results
        expected_data = {
            'organism_name': 'Escherichia coli',
            'species_taxid': '562',
            'expected_ungapped_length': 5000000,
            'minimum_ungapped_length': 4500000,
            'maximum_ungapped_length': 5500000
        }
        self.assertEqual(self.genome.qc_data['genome_size'], expected_data)

    @patch('os.path.isfile')
    @patch('pandas.read_csv')
    def test_get_assembly_size(self, mock_read_csv, mock_isfile):
        # Mock file existence
        mock_isfile.return_value = True
        # Mock DataFrame returned by pandas.read_csv
        mock_df = pd.DataFrame({
            'sample': [self.sample_name],
            'total_contig_length': [5000000]
        })
        mock_read_csv.return_value = mock_df

        # Call the method
        self.genome.get_assembly_size()

        # Assert expected results
        expected_data = {'total_length': 5000000}
        self.assertEqual(self.genome.qc_data['assembly_size'], expected_data)

    @patch('os.path.isfile')
    @patch('pandas.read_csv')
    def test_check_bracken(self, mock_read_csv, mock_isfile):
        # Mock file existence
        mock_isfile.return_value = True
        # Mock DataFrame for Bracken results
        mock_df = pd.DataFrame({
            'sample': [self.sample_name],
            'bracken_primary_species_abundance': [0.85],
            'bracken_primary_species': ['Escherichia coli'],  # Added this line
            'taxonomy_id': [562],
            'name': ['Escherichia coli'],
            'bracken_secondary_species': ['No secondary abundance > 1%'],
            'bracken_secondary_species_abundance': [0.0]
        })
        mock_read_csv.return_value = mock_df

        # Call the method
        self.genome.check_bracken(min_primary_abundance=0.80)

        # Assert expected results
        bracken_result = self.genome.qc_data['bracken']
        self.assertTrue(bracken_result['passed_bracken_QC'])
        self.assertEqual(bracken_result['primary_abundance_requirement'], 0.80)
        self.assertFalse(bracken_result['genus_conflict'])

    @patch('os.path.exists')
    @patch('os.path.isfile')
    @patch('pandas.read_csv')
    def test_check_mlst(self, mock_read_csv, mock_isfile, mock_exists):
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

        # Call the method
        self.genome.check_mlst(expected_genus='Escherichia')

        # Assert expected results
        mlst_result = self.genome.qc_data['mlst']
        self.assertTrue(mlst_result['passed_mlst'])



    def test_get_expected_genome_size_no_taxid(self):
        # Remove taxid and mock Bracken file existence
        self.genome.taxid = None
        with patch('os.path.isfile') as mock_isfile:
            mock_isfile.return_value = False
            # Expect a FileNotFoundError
            with self.assertRaises(FileNotFoundError):
                self.genome.get_expected_genome_size()

    @patch('os.walk')
    @patch('pandas.read_csv')
    def test_check_checkm(self, mock_read_csv, mock_walk):
        # Mock os.walk to find checkm.tsv files
        mock_walk.return_value = [
            ('/path/to/input', ['dir1'], ['file1']),
            ('/path/to/input/dir1', [], ['checkm.tsv'])
        ]

        # Mock DataFrame for CheckM results
        mock_df = pd.DataFrame({
            'Bin Id': [self.sample_name],
            'Completeness': [95.0],
            'Contamination': [2.0],
            'Sample': [self.sample_name]
        })
        mock_read_csv.return_value = mock_df

        # Call the method
        self.genome.check_checkm(min_completeness=80, max_contamination=10)

        # Assert expected results
        checkm_result = self.genome.qc_data['checkm']
        self.assertTrue(checkm_result['passed_checkm_QC'])
        self.assertTrue(checkm_result['passed_completeness'])
        self.assertTrue(checkm_result['passed_contamination'])

    @patch('os.path.isfile')
    @patch('pandas.read_csv')
    def test_check_assembly_scan(self, mock_read_csv, mock_isfile):
        # Mock file existence
        mock_isfile.return_value = True

        # Mock DataFrame for assembly scan results
        mock_df = pd.DataFrame({
            'sample': [self.sample_name],
            'total_contig_length': [5000000],
            'total_contig': [100],
            'n50_contig_length': [20000]
        })
        mock_read_csv.return_value = mock_df

        # Provide genome size data
        self.genome.qc_data['genome_size'] = {
            'organism_name': 'Escherichia coli',
            'species_taxid': '562',
            'expected_ungapped_length': 5000000,
            'minimum_ungapped_length': 4500000,
            'maximum_ungapped_length': 5500000
        }

        # Call the method
        self.genome.check_assembly_scan(maximum_contigs=500, minimum_n50=15000)

        # Assert expected results
        assembly_scan_result = self.genome.qc_data['assembly_scan']
        self.assertTrue(assembly_scan_result['passed_assembly_scan'])
        self.assertTrue(assembly_scan_result['passed_contigs'])
        self.assertTrue(assembly_scan_result['passed_N50'])
        self.assertTrue(assembly_scan_result['passed_genome_size'])

    @patch('builtins.open', new_callable=mock_open, read_data='{"summary": {"before_filtering": {"total_reads": 1000000, "total_bases": 160000000, "q30_rate": 0.95, "gc_content": 0.5}, "after_filtering": {"total_reads": 950000, "total_bases": 155000000, "q20_rate": 0.98, "q30_rate": 0.96, "gc_content": 0.5}}}')
    @patch('os.path.isfile')
    def test_check_fastp(self, mock_isfile, mock_file):
        # Mock file existence
        mock_isfile.return_value = True

        # Provide assembly size data
        self.genome.qc_data['assembly_size'] = {'total_length': 5000000}

        # Call the method
        self.genome.check_fastp(min_q30_bases=0.90, min_coverage=30)

        # Assert expected results
        fastp_result = self.genome.qc_data['fastp']
        self.assertTrue(fastp_result['passed_fastp_QC'])
        self.assertTrue(fastp_result['passed_q30_bases'])
        self.assertTrue(fastp_result['passed_coverage'])
        self.assertAlmostEqual(fastp_result['coverage'], 31.0, delta=0.1)


    @patch('pandas.DataFrame.to_csv')
    def test_get_qc_results(self, mock_to_csv):
        # Provide qc_results data
        self.genome.qc_results = {
            'sample': self.sample_name,
            'bracken': True,
            'mlst': True,
            'checkm': True,
            'assembly_scan': True,
            'fastp': True
        }

        # Provide qc_data with detected species
        self.genome.qc_data = {
            'bracken': {
                'bracken_primary_species': 'Escherichia coli'
            },
            'genome_size': {
                'organism_name': 'Escherichia coli'
            }
        }

        # Call the method
        results = self.genome.get_qc_results()

        # Assert expected results
        self.assertEqual(results, self.genome.qc_results)
        # Check that to_csv was called
        mock_to_csv.assert_called_once()

    def test_get_qc_thresholds(self):
        # Provide qc_requirements data
        self.genome.qc_requirements = {
            'sample': self.sample_name,
            'bracken': {'min_primary_abundance': 0.80},
            'mlst': {'expected_genus': 'Escherichia'},
            'checkm': {'max_contamination': 10, 'min_completeness': 80},
            'assembly_scan': {'maximum_contigs': 500, 'minimum_n50': 15000, 'minimum_ungapped_length': 4500000, 'maximum_ungapped_length': 5500000},
            'fastp': {'min_q30_bases': 0.90, 'min_coverage': 30}
        }

        # Call the method
        thresholds = self.genome.get_qc_thresholds()

        # Assert expected results
        self.assertEqual(thresholds, self.genome.qc_requirements)

if __name__ == '__main__':
    unittest.main()
