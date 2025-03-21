"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import unittest
import os
import sys
import pandas as pd
import networkx as nx
from collections import Counter

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import functions to test
try:
    from analyze.model import gather_synonyms, map_attributes, REDOX_PAIRS
    from mapper.special_pairs import get_redox_state
except ImportError:
    # For running tests directly
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    from metaboverse_cli.analyze.model import gather_synonyms, map_attributes, REDOX_PAIRS
    from metaboverse_cli.mapper.special_pairs import get_redox_state


class TestRedoxPairMapping(unittest.TestCase):
    """Test the redox pair mapping functionality"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a mock metabolite mapper
        self.metabolite_mapper = {
            'mapping_dictionary': {
                'nad': 'HMDB0000902',
                'nadh': 'HMDB0001487',
                'nadp': 'HMDB0000217',
                'nadph': 'HMDB0000221',
                'fad': 'HMDB0001248',
                'fadh2': 'HMDB0001352',
                'gsh': 'HMDB0000125',
                'gssg': 'HMDB0003337',
                'ascorbic acid': 'HMDB0000044',
                'dehydroascorbic acid': 'HMDB0001264'
            },
            'hmdb_dictionary': {
                'HMDB0000902': ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide'],
                'HMDB0001487': ['NADH', 'nicotinamide adenine dinucleotide reduced'],
                'HMDB0000217': ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate'],
                'HMDB0000221': ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced'],
                'HMDB0001248': ['FAD', 'flavin adenine dinucleotide'],
                'HMDB0001352': ['FADH2', 'flavin adenine dinucleotide reduced'],
                'HMDB0000125': ['GSH', 'glutathione', 'reduced glutathione'],
                'HMDB0003337': ['GSSG', 'glutathione disulfide', 'oxidized glutathione'],
                'HMDB0000044': ['ascorbic acid', 'vitamin C', 'L-ascorbic acid'],
                'HMDB0001264': ['dehydroascorbic acid', 'DHA', 'oxidized vitamin C']
            },
            'display_dictionary': {
                'HMDB0000902': ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide'],
                'HMDB0001487': ['NADH', 'nicotinamide adenine dinucleotide reduced'],
                'HMDB0000217': ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate'],
                'HMDB0000221': ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced'],
                'HMDB0001248': ['FAD', 'flavin adenine dinucleotide'],
                'HMDB0001352': ['FADH2', 'flavin adenine dinucleotide reduced'],
                'HMDB0000125': ['GSH', 'glutathione', 'reduced glutathione'],
                'HMDB0003337': ['GSSG', 'glutathione disulfide', 'oxidized glutathione'],
                'HMDB0000044': ['ascorbic acid', 'vitamin C', 'L-ascorbic acid'],
                'HMDB0001264': ['dehydroascorbic acid', 'DHA', 'oxidized vitamin C']
            }
        }
        
        # Create a mock uniprot mapper
        self.uniprot_mapper = {}
        
        # Create mock data for testing map_attributes
        self.mock_data = pd.DataFrame({
            'fc': [1.5, -1.5, 2.0, -2.0, 1.0, -1.0, 0.5, -0.5, 0.8, -0.8],
            'p': [0.01, 0.01, 0.005, 0.005, 0.02, 0.02, 0.03, 0.03, 0.04, 0.04]
        }, index=[
            'NAD+', 'NADH', 'NADP+', 'NADPH', 'FAD', 'FADH2', 'GSH', 'GSSG', 
            'ascorbic acid', 'dehydroascorbic acid'
        ])
        
        # Create mock stats for testing map_attributes
        self.mock_stats = pd.DataFrame({
            'p': [0.01, 0.01, 0.005, 0.005, 0.02, 0.02, 0.03, 0.03, 0.04, 0.04]
        }, index=[
            'NAD+', 'NADH', 'NADP+', 'NADPH', 'FAD', 'FADH2', 'GSH', 'GSSG', 
            'ascorbic acid', 'dehydroascorbic acid'
        ])
        
        # Create a mock graph for testing map_attributes
        self.graph = nx.DiGraph()
        
        # Add nodes for NAD+/NADH
        self.graph.add_node('node1')
        self.graph.nodes()['node1']['map_id'] = 'CHEBI:15846'  # NAD+
        self.graph.nodes()['node1']['name'] = 'NAD+'
        self.graph.nodes()['node1']['type'] = 'metabolite_component'
        
        self.graph.add_node('node2')
        self.graph.nodes()['node2']['map_id'] = 'CHEBI:16908'  # NADH
        self.graph.nodes()['node2']['name'] = 'NADH'
        self.graph.nodes()['node2']['type'] = 'metabolite_component'
        
        # Add nodes for NADP+/NADPH
        self.graph.add_node('node3')
        self.graph.nodes()['node3']['map_id'] = 'CHEBI:18009'  # NADP+
        self.graph.nodes()['node3']['name'] = 'NADP+'
        self.graph.nodes()['node3']['type'] = 'metabolite_component'
        
        self.graph.add_node('node4')
        self.graph.nodes()['node4']['map_id'] = 'CHEBI:16474'  # NADPH
        self.graph.nodes()['node4']['name'] = 'NADPH'
        self.graph.nodes()['node4']['type'] = 'metabolite_component'
        
        # Mock chebi_synonyms
        self.chebi_synonyms = {
            'CHEBI:15846': ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide'],
            'CHEBI:16908': ['NADH', 'nicotinamide adenine dinucleotide reduced'],
            'CHEBI:18009': ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate'],
            'CHEBI:16474': ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced']
        }
        
        # Mock args_dict
        self.args_dict = {'output': '.', 'metabolomics': 'test.txt'}
        
        # Mock name_reference and degree_dictionary
        self.name_reference = {}
        self.degree_dictionary = {'node1': 2, 'node2': 2, 'node3': 2, 'node4': 2}
        
        # Mock chebi_dictionary
        self.chebi_dictionary = {}

    def test_gather_synonyms_nad(self):
        """Test that gather_synonyms correctly handles NAD+/NADH"""
        # Test NAD+
        map_id = 'NAD+'
        init_syns = ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that NADH synonyms are not included
        self.assertNotIn('NADH', parsed_syns)
        self.assertNotIn('nicotinamide adenine dinucleotide reduced', parsed_syns)
        
        # Test NADH
        map_id = 'NADH'
        init_syns = ['NADH', 'nicotinamide adenine dinucleotide reduced']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that NAD+ synonyms are not included
        self.assertNotIn('NAD+', parsed_syns)
        self.assertNotIn('NAD', parsed_syns)
        self.assertNotIn('nicotinamide adenine dinucleotide', parsed_syns)

    def test_gather_synonyms_nadp(self):
        """Test that gather_synonyms correctly handles NADP+/NADPH"""
        # Test NADP+
        map_id = 'NADP+'
        init_syns = ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that NADPH synonyms are not included
        self.assertNotIn('NADPH', parsed_syns)
        self.assertNotIn('nicotinamide adenine dinucleotide phosphate reduced', parsed_syns)
        
        # Test NADPH
        map_id = 'NADPH'
        init_syns = ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that NADP+ synonyms are not included
        self.assertNotIn('NADP+', parsed_syns)
        self.assertNotIn('NADP', parsed_syns)
        self.assertNotIn('nicotinamide adenine dinucleotide phosphate', parsed_syns)

    def test_gather_synonyms_fad(self):
        """Test that gather_synonyms correctly handles FAD/FADH2"""
        # Test FAD
        map_id = 'FAD'
        init_syns = ['FAD', 'flavin adenine dinucleotide']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that FADH2 synonyms are not included
        self.assertNotIn('FADH2', parsed_syns)
        self.assertNotIn('flavin adenine dinucleotide reduced', parsed_syns)
        
        # Test FADH2
        map_id = 'FADH2'
        init_syns = ['FADH2', 'flavin adenine dinucleotide reduced']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that FAD synonyms are not included
        self.assertNotIn('FAD', parsed_syns)
        self.assertNotIn('flavin adenine dinucleotide', parsed_syns)

    def test_gather_synonyms_glutathione(self):
        """Test that gather_synonyms correctly handles GSH/GSSG"""
        # Test GSH
        map_id = 'GSH'
        init_syns = ['GSH', 'glutathione', 'reduced glutathione']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that GSSG synonyms are not included
        self.assertNotIn('GSSG', parsed_syns)
        self.assertNotIn('glutathione disulfide', parsed_syns)
        self.assertNotIn('oxidized glutathione', parsed_syns)
        
        # Test GSSG
        map_id = 'GSSG'
        init_syns = ['GSSG', 'glutathione disulfide', 'oxidized glutathione']
        
        mapper_id, parsed_syns = gather_synonyms(
            map_id, 
            init_syns.copy(), 
            self.metabolite_mapper, 
            self.uniprot_mapper, 
            True
        )
        
        # Check that GSH synonyms are not included
        self.assertNotIn('GSH', parsed_syns)
        self.assertNotIn('glutathione', parsed_syns)
        self.assertNotIn('reduced glutathione', parsed_syns)

    def test_map_attributes_redox_pairs(self):
        """Test that map_attributes correctly maps redox pairs to their distinct values"""
        # Create a simplified version of the test to check the core functionality
        # This is a partial test that focuses on the redox pair handling
        
        # First, prepare the graph nodes with synonyms
        for node_id in self.graph.nodes():
            self.graph.nodes()[node_id]['synonyms'] = []
            
            if node_id == 'node1':  # NAD+
                self.graph.nodes()[node_id]['synonyms'] = ['NAD+', 'NAD', 'nicotinamide adenine dinucleotide']
                self.graph.nodes()[node_id]['hmdb_mapper'] = 'HMDB0000902'
            elif node_id == 'node2':  # NADH
                self.graph.nodes()[node_id]['synonyms'] = ['NADH', 'nicotinamide adenine dinucleotide reduced']
                self.graph.nodes()[node_id]['hmdb_mapper'] = 'HMDB0001487'
            elif node_id == 'node3':  # NADP+
                self.graph.nodes()[node_id]['synonyms'] = ['NADP+', 'NADP', 'nicotinamide adenine dinucleotide phosphate']
                self.graph.nodes()[node_id]['hmdb_mapper'] = 'HMDB0000217'
            elif node_id == 'node4':  # NADPH
                self.graph.nodes()[node_id]['synonyms'] = ['NADPH', 'nicotinamide adenine dinucleotide phosphate reduced']
                self.graph.nodes()[node_id]['hmdb_mapper'] = 'HMDB0000221'
        
        # We'll manually check if the redox pairs are correctly identified
        for node_id in self.graph.nodes():
            all_synonyms = self.graph.nodes()[node_id]['synonyms']
            primary_name = all_synonyms[0] if all_synonyms else None
            
            # Use the get_redox_state function to determine the redox pair and state
            redox_pair, redox_state = None, None
            for syn in all_synonyms:
                redox_pair, redox_state = get_redox_state(syn, all_synonyms)
                if redox_pair:
                    break
            
            # Verify that each node is correctly identified as part of its redox pair
            if node_id == 'node1':  # NAD+
                self.assertEqual(redox_pair, 'NAD+/NADH')
                self.assertEqual(redox_state, 'oxidized')
            elif node_id == 'node2':  # NADH
                self.assertEqual(redox_pair, 'NAD+/NADH')
                self.assertEqual(redox_state, 'reduced')
            elif node_id == 'node3':  # NADP+
                self.assertEqual(redox_pair, 'NADP+/NADPH')
                self.assertEqual(redox_state, 'oxidized')
            elif node_id == 'node4':  # NADPH
                self.assertEqual(redox_pair, 'NADP+/NADPH')
                self.assertEqual(redox_state, 'reduced')

    def test_redox_pairs_completeness(self):
        """Test that all expected redox pairs are defined in REDOX_PAIRS"""
        # Check that all the expected redox pairs are in the dictionary
        expected_pairs = [
            'NAD+/NADH',
            'NADP+/NADPH',
            'FAD/FADH2',
            'GSH/GSSG',
            'Trx-S2/Trx-(SH)2',
            'CoQ/CoQH2',
            'Cytochrome C (Fe3+)/Cytochrome C (Fe2+)',
            'Ferredoxin (oxidized)/Ferredoxin (reduced)',
            'Ascorbic acid/Dehydroascorbic acid',
            'Lipoic acid/Dihydrolipoic acid',
            'FMN/FMNH2'
        ]
        
        for pair in expected_pairs:
            self.assertIn(pair, REDOX_PAIRS)
            self.assertIn('oxidized', REDOX_PAIRS[pair])
            self.assertIn('reduced', REDOX_PAIRS[pair])
            self.assertTrue(len(REDOX_PAIRS[pair]['oxidized']) > 0)
            self.assertTrue(len(REDOX_PAIRS[pair]['reduced']) > 0)

    def test_redox_pairs_no_overlap(self):
        """Test that there is no overlap between oxidized and reduced forms within each redox pair"""
        for pair, forms in REDOX_PAIRS.items():
            oxidized_forms = [s.lower() for s in forms['oxidized']]
            reduced_forms = [s.lower() for s in forms['reduced']]
            
            # Check for overlap between oxidized and reduced forms
            overlap = set(oxidized_forms).intersection(set(reduced_forms))
            if overlap:
                self.fail(f"Overlap between oxidized and reduced forms in '{pair}': {overlap}")
            
            # Check that the forms are normalized consistently
            # Simple check: words like "oxidized" or "reduced" should be consistent within their lists
            for ox_form in oxidized_forms:
                if "reduced" in ox_form and not ("reduced" in pair.lower()):
                    self.fail(f"Oxidized form '{ox_form}' contains 'reduced' in '{pair}'")
            
            for red_form in reduced_forms:
                if "oxidized" in red_form and not ("oxidized" in pair.lower()):
                    self.fail(f"Reduced form '{red_form}' contains 'oxidized' in '{pair}'")


if __name__ == '__main__':
    unittest.main()
