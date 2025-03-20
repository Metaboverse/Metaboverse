#!/usr/bin/env python3
"""
Test redox pair mapping with sample data
"""
import os
import sys
import pandas as pd
import networkx as nx
from collections import Counter

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import functions to test
try:
    from analyze.model import gather_synonyms, REDOX_PAIRS
except ImportError:
    # For running tests directly
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    from metaboverse_cli.analyze.model import gather_synonyms, REDOX_PAIRS

def test_with_sample_data():
    """Test redox pair mapping with sample data"""
    # Load sample data
    data_file = os.path.join(os.path.dirname(__file__), 'redox_test_data.txt')
    data = pd.read_csv(data_file, sep='\t', index_col=0)
    
    print(f"Loaded sample data with {len(data)} metabolites")
    print(data.head())
    
    # Create a mock metabolite mapper
    metabolite_mapper = {
        'mapping_dictionary': {},
        'hmdb_dictionary': {},
        'display_dictionary': {}
    }
    
    # Create a mock uniprot mapper
    uniprot_mapper = {}
    
    # Test each redox pair in the data
    redox_pairs_found = {}
    
    for metabolite in data.index:
        # Find which redox pair this metabolite belongs to
        for pair_name, synonyms in REDOX_PAIRS.items():
            if metabolite in synonyms or metabolite.lower() in [s.lower() for s in synonyms]:
                if pair_name not in redox_pairs_found:
                    redox_pairs_found[pair_name] = []
                redox_pairs_found[pair_name].append(metabolite)
                break
    
    # Print results
    print("\nRedox pairs found in the data:")
    for pair, metabolites in redox_pairs_found.items():
        print(f"  {pair}: {', '.join(metabolites)}")
    
    # Test gather_synonyms for each redox pair
    print("\nTesting gather_synonyms for each redox pair:")
    for pair, metabolites in redox_pairs_found.items():
        for metabolite in metabolites:
            init_syns = [metabolite]
            mapper_id, parsed_syns = gather_synonyms(
                metabolite, 
                init_syns.copy(), 
                metabolite_mapper, 
                uniprot_mapper, 
                True
            )
            
            # Check if any synonyms from other redox pairs are included
            other_pairs = [p for p in REDOX_PAIRS.keys() if p != pair]
            other_synonyms = []
            for other_pair in other_pairs:
                other_synonyms.extend(REDOX_PAIRS[other_pair])
            
            overlap = [s for s in parsed_syns if s in other_synonyms]
            
            if overlap:
                print(f"  ❌ {metabolite}: Found overlap with other redox pairs: {overlap}")
            else:
                print(f"  ✓ {metabolite}: No overlap with other redox pairs")
    
    # Check if redox pairs have opposite regulation
    print("\nChecking regulation patterns:")
    for pair_name, synonyms in REDOX_PAIRS.items():
        # Find the oxidized and reduced forms in the data
        oxidized = None
        reduced = None
        
        # The first item in the dictionary is typically the oxidized form
        # and the second item is typically the reduced form
        oxidized_name = list(REDOX_PAIRS.keys())[list(REDOX_PAIRS.keys()).index(pair_name)]
        reduced_name = None
        
        # Find the corresponding reduced form
        if oxidized_name.endswith('+'):
            # For NAD+/NADH and NADP+/NADPH
            reduced_candidate = oxidized_name.replace('+', 'H')
            if reduced_candidate in REDOX_PAIRS:
                reduced_name = reduced_candidate
        elif oxidized_name == 'FAD':
            reduced_name = 'FADH2'
        elif oxidized_name == 'GSH':
            reduced_name = 'GSSG'
        elif oxidized_name == 'Ascorbic acid':
            reduced_name = 'Dehydroascorbic acid'
        elif oxidized_name == 'Lipoic acid':
            reduced_name = 'Dihydrolipoic acid'
        elif oxidized_name == 'FMN':
            reduced_name = 'FMNH2'
        elif oxidized_name == 'CoQ':
            reduced_name = 'CoQH2'
        
        # Skip if we couldn't identify the pair
        if not reduced_name:
            continue
        
        # Find these metabolites in the data
        oxidized_in_data = [m for m in data.index if m in REDOX_PAIRS[oxidized_name] or 
                           m.lower() in [s.lower() for s in REDOX_PAIRS[oxidized_name]]]
        reduced_in_data = [m for m in data.index if m in REDOX_PAIRS[reduced_name] or 
                          m.lower() in [s.lower() for s in REDOX_PAIRS[reduced_name]]]
        
        if oxidized_in_data and reduced_in_data:
            oxidized_fc = data.loc[oxidized_in_data[0], 'log2fc']
            reduced_fc = data.loc[reduced_in_data[0], 'log2fc']
            
            if (oxidized_fc > 0 and reduced_fc < 0) or (oxidized_fc < 0 and reduced_fc > 0):
                print(f"  ✓ {oxidized_in_data[0]}/{reduced_in_data[0]}: Opposite regulation (FC: {oxidized_fc:.2f}/{reduced_fc:.2f})")
            else:
                print(f"  ❌ {oxidized_in_data[0]}/{reduced_in_data[0]}: Same regulation (FC: {oxidized_fc:.2f}/{reduced_fc:.2f})")

if __name__ == '__main__':
    test_with_sample_data() 