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
from __future__ import print_function
import pandas as pd
import xml.etree.ElementTree as et
import requests
import zipfile
import pickle
import io
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

def get_project_root():
    """Get the path to the project root directory"""
    current_file = Path(__file__).resolve()
    for parent in current_file.parents:
        if parent.name == 'cli':
            return parent
        if parent.name == 'Metaboverse':
            return parent / 'cli'
    return current_file.parent.parent.parent

# Add project root to path
project_root = get_project_root()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
if str(project_root.parent) not in sys.path:
    sys.path.insert(0, str(project_root.parent))

try:
    # First try normal package imports
    from metaboverse_cli.utils import prepare_output, write_database, write_database_json
    from metaboverse_cli.mapper.special_pairs import REDOX_PAIRS
except ImportError:
    try:
        # Then try relative imports
        from utils import prepare_output, write_database, write_database_json
        from special_pairs import REDOX_PAIRS
    except ImportError:
        try:
            # Finally try direct imports
            import importlib.util
            def load_module(name, path):
                spec = importlib.util.spec_from_file_location(name, path)
                if spec is None:
                    raise ImportError(f"Could not find module at {path}")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                return module

            base_path = os.path.dirname(__file__)  # Get the directory containing this file
            parent_path = os.path.dirname(os.path.dirname(base_path))  # Go up two levels to reach cli directory

            utils = load_module("utils", os.path.join(parent_path, "metaboverse_cli", "utils.py"))
            prepare_output = utils.prepare_output
            write_database = utils.write_database
            write_database_json = utils.write_database_json

            special_pairs = load_module("special_pairs", os.path.join(base_path, "special_pairs.py"))
            REDOX_PAIRS = special_pairs.REDOX_PAIRS

        except ImportError as e:
            print(f"Error importing dependencies: {e}")
            print(f"Current sys.path: {sys.path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"File location: {__file__}")
            print(f"Looking for special_pairs.py in: {os.path.join(base_path, 'special_pairs.py')}")
            raise


def update_redox_pairs(hmdb_dictionary, redox_pairs):
    """Check HMDB data for potential new redox pair synonyms
    
    Args:
        hmdb_dictionary (dict): Dictionary of HMDB metabolite data
        redox_pairs (dict): Current redox pairs dictionary
        
    Returns:
        dict: Dictionary of potential new synonyms for each redox pair
    """
    potential_additions = {}
    
    # Create a mapping of current synonyms to their redox pairs and states
    current_synonyms = {}
    for pair_name, forms in redox_pairs.items():
        for state in ['oxidized', 'reduced']:
            for syn in forms[state]:
                current_synonyms[syn.lower()] = (pair_name, state)
    
    # Check each HMDB entry for potential matches
    for hmdb_id, synonyms in hmdb_dictionary.items():
        for syn in synonyms:
            syn_lower = syn.lower()
            
            # Skip if we already know about this synonym
            if syn_lower in current_synonyms:
                continue
                
            # Check if this synonym might belong to any redox pair
            for pair_name, forms in redox_pairs.items():
                # Check if the synonym contains key terms from either form
                for state in ['oxidized', 'reduced']:
                    for known_syn in forms[state]:
                        # Look for key terms that indicate redox state
                        if state == 'oxidized':
                            if any(term in syn_lower for term in ['oxidized', 'oxidation', 'oxo', 'oxo-']):
                                if known_syn.lower() in syn_lower or syn_lower in known_syn.lower():
                                    if pair_name not in potential_additions:
                                        potential_additions[pair_name] = {'oxidized': [], 'reduced': []}
                                    potential_additions[pair_name]['oxidized'].append(syn)
                                    break
                        else:  # reduced
                            if any(term in syn_lower for term in ['reduced', 'reduction', 'hydro', 'dihydro']):
                                if known_syn.lower() in syn_lower or syn_lower in known_syn.lower():
                                    if pair_name not in potential_additions:
                                        potential_additions[pair_name] = {'oxidized': [], 'reduced': []}
                                    potential_additions[pair_name]['reduced'].append(syn)
                                    break
    
    return potential_additions

def validate_redox_mapping(metabolite_mapper):
    """
    Validate that redox pairs are properly separated in the metabolite mapper.
    """
    issues = []
    
    # Create reverse mapping from synonyms to their redox pair
    redox_syn_to_pair = {}
    for pair_name, synonyms in REDOX_PAIRS.items():
        for syn in synonyms:
            redox_syn_to_pair[syn.lower()] = pair_name
    
    # Check each mapping in the metabolite mapper
    for syn, mapped_to in metabolite_mapper['mapping_dictionary'].items():
        syn_lower = syn.lower()
        mapped_lower = mapped_to.lower()
        
        # If this is a redox pair synonym
        if syn_lower in redox_syn_to_pair:
            correct_pair = redox_syn_to_pair[syn_lower]
            
            # Check if it's mapped to a different redox pair
            for other_syn in metabolite_mapper['hmdb_dictionary'].get(mapped_to, []):
                other_lower = other_syn.lower()
                if other_lower in redox_syn_to_pair:
                    mapped_pair = redox_syn_to_pair[other_lower]
                    if mapped_pair != correct_pair:
                        issues.append(f"Cross-mapping detected: {syn} (from {correct_pair}) mapped to {other_syn} (from {mapped_pair})")
    
    return issues

def get_redox_pair_for_metabolite(name, synonyms):
    """
    Determine which redox pair a metabolite belongs to based on its name and synonyms.
    Returns (pair_key, is_oxidized) or (None, None) if not part of a redox pair.
    """
    name_lower = name.lower()
    synonyms_lower = {s.lower() for s in synonyms}
    
    # Check each redox pair
    for pair_key, pair_syns in REDOX_PAIRS.items():
        pair_syns_lower = {s.lower() for s in pair_syns}
        
        # Check if this metabolite matches any synonyms in the pair
        if name_lower in pair_syns_lower or any(s in pair_syns_lower for s in synonyms_lower):
            # Determine if this is the oxidized or reduced form
            is_oxidized = any(ox_pattern in pair_key.lower() for ox_pattern in ['oxidized', '+', '(ox)'])
            return pair_key, is_oxidized
            
    return None, None

def download_hmbd_reference(
        output_dir,
        url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
        file_name='hmdb_metabolites'):
    
    output_file = output_dir + file_name
    print('HMDB metabolite reference:', '\n\t', url)
    print('\tAttempting download...')
    
    try:
        hmdb_url = requests.get(url)
        if hmdb_url.ok:
            print('\tUnzipping...')
            hmdb_zip = zipfile.ZipFile(io.BytesIO(hmdb_url.content))
            hmdb_zip.extractall(output_dir)
            hmdb_zip = None
            print('\tDownload complete.')
            return os.path.join(output_dir, file_name + '.xml')
        else:
            print('\tWarning: Unable to download HMDB data. Will continue with ChEBI as primary source.')
            print('\tTo use HMDB data, please download manually from: https://hmdb.ca/downloads')
            print('\tand place the file in:', output_dir)
            return None
    except Exception as e:
        print('\tWarning: Error downloading HMDB data:', str(e))
        print('\tWill continue with ChEBI as primary source.')
        print('\tTo use HMDB data, please download manually from: https://hmdb.ca/downloads')
        print('\tand place the file in:', output_dir)
        return None

def parse_hmdb_synonyms(
        output_file,
        xml_tag='{http://www.hmdb.ca}'):
    """Retrieve HMDB chemical entity synonyms with strict redox pair handling
    """
    
    print("Parsing HMDB metabolite records...")
    try:
        tree = et.parse(output_file)
        contents = tree.getroot()
        tree = None  # Release the tree object
    except Exception as e:
        print(f"Error parsing XML file: {e}")
        raise

    # Create a mapping from synonyms to their redox pair
    redox_synonym_to_pair = {}
    for pair_key, synonyms in REDOX_PAIRS.items():
        for syn in synonyms:
            redox_synonym_to_pair[syn.lower()] = pair_key
            # Also add alphanumeric version
            redox_synonym_to_pair[''.join(c.lower() for c in syn if c.isalnum())] = pair_key

    hmdb_dictionary = {}
    display_dictionary = {}
    mapping_dictionary = {}
    
    # First pass to build initial dictionaries
    temp_synonyms = defaultdict(set)
    temp_display = defaultdict(set)

    for x in contents:
        if x.tag == xml_tag + 'metabolite':
            name = ''
            all_synonyms = set()
            display_synonyms = set()
            
            # First get all possible names/synonyms
            for y in x:
                if y.text != None:
                    simple_string = ''.join(
                        str(c).lower() for c in y.text if c.isalnum()
                    )
                else:
                    continue

                if y.tag == xml_tag + 'name':
                    name = simple_string
                    all_synonyms.add(simple_string)
                    all_synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

                elif y.tag == xml_tag + 'synonyms':
                    for child in y:
                        if child.text:
                            simple_child = ''.join(
                                str(c).lower() for c in child.text if c.isalnum()
                            )
                            all_synonyms.add(simple_child)
                            all_synonyms.add(str(child.text).lower())
                            display_synonyms.add(str(child.text))

                elif y.tag in [xml_tag + 'iupac_name', xml_tag + 'traditional_iupac']:
                    all_synonyms.add(simple_string)
                    all_synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))
                
                elif y.tag == xml_tag + 'accession':
                    all_synonyms.add(simple_string)
                    all_synonyms.add(str(y.text).lower())
                    all_synonyms.add(str(y.text))
                    display_synonyms.add(str(y.text))

            if name:
                # Determine which redox pair this metabolite belongs to
                pair_key, is_oxidized = get_redox_pair_for_metabolite(name, all_synonyms)
                
                if pair_key:
                    # Only keep synonyms that match this specific redox form
                    allowed_synonyms = {s.lower() for s in REDOX_PAIRS[pair_key]}
                    filtered_synonyms = {s for s in all_synonyms 
                                      if s.lower() in allowed_synonyms or 
                                      ''.join(c.lower() for c in s if c.isalnum()) in allowed_synonyms}
                    filtered_display = {s for s in display_synonyms 
                                     if s.lower() in allowed_synonyms or 
                                     ''.join(c.lower() for c in s if c.isalnum()) in allowed_synonyms}
                else:
                    # Not a redox pair, keep all synonyms
                    filtered_synonyms = all_synonyms
                    filtered_display = display_synonyms

                temp_synonyms[name].update(filtered_synonyms)
                temp_display[name].update(filtered_display)

    # Clear references to XML objects
    contents = None

    # Second pass to build final dictionaries
    for name, synonyms in temp_synonyms.items():
        hmdb_dictionary[name] = sorted(list(synonyms))
        display_dictionary[name] = sorted(list(temp_display[name]))
        for syn in synonyms:
            mapping_dictionary[syn] = name

    # Validate the mapping
    issues = validate_redox_mapping({'mapping_dictionary': mapping_dictionary, 'hmdb_dictionary': hmdb_dictionary})
    if issues:
        print("Warning: Found the following issues in redox pair mapping:")
        for issue in issues:
            print(f"  - {issue}")

    return hmdb_dictionary, display_dictionary, mapping_dictionary


def __main__(
        args_dict):
    """Build metabolite name mapping dictionary
    """

    output_file = download_hmbd_reference(
        output_dir=args_dict['output']
    )

    try:
        if output_file and os.path.exists(output_file):
            hmdb_dictionary, display_dictionary, mapping_dictionary = parse_hmdb_synonyms(
                output_file=output_file
            )
        else:
            print("Using ChEBI as primary metabolite mapping source...")
            hmdb_dictionary = {}
            display_dictionary = {}
            mapping_dictionary = {}

        # Check for potential new redox pair synonyms
        potential_additions = update_redox_pairs(hmdb_dictionary, REDOX_PAIRS)
        if potential_additions:
            print("\nPotential new redox pair synonyms found:")
            for pair_key, new_syns in potential_additions.items():
                if new_syns:
                    print(f"\n{pair_key}:")
                    for syn in new_syns:
                        print(f"  - {syn}")

        mapping_db = {
            'hmdb_dictionary': hmdb_dictionary,
            'display_dictionary': display_dictionary,
            'mapping_dictionary': mapping_dictionary
        }

        print('Writing database to file...')
        # First write the pickle file
        pickle_file = 'metabolite_mapping.pickle'
        write_database(
            output=args_dict['output'],
            file=pickle_file,
            database=mapping_db)
        
        # Create zip file containing the pickle file
        zip_path = os.path.join(args_dict['output'], 'metabolite_mapping.pickle.zip')
        pickle_path = os.path.join(args_dict['output'], pickle_file)
        
        print('Creating zip archive...')
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(pickle_path, pickle_file)
        
        print(f'Successfully created metabolite mapper at: {zip_path}')
        
        # Clean up temporary files
        try:
            os.remove(pickle_path)  # Remove uncompressed pickle file
        except Exception as e:
            print(f"Warning: Could not remove temporary pickle file: {e}")
            
        try:
            if output_file and os.path.exists(output_file):
                os.remove(output_file)  # Remove downloaded HMDB file
        except Exception as e:
            print(f"Warning: Could not remove HMDB XML file: {e}")
            print("You may need to manually delete this file later.")

    except Exception as e:
        print(f"Error during metabolite mapper creation: {e}")
        raise


def test():
    args_dict = {}
    args_dict['cmd'] = 'electrum'
    args_dict['output'] = 'C:\\Users\\jorda\\Desktop\\'
    output_dir = args_dict['output']
    __main__(args_dict)
