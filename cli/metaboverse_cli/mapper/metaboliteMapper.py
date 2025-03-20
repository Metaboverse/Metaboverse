"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah
"""
from __future__ import print_function
import os
import sys
import pickle
import tempfile
import zipfile
from pathlib import Path

def get_project_root():
    """Get the path to the project root directory"""
    current_file = Path(__file__).resolve()
    # Go up until we find the cli directory or reach the root
    for parent in current_file.parents:
        if parent.name == 'cli':
            return parent
        # Also check for the repository root as fallback
        if parent.name == 'Metaboverse':
            return parent / 'cli'
    return current_file.parent.parent.parent

# Add the project root and its parent to the Python path
project_root = get_project_root()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
if str(project_root.parent) not in sys.path:
    sys.path.insert(0, str(project_root.parent))

# Try different import strategies
try:
    # First try normal package import
    from metaboverse_cli.mapper.__main__ import parse_hmdb_synonyms, download_hmbd_reference, update_redox_pairs
    from metaboverse_cli.mapper.special_pairs import REDOX_PAIRS
    from metaboverse_cli.utils import write_database
except ImportError:
    try:
        # Then try relative import
        from mapper.__main__ import parse_hmdb_synonyms, download_hmbd_reference, update_redox_pairs
        from special_pairs import REDOX_PAIRS
        from utils import write_database
    except ImportError:
        try:
            # Finally try direct import from current directory
            import importlib.util
            def load_module(name, path):
                spec = importlib.util.spec_from_file_location(name, path)
                if spec is None:
                    raise ImportError(f"Could not find module at {path}")
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                return module

            base_path = os.path.dirname(__file__)
            parent_path = os.path.dirname(os.path.dirname(base_path))

            main = load_module("__main__", os.path.join(base_path, "__main__.py"))
            parse_hmdb_synonyms = main.parse_hmdb_synonyms
            download_hmbd_reference = main.download_hmbd_reference
            update_redox_pairs = main.update_redox_pairs

            special_pairs = load_module("special_pairs", os.path.join(base_path, "special_pairs.py"))
            REDOX_PAIRS = special_pairs.REDOX_PAIRS

            utils = load_module("utils", os.path.join(parent_path, "metaboverse_cli", "utils.py"))
            write_database = utils.write_database

        except ImportError as e:
            print(f"Error importing dependencies: {e}")
            print(f"Current sys.path: {sys.path}")
            print(f"Current directory: {os.getcwd()}")
            print(f"File location: {__file__}")
            raise

def ensure_metabolite_mapper(output_dir=None, force_rebuild=False):
    """
    Ensures a metabolite mapper exists, creating it if necessary.
    Returns the path to the metabolite mapper file.
    
    Args:
        output_dir (str): Directory to store the mapper. If None, uses a temp directory
        force_rebuild (bool): If True, rebuilds the mapper even if it exists
    
    Returns:
        str: Path to the metabolite mapper file
    """
    if output_dir is None:
        output_dir = tempfile.gettempdir()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    mapper_file = os.path.join(output_dir, 'metabolite_mapping.pickle')
    zip_file = mapper_file + '.zip'
    
    # Check if zipped mapper exists and we don't need to rebuild
    if os.path.exists(zip_file) and not force_rebuild:
        return mapper_file
    
    # Create the mapper
    print("Building metabolite mapper...")
    
    # Try to download HMDB data
    try:
        output_file = download_hmbd_reference(output_dir=output_dir)
        if output_file is None or not os.path.exists(output_file):
            print("Warning: HMDB download failed or file not found. Will continue with ChEBI as primary source.")
            hmdb_dictionary = {}
            display_dictionary = {}
            mapping_dictionary = {}
        else:
            print("Successfully downloaded HMDB data. Parsing metabolite records...")
            hmdb_dictionary, display_dictionary, mapping_dictionary = parse_hmdb_synonyms(
                output_file=output_file
            )
            # Clean up HMDB file if it exists
            if os.path.exists(output_file):
                os.remove(output_file)
    except Exception as e:
        print(f"Warning: Error processing HMDB data: {str(e)}")
        print("Will continue with ChEBI as primary source.")
        hmdb_dictionary = {}
        display_dictionary = {}
        mapping_dictionary = {}

    # Check for potential new redox pair synonyms
    if hmdb_dictionary:  # Only check if we have HMDB data
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

    print('Writing metabolite mapper to file...')
    # First write the pickle file
    write_database(
        output=output_dir,
        file='metabolite_mapping.pickle',
        database=mapping_db)
    
    # Create zip file containing the pickle file
    print('Creating zip archive...')
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(mapper_file, 'metabolite_mapping.pickle')
    
    # Remove the unzipped pickle file
    if os.path.exists(mapper_file):
        os.remove(mapper_file)
    
    return mapper_file

def __main__(args_dict):
    """
    Main function for the metaboliteMapper module.
    """
    # Ensure output directory exists
    output_dir = args_dict.get('output', os.getcwd())
    os.makedirs(output_dir, exist_ok=True)
    
    # Force rebuild mapper in output directory
    mapper_file = ensure_metabolite_mapper(output_dir=output_dir, force_rebuild=True)
    print(f"Metabolite mapper created successfully at: {mapper_file}")
    return 0 