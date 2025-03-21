from __future__ import print_function

"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah
"""

__version__ = '0.11.6'

import os
import sys
from pathlib import Path

def get_package_root():
    """Get the path to the package root directory"""
    if getattr(sys, 'frozen', False):
        # Running in a PyInstaller bundle
        return Path(sys._MEIPASS) / 'metaboverse_cli'
    else:
        # Running in development
        return Path(__file__).parent

def get_source_url():
    """Get the source URL from the source_url.txt file"""
    try:
        # First try to read from the package
        source_url_path = get_package_root() / 'source_url.txt'
        if source_url_path.exists():
            return source_url_path.read_text().strip()
        
        # Fallback to development path
        dev_path = Path(__file__).parent / 'source_url.txt'
        if dev_path.exists():
            return dev_path.read_text().strip()
            
        raise FileNotFoundError("Could not find source_url.txt in any expected location")
    except Exception as e:
        print(f"Error reading source_url.txt: {e}")
        print(f"Current directory: {os.getcwd()}")
        print(f"Package root: {get_package_root()}")
        print(f"Development path: {Path(__file__).parent}")
        raise

# Initialize global variables
SOURCE_URL = get_source_url()

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
