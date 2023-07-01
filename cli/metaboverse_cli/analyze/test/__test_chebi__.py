"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) 2022 Metaboverse

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
import json
import zipfile
import importlib.util
import numpy as np
import networkx as nx
import xml.etree.ElementTree as et
import pandas as pd
import pickle
import os

print("Testing prepare_data.py")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/prepare_data.py"))
prepare_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare_data)

zipped_net = os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.zip')
with zipfile.ZipFile(zipped_net, 'r') as zip_file:
    zip_file.extractall(
        os.path.abspath(
            './metaboverse_cli/analyze/test'))

network_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/HSA.mvdb")
with open(network_url, 'rb') as network_file:
    network = pickle.load(network_file)

# CHEBI mapping
print('Testing analyze/__main__.py for CHEBI mapping...')
metabolomics_url = os.path.abspath(
    "./metaboverse_cli/analyze/test/mixed_chebi_ids.txt")
__main__ = prepare_data.__main__
data, stats, unmapped = __main__(
    network=network,
    transcriptomics_url="None",
    proteomics_url="None",
    metabolomics_url=metabolomics_url)
assert data.shape == (8, 1), "Mixed CHEBI ID mapping failed"
"""
CHEBI:57972	                -0.187287     ->  "L-Ala"
CHEBI:18012	                -0.200046     ->  "FUMA"
CHEBI:30797	                -0.108502     ->  "MAL"
D-glucose	                -0.424363     ->  "Glc"
Fructose-6-phosphate	     0.066205     ->  "Fru(6)P"
L-LacTic Acid	            -0.000151     ->  "LACT"
SuCcINic acId	             0.068246     ->  "SUCCA"
gibberish	                -0.110443     ->  N/A
"""

args_chebi = {
    'database_source': 'reactome',
    'curation': 'HSA.mvdb',
    'organism_curation_file': 'None',
    'organism_id': 'HSA',
    'transcriptomics': 'None',
    'proteomics': 'none',
    'metabolomics': metabolomics_url,
    'output': os.path.abspath(
        './metaboverse_cli/analyze/test'),
    'output_file': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA_test.mvrs'),
    'collapse_with_modifiers': False,
    'broadcast_genes': True,
    'broadcast_metabolites': True,
    'labels': '0',
    'blocklist': '',
    'force_new_curation': True,
    'collapse_threshold': 0.5,
    'session_data': os.path.abspath(
        './metaboverse_cli/analyze/test/session_data.txt'),
}
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/__main__.py"))
__main__ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(__main__)
test_modeling = __main__.__main__

# Download
test_modeling(args_chebi)

with open(args_chebi['output_file']) as f:
    chebi_json = json.load(f)

for n in chebi_json['nodes']:
    if n['name'] == 'L-Ala':
        assert n['values'] == [-0.187286934], "Mixed CHEBI mapping failed"
    if n['name'] == 'SUCCA':
        assert n['values'] == [0.068245646], "Mixed CHEBI mapping failed"
    if n['values'] == [-0.11044283]:
        raise Exception("Mixed CHEBI mapping failed")
    if n['name'] == 'Fru(6)P':
        assert n['values'] == [0.06620505], "Mixed CHEBI mapping failed"
    if n['name'] == 'Glc':
        assert n['values'] == [-0.424362985], "Mixed CHEBI mapping failed"
    if n['name'] == 'Glc':
        assert n['values'] == [-0.424362985], "Mixed CHEBI mapping failed"
    if n['name'] == 'FUMA':
        assert n['values'] == [-0.200045781], "Mixed CHEBI mapping failed"
    if n['name'] == 'MAL':
        assert n['values'] == [-0.10850223], "Mixed CHEBI mapping failed"
    if n['name'] == 'LACT':
        assert n['values'] == [-0.000150599], "Mixed CHEBI mapping failed"

os.remove(args_chebi['output_file'])
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.mvdb'))
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.nbdb'))
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA_template.mvrs'))

print('Tests completed')
