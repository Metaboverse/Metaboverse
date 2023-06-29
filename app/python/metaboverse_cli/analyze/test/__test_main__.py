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
import zipfile
import importlib.util
import numpy as np
import networkx as nx
import xml.etree.ElementTree as et
import pandas as pd
import pickle
import os

# Run full test on data
print("Testing analyze/__main__.py for modeling data")
spec = importlib.util.spec_from_file_location(
    "", os.path.abspath("./metaboverse_cli/analyze/__main__.py"))
__main__ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(__main__)
test_modeling = __main__.__main__

args_dict = {
    'database_source': 'reactome',
    'curation': 'HSA.mvdb',
    'organism_id': 'HSA',
    'transcriptomics': os.path.abspath(
        './metaboverse_cli/analyze/test/rna_mapping_test.txt'),
    'proteomics': 'none',
    'metabolomics': os.path.abspath(
        './metaboverse_cli/analyze/test/metabolite_mapping_test.txt'),
    'output': os.path.abspath(
        './metaboverse_cli/analyze/test'),
    'output_file': os.path.abspath(
        './metaboverse_cli/analyze/test/HSA_test.mvrs'),
    'collapse_with_modifiers': False,
    'broadcast_genes': True,
    'labels': '0',
    'blocklist': '',
    'organism_curation_file': 'None',
    'force_new_curation': True,
    'collapse_threshold': 0.4,
    'session_data': os.path.abspath(
        './metaboverse_cli/analyze/test/session_data.txt'),
}

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

test_modeling(args_dict)

rna_unmapped = os.path.abspath(
    './metaboverse_cli/analyze/test/rna_mapping_test_unmapped.txt'
)
rna = pd.read_csv(rna_unmapped, sep='\t', index_col=0)
assert len(rna.index.tolist()) == 7028, 'RNA mapping experienced error'
os.remove(rna_unmapped)

metabolite_unmapped = os.path.abspath(
    './metaboverse_cli/analyze/test/metabolite_mapping_test_unmapped.txt'
)
met = pd.read_csv(metabolite_unmapped, sep='\t', index_col=0)
assert met.index.tolist() == ['bMethyl.2.oxovalerate', 'DSS',
                              'Phenylacetylglycine'], 'Metabolite mapping experienced error'
os.remove(metabolite_unmapped)

os.remove(args_dict['output_file'])
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.mvdb'))
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA.nbdb'))
os.remove(os.path.abspath(
    './metaboverse_cli/analyze/test/HSA_template.mvrs'))

print('Tests completed')
