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
import pandas as pd

COLUMN_REFERENCE = {
    'Metabolite': 0,
    'Query_protein': 1,
    'Protein_name': 2,
    'Gene_ID': 3,
    'Uniprot_ID': 4,
    'Protein_complex': 5,
    'log2_abundance': 6,
    'log2_abundance_corrected': 7,
    'met_mean': 8,
    'met_sd': 9,
    'p_value': 10,
    'q_value': 11,
    'Common_metabolite_name': 12,
    'MIDAS_ID': 13,
    'Pool': 14,
    'Screened_concentration_µM': 15,
    'KEGG_ID': 16,
    'HMDB_ID': 17,
    'SMILES_metabolite': 18,
    'KEGG_pathway_association': 19
}


def import_midas(
        filename,
        delimiter='\t'):
    """Import MIDAS database from flat format file
    """

    data = pd.read_csv(
        filename,
        sep=delimiter,
        encoding='unicode_escape')

    if len(data.columns.tolist()) != len(list(COLUMN_REFERENCE.keys())):
        raise Exception('Improperly formatted datatable provided: ', filename)

    return data, COLUMN_REFERENCE
