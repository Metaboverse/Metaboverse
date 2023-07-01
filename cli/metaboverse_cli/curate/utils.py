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
import os


def get_table(
        output_dir,
        url,
        column_names,
        organism='Homo sapiens',
        organism_key='organism'):
    """Get reactome table from web
    """

    # chebi_reactome_reactions
    file = unpack_table(
        url=url,
        output_dir=output_dir)

    if isinstance(column_names, list):
        header_type = None
    else:
        header_type = column_names

    data = pd.read_csv(
        file,
        sep='\t',
        header=header_type,
        low_memory=False)

    if isinstance(column_names, list) \
            or organism == None:
        data.columns = column_names
        data_organism = data.loc[data[organism_key] == organism]

    else:
        data_organism = data

    return data_organism


"""Open reactome table from web
"""


def unpack_table(
        url,
        output_dir='./'):

    file = output_dir + url.split('/')[-1]
    os.system('curl -kL ' + url + ' -o "' + file + '"')

    return file
