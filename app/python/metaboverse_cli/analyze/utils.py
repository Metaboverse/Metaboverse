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


def file_path(
        input):
    # Check input is contains full path address

    return os.path.abspath(input)


def check_suffix(
        file):
    # Get file suffix

    if file.split('.')[-1] == 'csv':
        suffix = ','
    elif file.split('.')[-1] == 'tsv':
        suffix = '\t'
    elif file.split('.')[-1] == 'txt':
        suffix = '\t'
    else:
        raise Exception(
            'Invalid data file provided. Expected a tab- or comma-delimited file')

    return suffix


def add_data(
        file):
    # Input data type
    # Check that file has full path
    file = file_path(file)

    # Figure out file type
    suffix = check_suffix(file)

    # Import dataframe
    data = pd.read_csv(
        file,
        sep=suffix,
        header=0,
        index_col=0,
        low_memory=False)

    return data


def convert_rgba(
        rgba_tuples,
        N=255):
    """Convert python RGBA tuple to web-friendly tuple for later viz
    """

    js = []
    for x in rgba_tuples:

        rgba_list = list(x)
        rgba_new = []
        for x in rgba_list[:3]:
            rgba_new.append(int(x * N))

        rgba_new.append(rgba_list[3])

        js.append(tuple(rgba_new))

    return js


def remove_defective_reactions(
        network):
    """
    """
    no_defective_reactions = {}
    for key in network['reaction_database'].keys():
        rxn_name = network['reaction_database'][key]['name'].lower()
        if 'defective' not in rxn_name \
                and 'mutant' not in rxn_name:
            no_defective_reactions[key] = network['reaction_database'][key]

    return no_defective_reactions
