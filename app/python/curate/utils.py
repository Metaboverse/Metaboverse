"""License Information
Metaboverse:
    A toolkit for navigating and analyzing gene expression datasets
    alias: metaboverse
    Copyright (C) 2019  Jordan A. Berg
    jordan <dot> berg <at> biochem <dot> utah <dot> edu

    This program is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License along with
    this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import print_function

"""Import dependencies
"""
import os
import sys
import pandas as pd

"""Get reactome table from web
"""
def get_table(
        output_dir,
        url,
        column_names,
        organism='Homo sapiens',
        organism_key='organism'):

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

    if output_dir[-1] != '/':
        output_dir = output_dir + '/'

    file = output_dir + url.split('/')[-1]

    os.system('curl -L ' + url + ' -o ' + file)

    return file
