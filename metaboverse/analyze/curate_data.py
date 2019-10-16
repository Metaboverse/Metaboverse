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
import sys
import pandas as pd
import numpy as np
from pandas.core.base import PandasObject
from sklearn import preprocessing

"""Import internal dependencies
"""
from metaboverse.analyze.utils import import file_path
from metaboverse.analyze.utils import check_suffix
from metaboverse.analyze.utils import add_data
from metaboverse.analyze.utils import format_data
from metaboverse.analyze.utils import format_times
from metaboverse.analyze.utils import even_spacing
from metaboverse.analyze.utils import ratio_spacing
from metaboverse.analyze.utils import sort_columns

# Generate dataframe collection
# Dataframe is assumed to have columns as sample names and analytes as rows
# Dataframe is assumed to be library, etc. normalized previously
# If using riboseq data, TE must already be calculated
# Base treatment is considered
class dataContainer():

    # Initialize data container with possible dataframes
    def __init__(
            self):

        self.transcriptome = pd.DataFrame
        self.proteome = pd.DataFrame
        self.metabolome = pd.DataFrame
        self.metadata = pd.DataFrame
        self.combined = pd.DataFrame

    # Add data objects
    def add_transcriptome(
            self,
            file):

        self.transcriptome = add_data(file)

    def add_proteome(
            self,
            file):

        self.proteome = add_data(file)

    def add_metabolome(
            self,
            file):

        self.metabolome = add_data(file)

    def add_metadata(
            self,
            file):

        self.metadata = add_data(file)

    # Format data objects
    def format_transcriptome(
            self):

        self.transcriptome = format_data(self.transcriptome, self.metadata)

    def format_proteome(
            self):

        self.proteome = format_data(self.proteome, self.metadata)

    def format_metabolome(
            self):

        self.metabolome = format_data(self.metabolome, self.metadata)

    # Normalize dataframes
    def normalize(
            df,
            method=None):

        # Z-score
        if method == 'standard':
            df[df.columns] = preprocessing.scale(df[df.columns], axis=1)

        # Log2 Fold Change compared to base condition (WT or time 0)
        else:
            df[df.columns] = df[df.columns].div(df[df.columns].iloc[:,0], axis=0)
            df[df.columns] = np.log2(df[df.columns] + 1e-7)
            df[df.columns] = df[df.columns].subtract(df[df.columns].iloc[:,0], axis=0)

    PandasObject.normalize = normalize

    # Impute intermediate time-points for time-course data
    def impute(
            self,
            even_spaced=True,
            total_timepoints=100):

        # Get some metadata
        #df = sort_columns(df)
        cols_numeric = [format_times(name) for name in self.columns.tolist()]
        tp_numbers = len(self.columns)
        tp_total = total_timepoints - tp_numbers
        self.columns = cols_numeric

        # Evenly spaced
        if even_spaced == True:
            print('Calculating the same number of imputations between time points...')
            self = even_spacing(
                self,
                cols_numeric,
                tp_numbers,
                tp_total)

        # Ratio spaced
        else:
            print('Calculating a proportional number of imputations between time points...')
            self = ratio_spacing(
                self,
                cols_numeric,
                tp_numbers,
                tp_total)

    PandasObject.impute = impute

def prep_metadata(
        args_dict,
        data):

    if 'metadata' in args_dict \
    and args_dict['rnaseq'] != None:
        data.add_metadata(args_dict['metadata'])

    else:
        print('Metadata required. Exiting...')
        sys.exit(1)

    return data

def prep_rnaseq(
        args_dict,
        data):

    if 'rnaseq' in args_dict \
    and args_dict['rnaseq'] != None:
        data.add_transcriptome(args_dict['rnaseq'])
        data.format_transcriptome()

    return data

def prep_proteomics(
        args_dict,
        data):

    if 'proteomics' in args_dict \
    and args_dict['proteomics'] != None:
        data.add_proteomics(args_dict['proteomics'])
        data.format_proteome()

    return data

def prep_metabolomics(
        args_dict,
        data):

    if 'metabolomics' in args_dict \
    and args_dict['metabolomics'] != None:
        data.add_metabolomics(args_dict['metabolomics'])
        data.format_metabolome()

    return data

def combine_data(
        args_dict,
        data):

    # Combine data
    frames = [
        data.transcriptome,
        data.proteome,
        data.metabolome]
    data.combined = pd.concat(frames)

    # Normalize
    if 'normalize' in args_dict \
    and args_dict['normalize'] == True:
        data.combined.normalize(
            method='standard')
    else:
        data.combined.normalize(
            method='fold')

    if 'time' in args_dict \
    and args_dict['time'] != None:
        data.combined.impute()
        data.combined = sort_columns(data.combined)

    return data.combined

def __main__(
        metadata,
        transcriptomics,
        proteomics,
        metabolomics):

    # Init data container
    data = dataContainer()

    data = prep_metadata(
        args_dict=args_dict,
        data=data)

    # Import RNA-seq/ribo-seq data, if exists
    data = prep_rnaseq(
        args_dict=args_dict,
        data=data)

    # Import proteomics data, if exists
    data = prep_proteomics(
        args_dict=args_dict,
        data=data)

    # Import metabolomics data, if exists
    data = prep_metabolomics(
        args_dict=args_dict,
        data=data)

    # Compile data to single table
    data_combined = combine_data(
        args_dict=args_dict,
        data=data)

    return data_combined
