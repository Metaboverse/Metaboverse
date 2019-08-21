"""
BioNet-Analyzer
A toolkit for navigating and analyzing gene expression datasets
alias: bionetter
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

import pandas as pd
import numpy as np
from pandas.core.base import PandasObject
from sklearn import preprocessing

from bionetter.utils import file_path, check_suffix, add_data, format_data, format_times, even_spacing, ratio_spacing, sort_columns

# Generate dataframe collection
# Dataframe is assumed to have columns as sample names and analytes as rows
# Dataframe is assumed to be library, etc. normalized previously
# If using riboseq data, TE must already be calculated
# Base treatment is considered
class dataContainer():

    # Initialize data container with possible dataframes
    def __init__(self):

        self.transcriptome = pd.DataFrame
        self.proteome = pd.DataFrame
        self.metabolome = pd.DataFrame
        self.metadata = pd.DataFrame

    # Add data objects
    def add_transcriptome(self, file):

        self.transcriptome = add_data(file)

    def add_proteome(self, file):

        self.proteome = add_data(file)

    def add_metabolome(self, file):

        self.metabolome = add_data(file)

    def add_metadata(self, file):

        self.metadata = add_data(file)

    # Format data objects
    def format_transcriptome(self):

        self.transcriptome = format_data(self.transcriptome, self.metadata)

    def format_proteome(self):

        self.proteome = format_data(self.proteome, self.metadata)

    def format_metabolome(self):

        self.metabolome = format_data(self.metabolome, self.metadata)

    # Normalize dataframes
    def normalize(df, method=None):

        if method == None:
            pass

        # Z-score
        elif method == 'standard':
            df[df.columns] = preprocessing.scale(df[df.columns], axis=1)

        # Log2 Fold Change compared to base condition (WT or time 0)
        elif method == 'fold':
            df[df.columns] = df[df.columns].div(df[df.columns].iloc[:,0], axis=0)
            df[df.columns] = np.log2(df[df.columns] + 1e-7)
            df[df.columns] = df[df.columns].subtract(df[df.columns].iloc[:,0], axis=0)

        else:
            raise Exception('Invalid normalization method provided')

    PandasObject.normalize = normalize

    # Impute intermediate time-points for time-course data
    def impute(self, even_spaced=True, total_timepoints=100):

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




### TEST SPACE ###
data = dataContainer()
data.add_transcriptome('/Users/jordan/Desktop/BioNet-Analyzer/tests/test_data/rnaseq_test.txt')
data.add_metadata('/Users/jordan/Desktop/BioNet-Analyzer/tests/test_data/rnaseq_test_meta.txt')
data.format_transcriptome()
data.transcriptome.normalize(method='fold')
data.transcriptome.impute()
data.transcriptome = sort_columns(data.transcriptome)

data.transcriptome.head()
