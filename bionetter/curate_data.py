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
from sklearn import preprocessing

from bionetter.utils import file_path, check_suffix

# Generate list of dataframes
# Dataframe is assumed to have columns as sample names and
class dataContainer:

    # Initialize data container with possible dataframes
    def __init__(self):

        self.transcriptome = pd.DataFrame
        self.proteome = pd.DataFrame
        self.metabolome = pd.DataFrame
        self.metadata = pd.DataFrame

    # Input data type
    def add_data(self, file, type):

        # Check that file has full path
        file = file_path(file)

        # Figure out file type
        suffix = check_suffix(file)
        print(suffix)
        # Import dataframe
        data = pd.read_csv(
            file,
            sep = suffix,
            header = 0,
            index_col = 0,
            low_memory = False
        )

        # Populate dataframe in container based on omics type
        if type.lower() == 'transcriptome':
            self.transcriptome = data

        elif type.lower() == 'proteome':
            self.proteome = data

        elif type.lower() == 'metabolome':
            self.metabolome = data

        elif type.lower() == 'metadata':
            self.metadata = data

        else:
            raise Exception('Invalid data type specified')

    # Organize data in order and group replicates
    def organize(self):

        print('working on it ')

    # Normalize dataframes
    def normalize(self, method=None):

        if method == None:
            pass

        elif method == 'z-score':
            self[self.columns] = preprocessing.scale(self[self.columns], axis=1)

        elif method == 'log2fc':
            self
            print('working on it ')
        else:
            raise Exception('Invalid normalization method provided')

    # Impute intermediate time-points for time-course data
    def impute(self, even_spacing=True, iterations=100):
        print('working on it ')




data = dataContainer()
data.add_data('/Users/jordan/Desktop/cccp_riboseq/cccp_count_table.tsv', 'transcriptome')
data.transcriptome.normalize
