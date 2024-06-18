"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) Jordan A. Berg, The University of Utah

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
import os
import pandas as pd


def file_path(input):
    """Check if input contains full path address and return the absolute path."""
    return os.path.abspath(input)


def check_suffix(file):
    """Get file suffix and return the appropriate delimiter."""
    if file.endswith('.csv'):
        return ','
    elif file.endswith(('.tsv', '.txt')):
        return '\t'
    else:
        raise Exception('Invalid data file provided. Expected a tab- or comma-delimited file')


def add_data(file):
    """Read a file and return a pandas DataFrame."""
    file = file_path(file)
    suffix = check_suffix(file)
    data = pd.read_csv(file, sep=suffix, header=0, index_col=0, low_memory=False)
    return data






