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

"""Import internal dependencies
"""

"""Per a omics type, generate fold changes and p values based on input data
and metadata
"""

def get_type(
        type):
    """Extract omics type from user option
    - Flexible selection based on first letter alone
    """

    if type[0].lower() == 't':
        print('Provided transcriptomics data')
        type = 'transcriptomics'
    if type[0].lower() == 'p':
        print('Provided proteomics data')
        type = 'transcriptomics'
    if type[0].lower() == 'm':
        print('Provided metabolomics data')
        type = 'transcriptomics'
    else:
        print('Invalid omics type provided. Exiting...')
        sys.exit(1)

    return type

def __main__(
        args_dict):
    """
    """

    type = get_type(
        type=args_dict['type'])
