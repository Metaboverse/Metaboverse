"""License Information
Metabo-verse:
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

"""Specify output sub-directory paths
"""
def set_paths(
        args_dict):

    args_dict['measurements'] = args_dict['output'] + 'measurements/'

    # Step 1

    # Step 2

    # Step 3

    # ...

    return args_dict

"""Preprocess data
"""
def __main__(
        args_dict):

    args_dict = set_paths(args_dict)

    # Set standards for each omic datatype for input
    # - transcriptomics = counts
    # - translatomics = counts formatted for easy TE calculation
    # - proteomics = raw quant
    # - metabolomics = raw quant

    # Run basic quality control on data after normalized

    return args_dict
