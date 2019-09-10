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
import os 
import pickle

"""Import internal dependencies
"""
from metaboverse.metaboverse_curate.load_reactions_db import __main__ as load_reactions

"""Plan
"""
# This file will run loading of reactions database, chebi, ensemble, uniprot, complex, etc.
# Will then create interface dictionary for metabolites, proteins, etc relation info, name to id, etc.
# Output total network as pickle

"""Write reactions database to pickle file
"""
def write_database(
        output,
        file,
        database):

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    dir = os.path.abspath(output) + '/'

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)

"""Curate reactome database
"""
def __main__(
        args_dict):

    # Load reactions
    reactions_database = load_reactions(
        species_id=args_dict['species'],
        output_dir=args_dict['output'])

    # Add interfacing information to reactions database

    # Write database to file
    write_database(
        output=args_dict['output'],
        file=args_dict['species'] + '_metaboverse_db.pickle',
        database=reactions_database)
