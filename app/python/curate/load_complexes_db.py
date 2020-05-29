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
import pandas as pd

"""Import internal dependencies
"""
from curate.utils import get_table

"""Get tables
"""
def __main__(
        output_dir):

    complex_participants = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/ComplexParticipantsPubMedIdentifiers_human.txt',
        column_names=0)
    os.remove(output_dir + 'ComplexParticipantsPubMedIdentifiers_human.txt')

    complex_pathway = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/Complex_2_Pathway_human.txt',
        column_names=0)
    os.remove(output_dir + 'Complex_2_Pathway_human.txt')

    return {
        'complex_participants': complex_participants,
        'complex_pathway': complex_pathway}
