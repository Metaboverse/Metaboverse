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
import pandas as pd

"""Import internal dependencies
"""
from metaboverse.metaboverse_curate.utils import get_table

"""Get tables
"""
def __main__(
        output_dir):

    all_levels = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'url',
            'reaction_id',
            'go_evidence', #TAS = traceable author statement, IEA = electrong annotation not manually reviewed
            'organism'])

    pe_all_levels = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2Reactome_PE_All_Levels.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    pe_pathways = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2Reactome_PE_Pathway.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    pe_reactions = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2Reactome_PE_Reactions.txt',
        column_names=[
            'source_id',
            'analyte_id',
            'analyte_name',
            'reaction_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    reactome = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2Reactome.txt',
        column_names=[
            'source_id',
            'process_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    reactome_reactions = get_table(
        output_dir=output_dir,
        url='https://reactome.org/download/current/UniProt2ReactomeReactions.txt',
        column_names=[
            'source_id',
            'process_id',
            'url',
            'reaction_name',
            'go_evidence',
            'organism'])

    return {
        'uniprot_all_levels': all_levels,
        'uniprot_pe_all_levels': pe_all_levels,
        'uniprot_pe_pathways': pe_pathways,
        'uniprot_pe_reactions': pe_reactions,
        'uniprot_reactome': reactome,
        'uniprot_reactome_reactions': reactome_reactions}
