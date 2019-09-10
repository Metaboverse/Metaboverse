import os
import xml.etree.ElementTree as et
import networkx as nx
import pandas as pd
import pickle
import sys
import time

__path__ = '/Users/jordan/Desktop/reactome_test/'



# chebi_all
chebi_file_all = __path__ + 'ChEBI2Reactome_All_Levels.txt'
chebi_all = pd.read_csv(
    chebi_file_all,
    sep='\t',
    header=None)
chebi_all.columns = [
    '?',
    'reactome_id',
    'url',
    'process',
    'go_evidence', #TAS = traceable author statement, IEA = electrong annotation not manually reviewed
    'organism']
chebi_all.shape

chebi_all.loc[chebi_all['organism'] == 'Homo sapiens'].head()

set(chebi_all.process.tolist())

# chebi_pe_all
chebi_file_peall = __path__ + 'ChEBI2Reactome_PE_All_Levels.txt'
chebi_peall = pd.read_csv(
    chebi_file_peall,
    sep='\t',
    header=None)
chebi_peall.columns = [
    '?',
    'general_reactome_id',
    'metabolite_name',
    'reactome_id',
    'url',
    'process',
    'go_evidence',
    'organism']
chebi_peall.shape

chebi_peall.loc[chebi_peall['organism'] == 'Homo sapiens'].head()

# chebi_pe_path
chebi_file_path = __path__ + 'ChEBI2Reactome_PE_Pathway.txt'
chebi_path = pd.read_csv(
    chebi_file_path,
    sep='\t',
    header=None)
chebi_path.columns = [
    '?',
    'general_reactome_id',
    'metabolite_name',
    'path_id',
    'url',
    'process',
    'go_evidence',
    'organism']
chebi_path.shape

chebi_path.loc[chebi_path['organism'] == 'Homo sapiens'].head()

# chebi_pe_reactions
chebi_file_reaction = __path__ + 'ChEBI2Reactome_PE_Reactions.txt'
chebi_reaction = pd.read_csv(
    chebi_file_reaction,
    sep='\t',
    header=None)
chebi_reaction.columns = [
    '?',
    'general_reactome_id',
    'metabolite_name',
    'reaction_id',
    'url',
    'reaction',
    'go_evidence',
    'organism']
chebi_reaction.shape

chebi_reaction.loc[chebi_reaction['organism'] == 'Homo sapiens'].head()

chebi_reaction.loc[(chebi_reaction['organism'] == 'Homo sapiens') & (chebi_reaction['reaction'] == 'ITPA hydrolyses XTP to XMP')]

# chebi_reactome
chebi_file_react = __path__ + 'ChEBI2Reactome.txt'
chebi_react = pd.read_csv(
    chebi_file_react,
    sep='\t',
    header=None)
chebi_react.columns = [
    '?',
    'reactome_id',
    'url',
    'process',
    'go_evidence',
    'organism']
chebi_react.shape

chebi_react.loc[chebi_react['organism'] == 'Homo sapiens'].head()

# chebi_reactom_reactions
chebi_file_react_reactions = __path__ + 'ChEBI2ReactomeReactions.txt'
chebi_react_reactions = pd.read_csv(
    chebi_file_react_reactions,
    sep='\t',
    header=None)
chebi_react_reactions.columns = [
    '?',
    'reactome_id',
    'url',
    'reaction',
    'go_evidence',
    'organism']
chebi_react_reactions.shape

chebi_react_reactions.loc[chebi_react_reactions['organism'] == 'Homo sapiens'].head()


























import json

human_json = __path__ + 'Homo_sapiens.json'
with open(human_json) as json_file:
    human_model = json.load(json_file)

human_model
