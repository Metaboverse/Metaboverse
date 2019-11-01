import os
import pandas as pd
import numpy as np
import pickle


# Get list for drop-down selection of all analytes for targeted analysis
def get_list_all_analytes():

    return

# Allow for variable naming of analytes for Metaboverse plotting
# Use CAS database, or KEGG name as dictionary to map
# If not included, make note and have user manually convert name (once entered, check dictionary)
# If not converted, will be left out of analysis
# Check two studies from metabolomics workbench from different labs (one known, one unknown, and compare naming to see in DB requires standardization)
def map_name_to_standard(
        name):

    metaboverse_name = name

    return metaboverse_name

# UI will have drop-down menu with list of analytes to choose from (full name (abbreviation))
# Map entity name to all IDs
# Find all pathways with entity involved
# Return list of pathways with hyperlinks for user to then viz
def target_analyte(
        name,
        reference):

    # Get Metaboverse-safe entity name
    name = map_name_to_standard(
        name=name)

    # Get all relevant Reactome IDs
    target_ids = []

    for analyte_id in reference['master_reference'].keys():
        if str(reference['master_reference'][analyte_id]).lower() == str(target).lower():
            target_ids.append(analyte_id)

    # Get pathways
    pathway_ids = []
    pathway_names = []

    for pathway_id in reference['pathways'].keys():

        for reactions in reference['pathways'][pathway_id]['reactions'].keys():

            analytes = []
            for x in reference['pathways'][pathway_id]['reactions'][reactions]['reactants'].keys():
                analytes.append(x)

            for y in reference['pathways'][pathway_id]['reactions'][reactions]['products'].keys():
                analytes.append(y)

            for z in reference['pathways'][pathway_id]['reactions'][reactions]['modifiers'].keys():
                analytes.append(z)

            if len(set(analytes).intersection(set(target_ids))) > 0:
                pathway_ids.append(pathway_id)
                pathway_names.append(reference['pathways'][pathway_id]['pathway_name'])
                break

        continue

    return set(pathway_names), set(pathway_ids)


# Test entity mapping and selection
def test_target():

    output = '/Users/jordan/Desktop/Metaboverse/tests/analysis_tests/'
    with open(output + 'HSA_metaboverse_db.pickle', 'rb') as network_file:
        reference = pickle.load(network_file)

    target = 'H+'
    pathway_names, pathway_ids = target_analyte(
        name=target,
        reference=reference)

def test_map():

    test_file = '/Users/jordan/Desktop/test.txt'
    data = pd.read_csv(
        test_file,
        sep='\t',
        header=None,
        index_col=0)
    del data.index.name
    data

for x in reference['master_reference'].keys():
    if 'ENS' in reference['master_reference'][x]:
        print(x)
        print(reference['master_reference'][x])


reference['master_reference']
