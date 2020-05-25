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
import pickle
import xml.etree.ElementTree as et
import pandas as pd

def test():

    file_url = "/Users/jordan/Desktop/metaboverse_data/sce_mct1_omics/_data/metabolomics_mct1_timecourse.txt"
    file = pd.read_csv(file_url, sep='\t', index_col=0)
    file.index.name = None
    _idx = file.index.tolist()
    ignore_enantiomers = True

    node_dict = {}
    for index, row in file.iterrows():
        node_dict[index] = {}
        node_dict[index]['values'] = row
        node_dict[index]['mapping_synonyms'] = []
        node_dict[index]['display_synonyms'] = []
        node_dict[index]['chebi_ids'] = set()
        _i = ''.join(str(c).lower() for c in index if c.isalnum()).lower()

        try:
            __i = mapping_db['mapping_dictionary'][_i]
            node_dict[index]['mapping_synonyms'] = mapping_db['hmdb_dictionary'][__i]
            node_dict[index]['display_synonyms'] = mapping_db['display_dictionary'][__i]
        except:

            try:
                # remove last letter
                __i = mapping_db['mapping_dictionary'][_i[:-1]]
                node_dict[index]['mapping_synonyms'] = mapping_db['hmdb_dictionary'][__i]
                node_dict[index]['display_synonyms'] = mapping_db['display_dictionary'][__i]
            except:
                if _i[0] == 'd' or _i[0] == 'l' and ignore_enantiomers == True:
                    try:
                        __i = mapping_db['mapping_dictionary'][_i[1:]]
                        node_dict[index]['mapping_synonyms'] = mapping_db['hmdb_dictionary'][__i]
                        node_dict[index]['display_synonyms'] = mapping_db['display_dictionary'][__i]
                    except:
                        print("Unable to map: ", index)

                else:
                    print("Unable to map: ", index)

    chebi_mapper = {**mapping_db['chebi_dictionary'], **mapping_db['uniprot_metabolites']}

    for k, v in node_dict.items():

        if k in chebi_mapper:
            node_dict[k]['chebi_ids'].add(chebi_mapper[k])


        for s in node_dict[k]['mapping_synonyms']:
            if s in chebi_mapper:
                node_dict[k]['chebi_ids'].add(chebi_mapper[s])

        node_dict[k]['chebi_ids'] = list(node_dict[k]['chebi_ids'])

def write_database(
        output,
        file,
        database):
    """Write reactions database to pickle file
    """

    # Check provided path exists
    if not os.path.isdir(output):
        os.makedirs(output)

    # Clean up path
    dir = os.path.abspath(output) + '/'

    # Write information to file
    with open(dir + file, 'wb') as file_product:
        pickle.dump(database, file_product)

def parse_hmdb_synonyms(
        output_dir,
        url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
        file_name='hmdb_metabolites',
        xml_tag = '{http://www.hmdb.ca}'):
    """Retrieve HMDB chemical entity synonyms
    """

    output_file = output_dir + file_name
    print("Downloading HMDB metabolite reference...")
    os.system('curl -L ' + url + ' -o ' + output_file + '.zip')
    print("Unzipping HMDB metabolite reference...")
    os.system('unzip ' + output_file + '.zip -d ' + output_dir)
    print("Parsing HMDB metabolite records...")
    hmdb_contents = et.parse(output_file + '.xml')
    contents = hmdb_contents.getroot()
    os.remove(output_file + '.zip')
    os.remove(output_file + '.xml')

    hmdb_dictionary = {}
    display_dictionary = {}
    mapping_dictionary = {}

    for x in contents:

        if x.tag == xml_tag + 'metabolite':

            name = ''
            synonyms = set()
            display_synonyms = set()

            for y in x:

                if y.text != None:
                    simple_string = ''.join(
                        str(c).lower() for c in y.text if c.isalnum()
                    )
                else:
                    simple_string = 'None'

                if y.tag == xml_tag + 'name':
                    name = simple_string
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

                if y.tag == xml_tag + 'synonyms':
                    for child in y:
                        simple_child = ''.join(
                            str(c).lower() for c in child.text if c.isalnum()
                        )
                        synonyms.add(simple_child.lower())
                        synonyms.add(str(child.text).lower())
                        display_synonyms.add(str(child.text))

                if y.tag == xml_tag + 'iupac_name':
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

                if y.tag == xml_tag + 'traditional_iupac':
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    display_synonyms.add(str(y.text))

            if name != '':
                hmdb_dictionary[name] = sorted(list(synonyms))
                display_dictionary[name] = sorted(list(display_synonyms))
                for l in hmdb_dictionary[name]:
                    mapping_dictionary[l] = name

    return hmdb_dictionary, display_dictionary, mapping_dictionary

def __main__(
        output_dir):
    """Build metabolite name mapping dictionary
    """

    hmdb_dictionary, display_dictionary, mapping_dictionary = parse_hmdb_synonyms(
        output_dir=output_dir
    )

    mapping_db = {
        'hmdb_dictionary': hmdb_dictionary,
        'display_dictionary': display_dictionary,
        'mapping_dictionary': mapping_dictionary
    }

    print('Writing database to file...')
    write_database(
        output=output_dir,
        file='metabolite_mapping.pickle',
        database=mapping_db
    )
