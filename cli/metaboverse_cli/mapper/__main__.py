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
import pandas as pd
import xml.etree.ElementTree as et
import requests
import zipfile
import pickle
import io
import os

"""Import internal dependencies
"""
try:
    from utils import prepare_output, write_database, write_database_json
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "", os.path.abspath("./metaboverse_cli/utils.py"))
    utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utils)
    prepare_output = utils.prepare_output
    write_database = utils.write_database
    write_database_json = utils.write_database_json


def download_hmbd_reference(
        output_dir,
        url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
        file_name='hmdb_metabolites'):
    
    output_file = output_dir + file_name
    print('HMDB metabolite reference:', '\n\t', url)
    print('\tDownloading...')
    hmdb_url = requests.get(url)
    if hmdb_url.ok:
        print('\tUnzipping...')
        hmdb_zip = zipfile.ZipFile(io.BytesIO(hmdb_url.content))
        hmdb_zip.extractall(output_dir)
        hmdb_zip = None
    else:
        raise Exception("Unable to download file at: " + url)
    print('\tDownload complete.')

    return os.path.join(output_dir, file_name + '.xml')

def parse_hmdb_synonyms(
        output_file,
        xml_tag='{http://www.hmdb.ca}'):
    """Retrieve HMDB chemical entity synonyms
    """

    print("Parsing HMDB metabolite records...")
    hmdb_contents = et.parse(output_file)
    contents = hmdb_contents.getroot()

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
                
                if y.tag == xml_tag + 'accession':
                    synonyms.add(simple_string)
                    synonyms.add(str(y.text).lower())
                    synonyms.add(str(y.text))
                    display_synonyms.add(str(y.text))

            if name != '':
                hmdb_dictionary[name] = sorted(list(synonyms))
                display_dictionary[name] = sorted(list(display_synonyms))
                for l in hmdb_dictionary[name]:
                    mapping_dictionary[l] = name

    return hmdb_dictionary, display_dictionary, mapping_dictionary


def __main__(
        args_dict):
    """Build metabolite name mapping dictionary
    """

    output_file = download_hmbd_reference(
        output_dir=args_dict['output']
    )

    hmdb_dictionary, display_dictionary, mapping_dictionary = parse_hmdb_synonyms(
        output_file=output_file
    )

    mapping_db = {
        'hmdb_dictionary': hmdb_dictionary,
        'display_dictionary': display_dictionary,
        'mapping_dictionary': mapping_dictionary
    }

    print('Writing database to file...')
    write_database(
        output=args_dict['output'],
        file='metabolite_mapping.pickle',
        database=mapping_db)
    
    os.remove(output_file)


def test():

    args_dict = {}
    args_dict['cmd'] = 'electrum'
    args_dict['output'] = 'C:\\Users\\jorda\\Desktop\\'
    output_dir = args_dict['output']
    __main__(args_dict)
