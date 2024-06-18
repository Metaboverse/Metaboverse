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
import zipfile 
import pickle
import re 
import os

from utils import progress_feed


def load_references(
        args_dict,
        ensembl,
        uniprot,
        chebi,
        uniprot_metabolites):
    """Load and prepare reference databases"""
    reverse_genes = {v: k for k, v in ensembl.items()}
    protein_dictionary = uniprot_ensembl_reference(
        uniprot_reference=uniprot,
        ensembl_reference=reverse_genes)
    progress_feed(args_dict, "graph", 1)

    chebi_dictionary = build_chebi_reference(
        chebi=chebi,
        uniprot=uniprot_metabolites)

    name_reference = build_name_reference(
        ensembl=ensembl,
        uniprot=uniprot)

    uniprot_mapper = {}
    for k, v in uniprot_metabolites.items():
        uniprot_mapper[v] = k

    return reverse_genes, protein_dictionary, chebi_dictionary, \
        name_reference, uniprot_mapper


def build_name_reference(ensembl, uniprot):
    name_reference = {}
    for k, v in ensembl.items():
        if 'phospho-' in v and '-phospho-' not in v:
            v = v.replace('phospho-', '')
        if '(' in v and ')' in v:
            v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
        if '  ' in v:
            v = v.replace('  ', ' ')
        v = v.strip()
        name_reference[v] = k
        name_reference[k] = k
    for k, v in uniprot.items():
        if 'phospho-' in v and '-phospho-' not in v:
            v = v.replace('phospho-', '')
        if '(' in v and ')' in v:
            v = re.sub("[\(\[].[^a-zA-Z]+?[\)\]]", "", v)
        if '  ' in v:
            v = v.replace('  ', ' ')
        v = v.strip()
        name_reference[v] = k
        name_reference[k] = k

    return name_reference


def build_chebi_reference(chebi, uniprot):
    chebi_dictionary = {}
    for k, v in chebi.items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_dictionary[_k] = v
        chebi_dictionary[k] = v
    for k, v in uniprot.items():
        _k = ''.join(c.lower() for c in str(k) if c.isalnum())
        chebi_dictionary[_k] = v
        chebi_dictionary[k] = v
    return chebi_dictionary


def uniprot_ensembl_reference(uniprot_reference, ensembl_reference):
    """Build cross-referencing dictionary to convert uniprot ID to corresponding Ensembl ID"""
    new_dict = {}
    for k, v in uniprot_reference.items():
        try:
            new_dict[k] = ensembl_reference[v]
        except:
            pass
    return new_dict


def load_metabolite_synonym_dictionary(
        dir=os.path.join(os.path.dirname(__file__), 'data'),
        file='metabolite_mapping.pickle'):
    print("Reading metabolite mapper...")
    with zipfile.ZipFile(dir + os.path.sep + file + '.zip', 'r') as zip_ref:
        metabolite_mapper = pickle.load(
            zip_ref.open(file)
        )
    return metabolite_mapper


def gather_synonyms(
        map_id,
        init_syns,
        metabolite_mapper,
        uniprot_mapper,
        ignore_enantiomers):

    if map_id in uniprot_mapper:
        init_syns.append(uniprot_mapper[map_id])
        init_syns.append(uniprot_mapper[map_id].lower())
        _u = ''.join(
            c.lower() for c in str(uniprot_mapper[map_id]) if c.isalnum())
        init_syns.append(_u)

    parsed_syns = set()
    for s in init_syns:
        parsed_syns.add(s)
        parsed_syns.add(s.lower())
        _s = ''.join(c.lower() for c in str(s) if c.isalnum())
        parsed_syns.add(_s)

        right_splits = [
            ', ',
            ' monocation',
            ' anion',
            ' monoanion',
            ', potassium',
            ', aluminum']
        left_splits = [
            'hydrogen ']
        for r in right_splits:
            if r in s.lower():
                parsed_syns.add(s.lower().split(r)[0])
        for l in left_splits:
            if l in s.lower():
                parsed_syns.add(s.lower().split(l)[1])

        if ignore_enantiomers == True:
            # if _s[0] == 'l' or _s[0] == 'd':
            #    parsed_syns.add(_s[1:])
            if s[0:2] == 'L-' or s[0:2] == 'D-':
                parsed_syns.add(s[2:])
            if s.lower()[0:2] == 'l-' or s.lower()[0:2] == 'd-':
                parsed_syns.add(s.lower()[2:])

    # Step 1: If ignore_enantiomers, trim off from beginning in data table and dictionaries, but search both with and without
    mapper_id = None
    check_keys = []
    search_keys = []
    log_keys = []
    for p in list(parsed_syns):
        if len(p) <= 1:  # ignore electrons, etc in mapping
            pass
        elif p in metabolite_mapper['mapping_dictionary']:
            mapper_id = metabolite_mapper['mapping_dictionary'][p]
        elif map_id in uniprot_mapper:
            if uniprot_mapper[map_id] in metabolite_mapper['mapping_dictionary']:
                mapper_id = metabolite_mapper['mapping_dictionary'][uniprot_mapper[map_id]]
        else:
            pass
        if mapper_id != None:
            log_keys.append(mapper_id)
    if len(log_keys) > 0:
        mapper_id = [word for word, word_count in Counter(
            log_keys).most_common(1)][0]
    else:
        mapper_id = None

    parsed_syns_list = list(parsed_syns)
    if mapper_id != None:
        parsed_syns_list.append(mapper_id)
    if mapper_id in metabolite_mapper['hmdb_dictionary']:
        for m in metabolite_mapper['hmdb_dictionary'][mapper_id]:
            parsed_syns_list.append(m)

    return mapper_id, parsed_syns_list