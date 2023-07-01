"""License Information
metaboverse-cli
Back-end CLI Tool for Curating Metabolic Networks for Metaboverse
https://github.com/Metaboverse/metaboverse-cli/
alias: metaboverse-cli

MIT License

Copyright (c) 2022 Metaboverse

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
import importlib.util
import pickle
import os

spec = importlib.util.spec_from_file_location(
    "__main__", os.path.abspath("./metaboverse_cli/mapper/__main__.py"))
mapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mapper)

# Run
args_dict = {
    'output': os.path.abspath(
        os.path.join(".", "metaboverse_cli", "mapper", "test")) + os.path.sep}
mapper.__main__(
    args_dict)

# Check
mapper_url = os.path.abspath(
    os.path.join(".", "metaboverse_cli", "mapper", "test", "metabolite_mapping.pickle"))
with open(mapper_url, 'rb') as mapper_file:
    mapper = pickle.load(mapper_file)

assert list(mapper.keys()) == ['hmdb_dictionary', 'display_dictionary',
                               'mapping_dictionary'], 'metabolite mapper failed to generate dictionaries'

key0 = list(mapper['hmdb_dictionary'].keys())[0]
assert type(mapper['hmdb_dictionary'][key0]
            ) == list, 'HMDB dictionary improperly formatted'

key1 = list(mapper['display_dictionary'].keys())[0]
assert type(mapper['display_dictionary'][key1]
            ) == list, 'Display dictionary improperly formatted'

key2 = list(mapper['mapping_dictionary'].keys())[0]
assert type(mapper['mapping_dictionary'][key2]
            ) == str, 'Mapping dictionary improperly formatted'

# Clean
os.remove(
    os.path.abspath(
        os.path.join(".", "metaboverse_cli", "mapper", "test", "metabolite_mapping.pickle")))

print('Tests completed')
