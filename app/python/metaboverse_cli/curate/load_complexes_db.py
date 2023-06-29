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
import os

"""Import internal dependencies
"""
try:
    from curate.utils import get_table
except:
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "get_table", os.path.abspath("./metaboverse_cli/curate/utils.py"))
    get_table = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(get_table)
    get_table = get_table.get_table

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
