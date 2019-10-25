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
from curation_tests.unit_tests import __main__ as curation_tests_unit
from curation_tests.module_tests import __main__ as curation_tests_module
from preprocessing_tests.unit_tests import __main__ as preprocessing_tests_unit
from preprocessing_tests.module_tests import __main__ as preprocessing_tests_module
from analysis_tests.unit_tests import __main__ as analysis_tests_unit
from analysis_tests.module_tests import __main__ as analysis_tests_module
from viz_tests.unit_tests import __main__ as viz_tests_unit
from viz_tests.module_tests import __main__ as viz_tests_module

def run_curation_tests():

    curation_tests_unit()
    curation_tests_module()

def run_preprocessing_tests():

    preprocessing_tests_unit()
    preprocessing_tests_module()

def run_analysis_tests():

    analysis_tests_unit()
    analysis_tests_module()

def run_viz_tests():

    viz_tests_unit()
    viz_tests_module()

def __main__():

    run_curation_tests()
    run_preprocessing_tests()
    run_analysis_tests()
    run_viz_tests()

__main__()
