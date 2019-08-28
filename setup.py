"""
BioNet-Analyzer
A toolkit for navigating and analyzing gene expression datasets
alias: bionetter
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

"""IMPORT DEPENDENCIES"""
from setuptools import setup
import re
import os
from metabalyzer.__init__ import __version__, __dependencies__

__path__  =  os.path.dirname(os.path.realpath(__file__)) + '/'

"""Setup arguments"""
setup(
    name = 'MetaboNet-Analyzer',
    version = __version__,
    description = 'A toolkit for navigating and analyzing biological networks',
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    author = 'Jordan Berg',
    author_email = 'jordan.berg@biochem.utah.edu',
    url = 'https://github.com/j-berg/MetaboNet-Analyzer',
    packages = ['metabalyzer'],
    exclude= ['tests','docs','recipes'],
    package_dir = {'metabalyzer': 'metabalyzer'},
    license = 'GPL-3.0',
    zip_safe = False,

    install_requires = __dependencies__,

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)
