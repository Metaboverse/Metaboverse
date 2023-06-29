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

from setuptools import setup
import re
import os

__path__ = os.path.dirname(os.path.realpath(__file__))

"""Get version"""
with open(os.path.join(__path__, 'metaboverse_cli', '__init__.py'), 'r') as fd:
    __version__ = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)

"""Setup arguments"""
setup(
    name='metaboverse-cli',
    version=__version__,
    description='A toolkit for navigating and analyzing biological networks',
    author='Jordan A. Berg',
    author_email='jordan.berg@biochem.utah.edu',
    url='https://github.com/Metaboverse/metaboverse-cli',
    packages=['metaboverse_cli'],
    exclude=[
        os.path.join('metaboverse_cli', 'test'),
        os.path.join('metaboverse_cli', 'mapper', 'test'),
        os.path.join('metaboverse_cli', 'curate', 'test'),
        os.path.join('metaboverse_cli', 'analyze', 'test'),
        'docs'
    ],
    package_dir={'metaboverse_cli': '.'},
    license='GPL-3.0',
    zip_safe=False,
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'scikit-learn',
        'networkx'
    ],
    entry_points={
        "console_scripts": [
            "metaboverse = metaboverse_cli.__main__:main"
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
