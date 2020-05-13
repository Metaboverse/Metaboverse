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
import sys
import pickle

"""Import internal dependencies
"""
from python.__init__ import __version__
from python.__init__ import __dependencies__
from python.arguments import parse_arguments
from python.arguments import get_dependencies
from python.preprocess.__main__ import __main__ as preprocess
from python.curate.__main__ import __main__ as curate
from python.analyze.__main__ import __main__ as analyze
from python.utils import progress_feed, update_session

def check_dependencies():

    try:
        import pandas
    except:
        os.system('pip install pandas')
    try:
        import numpy
    except:
        os.system('pip install numpy')
    try:
        import scipy
    except:
        os.system('pip install scipy')
    try:
        import networkx
    except:
        os.system('pip install networkx')
    try:
        import pickle
    except:
        os.system('pip install pickle')
    try:
        import sklearn
    except:
        os.system('pip install scikit-learn')
    try:
        import argparse
    except:
        os.system('pip install argparse')
    try:
        import matplotlib
    except:
        os.system('pip install matplotlib<3.0.0,>=2.1.1')

"""Run metaboverse
"""
def main(
        args=None):

    check_dependencies()

    # Read in arguments
    args, args_dict = parse_arguments(
        args,
        __version__)

    if args_dict['cmd'] == 'preprocess':

        print('Preprocessing ' + args_dict['type'] + ' dataset...')
        preprocess(args_dict)
        sys.stdout.flush()
        sys.exit(1)

    # Run metaboverse-curate
    elif args_dict['cmd'] == 'curate':

        print(args_dict)
        if str(args_dict['organism_curation']) != 'None':
            progress_feed(
                args_dict=args_dict,
                process="curate",
                amount=50)
            # Update args_dict with path for network model
            with open(args_dict['organism_curation'], 'rb') as network_file:
                network = pickle.load(network_file)
                args_dict['species_id'] = network['species_id']
                args_dict['output_file'] = args_dict['output'] \
                    + args_dict['species_id'] \
                    + '_global_reactions.json'
                args_dict['network'] = args_dict['organism_curation']

                # add variables back to session data json file
                session_file = args_dict['session_data']
                update_session(
                    session_file=session_file,
                    key='species_id',
                    value=args_dict['species_id'])
                update_session(
                    session_file=session_file,
                    key='output_file',
                    value=args_dict['output_file'])
                update_session(
                    session_file=session_file,
                    key='database_url',
                    value=args_dict['output_file'])


            print('Skipping organism network modeling as one was provided by' \
            + ' the user...')
            sys.stdout.flush()
        else:
            print('Curating network model...')
            args_dict['network'] = args_dict['model_file']
            args_dict = curate(args_dict)
            sys.stdout.flush()

        print('Curating data onto the network model...')
        analyze(args_dict)
        sys.stdout.flush()
        sys.exit(1)

    # Print some error messaging
    else:
        raise Exception('Invalid sub-module selected')

    # Exit
    # Check log file for errors and exceptions
    #get_dependencies(args_dict)
    sys.stdout.flush()

"""Run main
"""
if __name__ == '__main__':

    sys.exit(main() or 0)
