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

"""Import internal dependencies
"""
from app.python.__init__ import __version__
from app.python.__init__ import __dependencies__
from app.python.arguments import parse_arguments
from app.python.arguments import get_dependencies
from app.python.preprocess.__main__ import __main__ as preprocess
from app.python.curate.__main__ import __main__ as curate
from app.python.utils import progress_feed

"""Run metaboverse
"""
def main(
        args=None):

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

        print('Curating network model...')
        args_dict = curate(args_dict)
        sys.stdout.flush()


        print('Analyzing data in context of network model...')
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
