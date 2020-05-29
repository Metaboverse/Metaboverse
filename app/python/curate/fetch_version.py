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
import re
from urllib.request import Request, urlopen

"""Import internal dependencies
"""
from curate.__init__ import __reactome_version__

"""Set globals
"""
search_string = '<div class=\"moduletable favstats\" ><h3><i class=\"fa fa-calendar\"></i>Version '
url = 'https://reactome.org/'

"""
"""
def check_reactome_version():

    # Check reactome release version
    req = Request(
        url,
        headers={
            'User-Agent': 'Mozilla/5.0'})
    webpage = urlopen(req).read()
    page = webpage.decode('utf-8').splitlines()

    for line in page:

        if search_string in line:
            snippet = re.sub('<[^<]+?>', '', line)
            result = re.search('Version (.*) released', snippet)
            parsed_version = int(result.group(1))

    if parsed_version:
        if parsed_version == __reactome_version__:
            pass
        else:
            print('Warning: Current version of Reactome database in Metaboverse behind most recent Reactome version')
    else:
        print('Warning: Unable to find current Reactome database version')

"""
"""
def __main__():

    check_reactome_version()
