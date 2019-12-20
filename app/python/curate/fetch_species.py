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

"""Set globals
"""
search_string_start = 'title=\"'
search_string_stop = '\">'
avoid_string_start_1 = 'title=\"Show Entries\">'
avoid_string_start_2 = 'title=\"Community library of icons'
avoid_string_start_3 = 'title=\"Back to Top'
stop_string = 'No entries found for this species'
url = 'https://www.reactome.org/content/schema/objects/Species?page='

"""
"""
def fetch_species():

    page_number = 1
    organisms = []

    # Check reactome release version
    while page_number != 0:

        req = Request(
            url + str(page_number),
            headers={
                'User-Agent': 'Mozilla/5.0'})
        webpage = urlopen(req).read()
        page = webpage.decode('utf-8').splitlines()

        for line in page:

            if stop_string in line:
                return organisms

            elif search_string_start in line:

                if avoid_string_start_1 in line \
                or avoid_string_start_2 in line \
                or avoid_string_start_3 in line:
                    pass

                else:
                    result = re.search(search_string_start + '(.*)' + search_string_stop, line)
                    parsed_name = result.group(1)
                    organisms.append(parsed_name)

            else:
                pass

        page_number += 1

    return organisms

"""
"""
def __main__():

    organisms = fetch_species()

    return organisms
