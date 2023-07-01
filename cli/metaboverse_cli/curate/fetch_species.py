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
from urllib.request import Request, urlopen
import re

"""Set globals
"""
search_string_start = 'title=\"'
search_string_stop = '\">'
avoid_string_start_1 = 'title=\"Show Entries\">'
avoid_string_start_2 = 'title=\"Community library of icons'
avoid_string_start_3 = 'title=\"Back to Top'
stop_string = 'No entries found for this species'
url = 'https://www.reactome.org/content/schema/objects/Species?page='


def fetch_species():

    page_number = 1
    organisms = []
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
                    result = re.search(search_string_start +
                                       '(.*)' + search_string_stop, line)
                    parsed_name = result.group(1)
                    organisms.append(parsed_name)

            else:
                pass

        page_number += 1

    organisms = [x for x in organisms if '<i class=' not in x]

    return organisms


def __main__():
    organisms = fetch_species()
    return organisms
