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

"""Import internal dependencies
"""
from metaboverse.metaboverse_curate.load_reactions_db import __main__ as load_reactions
from metaboverse.metaboverse_curate.load_chebi_db import __main__ as load_chebi
from metaboverse.metaboverse_curate.load_uniprot_db import __main__ as load_uniprot
from metaboverse.metaboverse_curate.load_ensembl_db import __main__ as load_ensembl
from metaboverse.metaboverse_curate.load_ncbi_db import __main__ as load_ncbi
from metaboverse.metaboverse_curate.load_mirbase_db import __main__ as load_mirbase
from metaboverse.metaboverse_curate.load_complexes_db import __main__ as load_complexes
from metaboverse.metaboverse_curate.__main__ import __main__ as curate
