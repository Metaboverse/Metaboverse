/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) Jordan A. Berg

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
*/
var assert = require('assert')
var path = require('path')

let motifs_loc = path.join(__dirname, '../js/motifs.js')
var motif_test = require(motifs_loc);
motif_test()

let graph_loc = path.join(__dirname, '../js/graph.js')
var graph_test = require(graph_loc);
graph_test()
