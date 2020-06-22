var assert = require('assert')
var path = require('path')

let motifs_loc = path.join(__dirname, '../js/motifs.js')
var motif_test = require(motifs_loc);
motif_test()

let graph_loc = path.join(__dirname, '../js/graph.js')
var graph_test = require(graph_loc);
graph_test()
