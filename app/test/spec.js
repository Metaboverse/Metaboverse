const Application = require('spectron').Application
const assert = require('assert')
const electronPath = require('electron') // Require Electron from the binaries included in node_modules.
const path = require('path')

let graph_loc = path.join(__dirname, '../js/graph.js')
var graph_test = require(graph_loc);
graph_test()

let motifs_loc = path.join(__dirname, '../js/motifs.js')
var motif_test = require(motifs_loc);
motif_test()
