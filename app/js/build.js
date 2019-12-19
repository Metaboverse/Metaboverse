/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/j-berg/Metaboverse/
alias: metaboverse

Copyright (C) 2019 Jordan A. Berg
Email: jordan<dot>berg<at>biochem<dot>utah<dot>edu

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

const {app, BrowserWindow} = require('electron')
const {ipcRenderer, ipcMain, remote} = require('electron')
const { dialog } = require('electron').remote
var fs = require('fs')

var $ = require('jquery')

var i = 0;
function move(_callback) {
  if (i == 0) {
    i = 1;
    var elem = document.getElementById("progressBar");
    var width = 10;
    var id = setInterval(frame, 10);
    function frame() {
      if (width >= 100) {
        clearInterval(id);
        i = 0;
        _callback();
      } else {
        width++;
        elem.style.width = width + "%";
        elem.innerHTML = width + "%";
      }
    }
  }
}

var transcriptomics = false;
var proteomics = false;
var metabolomics = false;

var t_val = getDefault("transcriptomics")
var p_val = getDefault("proteomics")
var m_val = getDefault("metabolomics")

if (t_val !== null) {
  transcriptomics = true;
}

if (p_val !== null) {
  proteomics = true;
}

if (m_val !== null) {
  metabolomics = true;
}

function curateNetwork() {

  // Replace with Metaboverse orbiting network with load bar in the middle
  move(() => displayOptions())
  // Run curation here with available session data

}

function displayOptions() {

  if (transcriptomics === true | proteomics === true | metabolomics === true) {

    $('#content').replaceWith('<a href="motif.html"><div id="continue"><font size="3">Run Motif Search</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>')

  } else {

    $('#content').replaceWith('<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>')

  }
}
