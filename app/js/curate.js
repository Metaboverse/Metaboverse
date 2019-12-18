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


Portions of the force graphing below based on or adapted from code from Mike Bostock
The original code is under the GNU General Public License v3.0, allowing for modification
and distribution
License and copyright notice: GNU General Public License v3.0
Changes:
  - Heavily modified and added to the style CSS for more flexibility in plotting
  - Adapted general D3 plotting functions and commands to work with input data and accept flexibility
  - Modified plotting functions to allow for the differential shading of nodes
  - All other components are original
Source:
http://bl.ocks.org/mbostock/1153292
https://bl.ocks.org/mbostock/1212215
*/

const {app, BrowserWindow} = require('electron')
const {ipcRenderer, ipcMain, remote} = require('electron')
const { dialog } = require('electron').remote
var fs = require('fs')


var $ = require('jquery')
var reactome_api = "https://reactome.org/ContentService/data/species/all";

$.getJSON(reactome_api, function(data) {

  // Get species name and ID from Reactome API
  var abbreviation_dict = {}
  data.forEach(function(datum) {

    abbreviation_dict[datum["displayName"]] = datum["abbreviation"]

  });

  // Get species names (keys) as list
  speciesList = Object.getOwnPropertyNames(
    abbreviation_dict
  ).map(function(k) {
    return k;
  });
  speciesList.unshift("Select an organism..."); // Add select prompt to menu bar

  // Generate drop-down menu for species select
  var menu = document.getElementById("speciesMenu");
  for (var i = 0; i < speciesList.length; i++) {
    var option = document.createElement("option");
    option.innerHTML = speciesList[i];
    option.value = speciesList[i];
    menu.appendChild(option);
  };

});

// Change user selection based on input
function selectOrganism() {

  var selection = document.getElementById("speciesMenu").value;
  console.log("User selected:", selection)
  update_session_info("organism", selection)

};

// Select output directory from pop-out menu
window.addEventListener('load', function(event) {

  document.getElementById("selectOutput").onclick = function (event) {

    filename = dialog.showSaveDialog(
      {
        "defaultPath": "./",
        "properties": ["createDirectory"],
        "filters": [
          {
            "name": "JSON",
            "extensions": ["json", "JSON"]
          }]
      }
    ).then(result => {

      filename = result.filePath;
      if (filename === undefined) {
        alert('File selection unsuccessful');
        return;
      }

      console.log(filename)
      update_session_info("database_url", filename)

      //path = filename.substring(0, filename.lastIndexOf("/"));
      update_session_info("output", "try")
      console.log("hello")
    }).catch(err => {
      console.log(err)
    })

  }
})
