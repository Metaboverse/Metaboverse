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

const { ipcRenderer, ipcMain, remote } = require("electron");
const { dialog } = require("electron").remote;

var $ = require("jquery");
var reactome_api = "https://reactome.org/ContentService/data/species/all";

var abbreviation_dict = {};
$.getJSON(reactome_api, function(data) {
  // Get species name and ID from Reactome API

  data.forEach(function(datum) {
    abbreviation_dict[datum["displayName"]] = datum["abbreviation"];
  });

  // Get species names (keys) as list
  speciesList = Object.getOwnPropertyNames(abbreviation_dict).map(function(k) {
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
  }
});

// Change user selection based on input
function selectOrganism() {
  var selection = document.getElementById("speciesMenu").value;
  console.log("User selected:", abbreviation_dict[selection]);
  update_session_info("organism", selection, (abbrev_dict = abbreviation_dict));
  species_change = true;
  check_changes();
}

// Select output directory from pop-out menu
window.addEventListener("load", function(event) {
  document.getElementById("output-input").onclick = function(event) {
    filename = dialog
      .showSaveDialog({
        defaultPath: "../../",
        properties: ["createDirectory"],
        filters: [
          {
            name: "JSON",
            extensions: ["json", "JSON"]
          }
        ]
      })
      .then(result => {
        filename = result.filePath;
        if (filename === undefined) {
          alert("File selection unsuccessful");
          return;
        }

        console.log(filename);
        update_session_info("database_url", filename);
        $('#selectedOutput').html('<font size="2">' + filename + '</font>');
      })
      .catch(err => {
        console.log(err);
      });

    output_change = true;
    check_changes();
  };
});

var output_change = false;
var species_change = false;

function check_changes() {
  if ((output_change === true) & (species_change === true)) {
    $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a>');
  }
}

// Drop pre-existing metabolic network curation for further analysis
window.addEventListener("load", function(event) {
  document.getElementById("curation-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("curation-input").value.split(".");

    if (inputVal[inputVal.length - 1] !== "pickle") {
      alert(
        "Input is not a .pickle file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        update_session_info("curation_url", f.path);
        $('#selectedFile').html('<font size="2">' + f.path + '</font>');

        $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a>');
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .pickle file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };
});

ipcRenderer.on("curation-input", (event, data) => {
  $("#curation-input").text(data);
});
