/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2021 Jordan A. Berg
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

var { ipcRenderer, ipcMain, remote } = require("electron");
var { dialog } = require("electron").remote;
var path = require("path");
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
  if (get_session_info("organism") !== null) {
    $('#speciesMenu').val(get_session_info("organism"));
  }
});



// Change user selection based on input
window.addEventListener("load", function(event) {
  event.preventDefault();
  event.stopPropagation();

  document.getElementById("speciesMenu").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var selection = document.getElementById("speciesMenu").value;
    console.log("User selected:", abbreviation_dict[selection]);
    update_session_info("organism", selection, (abbrev_dict = abbreviation_dict));
    species_change = true;
    check_changes();
  }

  // check if inputs already there
  if (get_session_info("organism") !== null && get_session_info("database_url") !== null) {
    $("#content").html(
      '<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
  } else if (get_session_info("curation_url") !== null || get_session_info("template_url") !== null && get_session_info("database_source") !== "reactome") {
    $("#content").html(
      '<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
  } else if (get_session_info("curation_url") !== null && get_session_info("database_source") !== "biomodels/bigg") {
    $("#content").html(
      '<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
  } else {
    $('#content').append('<a href=""></a>')
  }

  // Get curation reference file from user
  var curURL = get_session_info("curation_url");
  var defaultCuration = "No file selected";
  if (curURL !== null) {
    if (curURL.split('.').pop().toLowerCase() === 'mvdb') {
      defaultCuration = curURL;
      $("#content").html(
        '<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
      window.scrollTo(0,document.body.scrollHeight);
    }
  }
  $('#selectedFile').append('<font size="2">' + defaultCuration + '</font>');

  // Get network template file from user
  var templateURL = get_session_info("template_url");
  var defaultTemplate = "No file selected";
  if (templateURL !== null) {
    if (templateURL.split('.').pop().toLowerCase() === 'mvrs') {
      defaultTemplate = templateURL;
    }
  }
  $('#selectedTemplate').append('<font size="2">' + defaultTemplate + '</font>');

  // Get reaction neighbors dictionary from user
  var neighborsURL = get_session_info("neighbors_url");
  var defaultNeighbors = "No file selected";
  if (neighborsURL !== null) {
    if (neighborsURL.split('.').pop().toLowerCase() === 'mvrs') {
      defaultNeighbors = neighborsURL;
    }
  }
  $('#selectedNeighbors').append('<font size="2">' + defaultNeighbors + '</font>');

  var defaultSBML = "No file selected";
  if (curURL !== null) {
    if (curURL.split('.').pop().toLowerCase() === 'sbml' || curURL.split('.').pop().toLowerCase() === 'xml') {
      defaultSBML = curURL;
      $("#content").html(
        '<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
      window.scrollTo(0,document.body.scrollHeight);
    }
  }
  $('#selectedSBML').append('<font size="2">' + defaultSBML + '</font>');

  var outputURL = get_session_info("database_url");
  var defaultOutput = "No output selected";
  if (outputURL !== null) {
    defaultOutput = outputURL;
  }
  $('#selectedOutput').append('<font size="2">' + defaultOutput + '</font>');

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session()
  }

  document.getElementById("dropDatabaseOutput").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('output-input').click();
  }

  document.getElementById("dropDatabaseCuration").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('curation-input').click();
  }

  document.getElementById("dropGraphTemplate").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('template-input').click();
  }

  document.getElementById("dropNeighborDictionary").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('neighbors-input').click();
  }

  document.getElementById("dropDatabaseSBML").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('sbml-input').click();
  }

  // Select output directory from pop-out menu
  document.getElementById("output-input").onclick = function(event) {
    filename = dialog
      .showSaveDialog({
        defaultPath: ".." + path.sep + ".." + path.sep,
        properties: ["createDirectory"],
        filters: [
          { name: "Metaboverse-formatted database (*.mvrs)", extensions: ["mvrs"] }
        ]
      })
      .then(result => {
        let hasExtension = /\.[^\/\\]+$/.test(result.filePath);
        if (hasExtension === false) {
          result.filePath = `${ result.filePath }.${ "mvrs" }`;
        }
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

  // Drop pre-existing metabolic network curation for further analysis
  document.getElementById("curation-input").onchange = function(event) {
    console.log(document.getElementById("curation-input"))
    event.preventDefault();
    event.stopPropagation();

    if (document.getElementById("curation-input").value.split('.').pop().toLowerCase() !== "mvdb") {
      alert(
        "Input is not a .mvdb file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        update_session_info("curation_url", f.path);
        update_session_info("database_source", "reactome");
        $('#selectedFile').html('<font size="2">' + f.path + '</font>');
        $('#selectedSBML').html('<font size="2">No file selected</font>');
        $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
        document.getElementById("curation-input").value = '';
        window.scrollTo(0,document.body.scrollHeight);
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .mvdb file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  // Drop pre-existing organism reaction network template for further analysis
  document.getElementById("template-input").onchange = function(event) {
    console.log(document.getElementById("template-input"))
    event.preventDefault();
    event.stopPropagation();

    if (document.getElementById("template-input").value.split('.').pop().toLowerCase() !== "mvrs") {
      alert(
        "Input is not a .mvrs file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        update_session_info("template_url", f.path);
        update_session_info("database_source", "reactome");
        $('#selectedTemplate').html('<font size="2">' + f.path + '</font>');
        $('#selectedSBML').html('<font size="2">No file selected</font>');
        $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
        document.getElementById("template-input").value = '';
        window.scrollTo(0,document.body.scrollHeight);
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .mvrs file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  // Drop pre-existing reaction neighbors dictionary for further analysis
  document.getElementById("neighbors-input").onchange = function(event) {
    console.log(document.getElementById("neighbors-input"))
    event.preventDefault();
    event.stopPropagation();

    if (document.getElementById("neighbors-input").value.split('.').pop().toLowerCase() !== "nbdb") {
      alert(
        "Input is not a .nbdb file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        update_session_info("neighbors_url", f.path);
        update_session_info("database_source", "reactome");
        $('#selectedNeighbors').html('<font size="2">' + f.path + '</font>');
        $('#selectedSBML').html('<font size="2">No file selected</font>');
        document.getElementById("neighbors-input").value = '';
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .nbdb file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  // Drop pre-existing BiGG or BioModels network curation for further analysis
  document.getElementById("sbml-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("sbml-input").value.split(".");

    if (document.getElementById("sbml-input").value.split('.').pop().toLowerCase() !== "xml" && document.getElementById("sbml-input").value.split('.').pop().toLowerCase() !== "sbml") {
      alert(
        "Input is not a .xml or .sbml file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        update_session_info("curation_url", f.path);
        update_session_info("database_source", "biomodels/bigg");
        $('#selectedSBML').html('<font size="2">' + f.path + '</font>');
        $('#selectedFile').html('<font size="2">No file selected</font>');
        $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
        document.getElementById("sbml-input").value = '';
        window.scrollTo(0,document.body.scrollHeight);
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .xml or .sbml file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  if (get_session_info("forceNewCuration") === true) {
    document.getElementById("force_new_curation").checked = true;
  } else {
    document.getElementById("force_new_curation").checked = false;
  }
  document.getElementById("force_new_curation").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("forceNewCuration") === false) {
      update_session_info("forceNewCuration", true);
    } else {
      update_session_info("forceNewCuration", false);
    }
    console.log("Force fresh curation: ", get_session_info("forceNewCuration"))
  }




});

var output_change = false;
var species_change = false;

function check_changes() {
  if ((output_change === true) & (species_change === true)) {
    $('#content').html('<a href="variables.html"><div id="continue"><font size="3">Continue</font></div></a><br><br>');
    window.scrollTo(0,document.body.scrollHeight);
  }
}
