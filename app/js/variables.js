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

var fs = require("fs");
var $ = require("jquery");

window.addEventListener("load", function(event) {
  event.preventDefault();
  event.stopPropagation();

  // If user goes back to this page, force re-curation
  update_session_info("processed", false);

  // Transcriptomics
  var transcriptomicsURL = get_session_info("transcriptomics");
  var defaultTranscriptomics = "No file selected";
  if (transcriptomicsURL !== null) {
    defaultTranscriptomics = transcriptomicsURL;
  }
  $('#selectedTranscriptomics').append('<font size="2">' + defaultTranscriptomics + '</font>');

  // Proteomics
  var proteomicsURL = get_session_info("proteomics");
  var defaultProteomics = "No file selected";
  if (proteomicsURL !== null) {
    defaultProteomics = proteomicsURL;
  }
  $('#selectedProteomics').append('<font size="2">' + defaultProteomics + '</font>');

  // Metabolomics
  var metabolomicsURL = get_session_info("metabolomics");
  var defaultMetabolomics = "No file selected";
  if (metabolomicsURL !== null) {
    defaultMetabolomics = metabolomicsURL;
  }
  $('#selectedMetabolomics').append('<font size="2">' + defaultMetabolomics + '</font>');

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session()
  }

  console.log(get_session_info("collapseWithModifiers"))
  if (get_session_info("collapseWithModifiers") === true) {
    console.log('mod')
    document.getElementById("use_modifiers_in_collapse").checked = true;
  } else {
    console.log('mod_')
    document.getElementById("use_modifiers_in_collapse").checked = false;
  }
  document.getElementById("use_modifiers_in_collapse").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("collapseWithModifiers") === false) {
      update_session_info("collapseWithModifiers", true);
    } else {
      update_session_info("collapseWithModifiers", false);
    }
    console.log("Reaction collapse evaluation with modifiers: ", get_session_info("collapseWithModifiers"))
  }

  console.log(get_session_info("broadcastGeneExpression"))
  if (get_session_info("broadcastGeneExpression") === true) {
    console.log('broad')
    document.getElementById("broadcast_gene_expression").checked = true;
  } else {
    console.log('broad_')
    document.getElementById("broadcast_gene_expression").checked = false;
  }
  document.getElementById("broadcast_gene_expression").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("broadcastGeneExpression") === false) {
      update_session_info("broadcastGeneExpression", true);
    } else {
      update_session_info("broadcastGeneExpression", false);
    }
    console.log("Broadcast gene expression values: ", get_session_info("broadcastGeneExpression"))
  }

  console.log(get_session_info("broadcastMetabolites"))
  if (get_session_info("broadcastMetabolites") === true) {
    console.log('broadM')
    document.getElementById("broadcast_metabolite_expression").checked = true;
  } else {
    console.log('broadM_')
    document.getElementById("broadcast_metabolite_expression").checked = false;
  }
  document.getElementById("broadcast_metabolite_expression").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("broadcastMetabolites") === false) {
      update_session_info("broadcastMetabolites", true);
    } else {
      update_session_info("broadcastMetabolites", false);
    }
    console.log("Broadcast metabolite values: ", get_session_info("broadcastMetabolites"))
  }

  document.getElementById("transcriptomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('transcriptomics-input').click();
  }

  document.getElementById("transcriptomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document
      .getElementById("transcriptomics-input")
      .value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedTranscriptomics').html('<font size="2">' + f.path + '</font>');
        update_session_info("transcriptomics", f.path);

        transcriptomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  document.getElementById("proteomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('proteomics-input').click();
  }

  document.getElementById("proteomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("proteomics-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedProteomics').html('<font size="2">' + f.path + '</font>');
        update_session_info("proteomics", f.path);

        proteomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  document.getElementById("metabolomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('metabolomics-input').click();
  }

  document.getElementById("metabolomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("metabolomics-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedMetabolomics').html('<font size="2">' + f.path + '</font>');

        update_session_info("metabolomics", f.path);

        metabolomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };
  /*
  document.getElementById("reactions-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("reactions-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedReactions').html('<font size="2">' + f.path + '</font>');

        update_session_info("additional_reactions", f.path);

        reactions = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  }
  */

  if (get_session_info("experiment_name") !== null) {
    $('#updateExperimentName').val(get_session_info("experiment_name"));
  }

  if (get_session_info("experiment_type") !== null) {
    $('#updateExperiment').val(get_session_info("experiment_type"));

    if ((get_session_info("experiment_type") === "timecourse") | (get_session_info("experiment_type") === "multiple_conditions")) {
      $("#nameField").html(
        "<form>" +
        "Sample labels: " +
        "<button class='info' title='Enter the names for each condition or timepoint for you dataset in the order that they appear in the data table. Labels should be separated by a comma.'><i>i</i></button>" +
        "<br />" +
        "<br />" +
        "<input type='text' class='experimentName' id='updateExperimentLabels'></input>" +
        "</form>" +
        "<br />"
      );
      $('#updateExperimentLabels').val(get_session_info("labels"));

      document.getElementById("updateExperimentLabels").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();

        var inputVal = document.getElementById("updateExperimentLabels").value;

        try {
          console.log("Your provided labels: ", inputVal);

          update_session_info("labels", inputVal);
        } catch (error) {
          console.log(error);
          alert(
            "Labels are not valid."
          );
        }
      }

    } else {
      $("#nameField").html('');
      update_session_info("labels", "0");
    }
  }





  document.getElementById("updateExperiment").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var experiment_type = document.getElementById("updateExperiment").value;
    if (experiment_type === "null" || experiment_type === null) {
      experiment_type = null;
    } else {}

    try {
      update_session_info("experiment_type", experiment_type);
    } catch (error) {
      console.log(error);
      alert(error);
    }

    // If timecourse or multiple conditions, have user input labels in correct order as listed  in dataframe
    if ((experiment_type === "timecourse") | (experiment_type === "multiple_conditions")) {
      $("#nameField").html(
        "<form>" +
        "Sample labels: " +
        "<button class='info' title='Enter the names for each condition or timepoint for you dataset in the order that they appear in the data table. Labels should be separated by a comma.'><i>i</i></button>" +
        "<br />" +
        "<br />" +
        "<input type='text' class='experimentName' id='updateExperimentLabels'></input>" +
        "</form>" +
        "<br />"
      );

      document.getElementById("updateExperimentLabels").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();

        var inputVal = document.getElementById("updateExperimentLabels").value;

        try {
          console.log("Your provided labels: ", inputVal);

          update_session_info("labels", inputVal);
        } catch (error) {
          console.log(error);
          alert(
            "Labels are not valid."
          );
        }
      }

    } else {
      $("#nameField").html('');
      update_session_info("labels", "0");
    }
  };

  document.getElementById("updateExperimentName").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateExperimentName").value
    //.split(".");

    try {
      var replacer = String(String.fromCharCode(92) + " ");
      var inputVal = inputVal.replace(/ /g, replacer)
      console.log("Your provided experiment name: ", inputVal);
      update_session_info("experiment_name", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "Experiment name is not valid."
      );
    }
  }

  $('#updateBlocklist').val(get_session_info("blocklist"));

  document.getElementById("updateBlocklist").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateBlocklist").value;

    try {
      console.log("Your provided blocklisted entities: ", inputVal);

      update_session_info("blocklist", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "IDs are not valid."
      );
    }
  }
});
