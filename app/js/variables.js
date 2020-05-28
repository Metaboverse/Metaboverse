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

var fs = require("fs");
var $ = require("jquery");

var coll_mos = false;
var broadcast_gene = true;

function collapseWithModifiers() {
  if (coll_mos === false) {
    coll_mos = true;
    update_session_info("collapseWithModifiers", true);
  } else {
    coll_mos = false;
    update_session_info("collapseWithModifiers", false);
  }
  console.log("Reaction collapse evaluation with modifiers: ", coll_mos)
}

function broadcastGeneExpression() {
  if (broadcast_gene === false) {
    broadcast_gene = true;
    update_session_info("broadcastGeneExpression", true);
  } else {
    broadcast_gene = false;
    update_session_info("broadcastGeneExpression", false);
  }
  console.log("Broadcast gene expression values: ", broadcast_gene)
}

window.addEventListener("load", function(event) {
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
});

window.addEventListener("load", function(event) {
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
});

window.addEventListener("load", function(event) {
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
});

window.addEventListener("load", function(event) {

  update_session_info("labels", "none");
  document.getElementById("updateExperiment").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var experiment_type = document.getElementById("updateExperiment").value;
    if (experiment_type === "null") {
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
        "<form>"
        + "Sample labels: "
        + "<button class='info' title='Enter the names for each condition or timepoint for you dataset in the order that they appear in the data table. Labels should be separated by a comma.'><i>i</i></button>"
        + "<br />"
        + "<br />"
        + "<input type='text' class='experimentName' id='updateExperimentLabels'></input>"
        + "</form>"
        + "<br />"
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
});

window.addEventListener("load", function(event) {
  document.getElementById("updateExperimentName").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateExperimentName").value.split(".");

    try {
      console.log("Your provided experiment name: ", inputVal);

      update_session_info("experiment_name", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "Experiment name is not valid."
      );
    }
  }
});

window.addEventListener("load", function(event) {

  var inputVal = document.getElementById("updateBlacklist").value;
  update_session_info("blacklist", inputVal);

  document.getElementById("updateBlacklist").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateBlacklist").value;

    try {
      console.log("Your provided blacklisted entities: ", inputVal);

      update_session_info("blacklist", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "IDs are not valid."
      );
    }
  }


});
