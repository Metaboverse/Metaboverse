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
        $('#selectedTranscriptomics').replaceWith('<font size="2">' + f.path + '</font>');
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
        $('#selectedProteomics').replaceWith('<font size="2">' + f.path + '</font>');
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
        $('#selectedMetabolomics').replaceWith('<font size="2">' + f.path + '</font>');

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
  document.getElementById("metadata-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("metadata-input").value.split(".");

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
        $('#selectedMetadata').replaceWith('<font size="2">' + f.path + '</font>');

        update_session_info("metadata", f.path);
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
  document.getElementById("updateExperiment").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateExperiment").value;

    if (inputVal === "expType1") {
      experiment_type = null;
    } else if (inputVal === "expType2") {
      experiment_type = "timecourse";
    } else if (inputVal === "expType3") {
      experiment_type = "flux_balance";
    } else if (inputVal === "expType4") {
      experiment_type = "multiple_conditions";
    } else {}

    try {
      update_session_info("experiment", experiment_type);
    } catch (error) {
      console.log(error);
      alert(error);
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

      update_session_info("experiment", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "Experiment name is not valid."
      );
    }
  }
});
