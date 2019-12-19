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

window.addEventListener('load', function(event) {

  document.getElementById("dropTranscriptomics").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropTranscriptomics").value.split(".")

  if (inputVal[inputVal.length-1] !== "tsv" & inputVal[inputVal.length-1] !== "txt") {
    alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("transcriptomics", path)

      transcriptomics = true;

    } catch (error) {
      console.log(error)
      alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work.')
    }
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("dropProteomics").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropProteomics").value.split(".")

  if (inputVal[inputVal.length-1] !== "tsv" & inputVal[inputVal.length-1] !== "txt") {
    alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("proteomics", path)

      proteomics = true;

    } catch (error) {
      console.log(error)
      alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work.')
    }
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("dropMetabolomics").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropMetabolomics").value.split(".")

  if (inputVal[inputVal.length-1] !== "tsv" & inputVal[inputVal.length-1] !== "txt") {
    alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("metabolomics", path)

      metabolomics = true;

    } catch (error) {
      console.log(error)
      alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work.')
    }
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("dropMetadata").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropMetadata").value.split(".")

  if (inputVal[inputVal.length-1] !== "tsv" & inputVal[inputVal.length-1] !== "txt") {
    alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("metadata", path)

    } catch (error) {
      console.log(error)
      alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work.')
    }
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("dropBlacklist").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropBlacklist").value.split(".")

  if (inputVal[inputVal.length-1] !== "tsv" & inputVal[inputVal.length-1] !== "txt") {
    alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("blacklist", path)

    } catch (error) {
      console.log(error)
      alert('Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work.')
    }
  }

  }
})

var collapse = false;
window.addEventListener('load', function(event) {

  document.getElementById("updateCollapse").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()
  if (collapse === true) {
    collapse = false;
  } else {
    collapse = true;
  }

  try {

    update_session_info("collapse_missing_reactions", collapse)

  } catch (error) {
    console.log(error)
    alert(error)
  }

  }
})

var split = false;
window.addEventListener('load', function(event) {

  document.getElementById("updateSplit").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()
  if (split === true) {
    split = false;
  } else {
    split = true;
  }

  try {

    update_session_info("split_duplicate_nodes", split)

  } catch (error) {
    console.log(error)
    alert(error)
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("updateExperiment").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("updateExperiment").value

  if (inputVal === "expType1") {
    experiment_type = null;
  } else if (inputVal === "expType2") {
    experiment_type = "timecourse";
  } else if (inputVal === "expType3") {
    experiment_type = "flux_balance";
  } else if (inputVal === "expType4") {
    experiment_type = "multiple_conditions";
  } else {

  }

  try {

    update_session_info("experiment", experiment_type)

  } catch (error) {
    console.log(error)
    alert(error)
  }

  }
})

window.addEventListener('load', function(event) {

  document.getElementById("updateNormalization").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("updateNormalization").value

  if (inputVal === "normType1") {
    normalization_type = null;
  } else if (inputVal === "normType2") {
    normalization_type = "log";
  } else if (inputVal === "normType3") {
    normalization_type = "z-score";
  } else if (inputVal === "normType4") {

  } else {

  }

  try {

    update_session_info("normalization", normalization_type)

  } catch (error) {
    console.log(error)
    alert(error)
  }

  }
})
