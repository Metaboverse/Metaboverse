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

const {
  ipcRenderer
} = require("electron");
var $ = require("jquery");

window.addEventListener("load", function(event) {

  let session_items = {

    "spacer4": "Experiment metadata:",

    "Experiment name": "experiment_name",
    "Experiment type": "experiment_type",
    "Sample labels": "labels",


    "spacer5": "Input data:",

    "Transcriptomics file name": "transcriptomics",
    "Proteomics file name": "proteomics",
    "Metabolomics file name": "metabolomics",


    "spacer1": "Database:",

    "Output location": "output",

    "small_gap0": "",

    "Database file name": "output_file",
    "Metaboverse-cli version used": "model_version",
    "Date generated": "model_date",


    "spacer2": "Curation:",

    "Organism name": "organism",
    "Organism ID": "organism_id",
    "Organism database source": "database_source",
    "Organism database version": "database_version",
    "Organism curation": "curation_url",
    "Metaboverse-cli organism curation version": "curation_version",
    "Organism curation date generated": "curation_date",


    "spacer3": "Other curation sources:",

    "Reaction neighbors database": "neighbors_url",
    "Metaboverse-cli reaction neighbors database version": "neighbors_version",
    "Reaction neighbors database date generated": "neighbors_date",

    "small_gap2": "",

    "Reaction network template": "template_url",
    "Metaboverse-cli reaction network template version": "template_version",
    "Reaction network template date generated": "template_date",


    "spacer6": "Other information:",

    "Reaction collapse used modifiers?": "collapseWithModifiers",
    "Gene expression broadcast to missing proteins?": "broadcastGeneExpression",
    "Metabolites broadcast to protein complexes?": "broadcastMetabolites",
    "Blocklisted nodes": "blocklist",
    "Reaction collapse threshold": "collapse_threshold"
  };

  let display = "";
  for (item in session_items) {
    if (item.includes("spacer")) {
      display = display +
        "<h4>" + session_items[item] + "</h4>";
    } else if (item.includes("small_gap")) {
      display = display + "<br>";
    } else {
      let display_item = getArgument(session_items[item]);
      if (display_item === undefined) {
        display_item = "Unable to find this information.";
      } else if (typeof display_item === 'number') {
        display_item = display_item;
      } else if (typeof display_item === 'object') {
        display_item = display_item[0];
      } else if (display_item === true) {
        display_item = "True";
      } else if (display_item === false) {
        display_item = "False";
      } else if (item.includes("Experiment type") && display_item === "None") {
        display_item = "Standard/Two-condition";
      } else if (display_item === "" || display_item === "None") {
        display_item = "Not provided";
      } else {
        display_item = display_item[0].toUpperCase() + display_item.substring(1);
      }
      display = display +
        "&#8226;&nbsp;&nbsp;&nbsp;&nbsp;" +
        item + ":&nbsp;" +
        "<font color='#00008b'>" + display_item + "</font>" +
        "<br />";
    }
  }
  document.getElementById("display-session").innerHTML = display;
})
