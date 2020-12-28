/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2020 Jordan A. Berg
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

const { ipcRenderer } = require("electron");
var $ = require("jquery");

window.addEventListener("load", function(event) {

  let session_items = {

    "spacer1": "Database:",

    "Database file name": "database_url",
    "Curation database file name": "curation_url",
    "Output location": "output",
    "Date generated": "database_date",

    "spacer2": "Curation:",

    "Organism": "organism",
    "Organism Reactome ID": "organism_id",
    "Date curated": "curation_date",
    "Database version": "database_version",

    "spacer3": "Experiment information:",

    "Experiment": "experiment_name",
    "Type": "experiment_type",
    "Sample labels": "labels",

    "spacer4": "User provided data information:",

    "Transcriptomics file name": "transcriptomics",
    "Proteomics file name": "proteomics",
    "Metabolomics file name": "metabolomics",

    "spacer5": "Other:",

    "Reaction collapse used modifiers?": "collapseWithModifiers",
    "Gene expression broadcast to missing proteins?": "broadcastGeneExpression",
    "Metabolites broadcast to protein complexes?": "broadcastMetabolites",
    "Blocklisted nodes": "blocklist"
  };

  let display = "";
  for (item in session_items) {
    if (item.includes("spacer")) {
      display = display
        + "<h4>" + session_items[item] + "</h4>";
    } else {
      let display_item = getArgument(session_items[item]);
      if (display_item === undefined) {
        display_item = "Cannot find this information from the file. Sorry!";
      } else if (typeof display_item === 'number') {
        display_item = display_item;
      } else if (typeof display_item === 'object') {
        display_item = display_item[0];
      } else if (display_item === true) {
        display_item = "True";
      } else if (display_item === false) {
        display_item = "False";
      } else if (display_item === "" || display_item === "None") {
        display_item = "N/A";
      } else {
        display_item = display_item[0].toUpperCase() + display_item.substring(1);
      }
      display = display
        + "&#8226;&nbsp;&nbsp;&nbsp;&nbsp;"
        + item + ":&nbsp;"
        + "<font color='#00008b'>" + display_item + "</font>"
        + "<br />";
    }
  }
  document.getElementById("display-session").innerHTML = display;
})
