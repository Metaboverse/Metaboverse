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

const { ipcRenderer } = require("electron");
var $ = require("jquery");

function loadSessionData() {

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
    "Reactome version": "reactome_version",

    "spacer3": "Experiment information:",

    "Experiment": "experiment_name",
    "Type": "experiment_type",
    "Sample labels": "labels",

    "spacer4": "User provided data information:",

    "Transcriptomics file name": "transcriptomics",
    "Proteomics file name": "proteomics",
    "Metabolomics file name": "metabolomics",
    "Maximum absolute expression/abundance value": "max_value",
    "Maximum absolute statistical value": "max_stat",

    "spacer5": "Other:",

    "Reaction collapse used modifiers?": "collapseWithModifiers",
    "Gene expression broadcast to missing proteins?": "broadcastGeneExpression",
    "Blacklisted nodes": "blacklist"
  };

  let display = "";
  for (item in session_items) {
    console.log(item)
    if (item.includes("spacer")) {
      display = display
        + "<h4>" + session_items[item] + "</h4>";
    } else {
      let display_item = getArgument(session_items[item]);
      if (display_item === undefined) {
        display_item = "Not provided";
      } else if (typeof display_item == 'number') {
        display_item = display_item;
      } else if (typeof display_item == 'object') {
        display_item = display_item[0];
      } else if (display_item == true) {
        display_item = "True";
      } else if (display_item == false) {
        display_item = "False";
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
  document.getElementById("session_data").innerHTML = display;
}
