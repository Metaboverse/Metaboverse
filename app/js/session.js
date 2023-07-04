/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) 2022 Metaboverse

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

var { ipcRenderer } = require("electron");
var $ = require("jquery");


let session_items = {

  "Displaying session data for file": "database_url",

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


async function load_session_data() {  
  let display = "";
  for (item in session_items) {
    if (item === "Displaying session data for file") {
      let display_item = await get_argument_async(session_items[item]);
      display += "<center>" + 
          "<b>" + item + "</b>: " + 
          "<font color='#00008b'>" + display_item + "</font>" +
        "</center>" + 
        "<br />";
    } else if (item.includes("spacer")) {
      display +=
        "<h4>" + session_items[item] + "</h4>";
    } else if (item.includes("small_gap")) {
      display += "<br>";
    } else {
      let display_item = await get_argument_async(session_items[item]);
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
      display_item = display_item.toString().replace(/\\\\ /g, ' ').replace(/\\\\/g, '\\');
      display +=
        "&#8226;&nbsp;&nbsp;&nbsp;&nbsp;" +
        item + ":&nbsp;" +
        "<font color='#00008b'>" + display_item + "</font>" +
        "<br />";
    }
  }
  document.getElementById("display-session").innerHTML = display;
}


window.addEventListener("load", function(event) {
  load_session_data();
})
