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

const exec = require("child_process").exec;
var path = require("path");
var fs = require("fs");
var $ = require("jquery");

var userDataPath = app.getPath("userData");
var session_file = userDataPath + "/session_data.json";

// do something in case there are multiple spaces in the session data path
//if (" " in session_file) {
//  for
//}

var session_format = session_file.split(" ")[0]
  + String.fromCharCode(92)
  + " "
  + session_file.split(" ")[1]

var progressFile = "data/progress_log.json";
var scriptFilename = path.join(__dirname, "../python", "__main__.py");

fs.copyFile(
  "data/progress_log_template.json",
  "data/progress_log.json",
  err => {
    if (err) throw err;
    console.log("Progress log file was copied for this session");
  }
);

fs.watch(progressFile, function(event, filename) {
  var elem = document.getElementById("progressBar");
  var sum_values = 0;
  var session = JSON.parse(fs.readFileSync(progressFile).toString());
  for (j in session) {
    //loop through the array
    sum_values += session[j]; //Do the math!
  }
  elem.style.width = sum_values + "%";
  elem.innerHTML = sum_values + "%";

  if (sum_values >= 100) {
    sum_values = 0;
    displayOptions();
  }
});

var transcriptomics = false;
var proteomics = false;
var metabolomics = false;

var t_val = getDefault("transcriptomics");
var p_val = getDefault("proteomics");
var m_val = getDefault("metabolomics");

if (t_val !== null) {
  transcriptomics = true;
}

if (p_val !== null) {
  proteomics = true;
}

if (m_val !== null) {
  metabolomics = true;
}

// Run Python CLI
// Assemble parameters based on available session data
// Pass progress_log file for updates for each section of processing
function parseCommand(args_dict) {
  var commandString = "";
  for (key in args_dict) {
    if (args_dict[key] === false) {
    } else if (args_dict[key] === true) {
      commandString = commandString + " --" + key;
    } else {
      commandString = commandString + " --" + key + " " + args_dict[key];
    }
  }
  return commandString;
}

function callback() {
  console.log("finished");
}

function execute(command, callback) {
  exec(command, (error, stdout, stderr) => {
    callback(stdout);
    callback(stderr);
  });
}

runBuild = function(_callback) {

  if (get_session_info("processed") === true) {
    var elem = document.getElementById("progressBar");
    elem.style.width = "100%";
    elem.innerHTML = "100%";
    $("#content").replaceWith(
      '<a href="motif.html"><div id="continue"><font size="3">Run Motif Search</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>'
    );
  } else {
    curated = getArgument("curation_url")

    if (curated !== "None") {
      graphDictionary = {
        output: getArgument("output"),
        output_file: "find",
        species_id: "find",
        organism_curation: curated,
        metadata: getArgument("metadata"),
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        experiment: getArgument("experiment"),
        progress_log: path.resolve("data/progress_log.json"),
        session_data: session_format
      }
    } else {
      graphDictionary = {
        output: getArgument("output"),
        output_file: getArgument("database_url"),
        species_id: getArgument("organism_id"),
        model_file:
          getArgument("output") +
          getArgument("organism_id") +
          "_metaboverse_db.pickle",
        metadata: getArgument("metadata"),
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        experiment: getArgument("experiment"),
        progress_log: path.resolve("data/progress_log.json"),
        session_data: session_format
      }
    }

    var cmd = parseCommand(graphDictionary);
    console.log("Running: python " + scriptFilename + " curate " + cmd);
    execute("python " + scriptFilename + " curate " + cmd, output => {
      console.log(output);
      update_session_info("processed", true);
    });
    return _callback;
  }
};

function displayOptions() {
  if (
    (transcriptomics === true) |
    (proteomics === true) |
    (metabolomics === true)
  ) {
    $("#content").replaceWith(
      '<a href="motif.html"><div id="continue"><font size="3">Run Motif Search</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>'
    );
  } else {
    $("#content").replaceWith(
      '<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>'
    );
  }
}
