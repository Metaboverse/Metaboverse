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

var app = require("electron").remote.app;
var basePath = app.getAppPath();

var userDataPath = app.getPath("userData");
var session_file = userDataPath + "/session_data.json";
var progress_file = userDataPath + "/progress_log.json"

// do something in case there are multiple spaces in the session data path
//if (" " in session_file) {
//  for
//}
var replacer = String(String.fromCharCode(92) + " ");
var session_format = session_file.replace(" ", replacer)
var progress_format = progress_file.replace(" ", replacer)

var scriptFilename;
if (navigator.appVersion.indexOf("Win") != -1) {
  scriptFilename = path.join(__dirname, "../python", "metaboverse-win");
} else if (navigator.appVersion.indexOf("Mac") != -1) {
  scriptFilename = path.join(__dirname, "../python", "metaboverse-mac");
} else if (navigator.appVersion.indexOf("X11") != -1) {
  scriptFilename = path.join(__dirname, "../python", "metaboverse-mac");
} else if (navigator.appVersion.indexOf("Linux") != -1) {
  scriptFilename = path.join(__dirname, "../python", "metaboverse-linux");
} else {
  scriptFilename = path.join(__dirname, "../python", "__main__.py");
}

fs.copyFile(
  basePath + "/data/progress_log_template.json",
  progress_file,
  err => {
    if (err) throw err;
    console.log("Progress log file was copied for this session");
  }
);

fs.watch(progress_file, function(event, filename) {
  var elem = document.getElementById("progressBar");
  var sum_values = 0;
  var session = JSON.parse(fs.readFileSync(progress_file).toString());
  for (j in session) {
    //loop through the array
    if (sum_values += session[j] < 100) {
      sum_values += session[j]; //Do the math!
    }
    if (sum_values >= 100) {
      sum_values = 100;
      break;
    }
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

function execute(command, callback) {
  exec(command, (error, stdout, stderr) => {
    callback(stdout);
    if (stdout.toLowerCase().includes("exception")
        | stdout.toLowerCase().includes("error")
        | stderr.toLowerCase().includes("exception")
        | stderr.toLowerCase().includes("error")) {
      alert("Metaboverse Build encountered an error             \n\n" + stderr);
      callback(stderr);
    }
  });
}

runBuild = function(_callback) {

  if (get_session_info("processed") === true) {
    var elem = document.getElementById("progressBar");
    elem.style.width = "100%";
    elem.innerHTML = "100%";
    displayOptions();
  } else {
    curated = getArgument("curation_url")
    let labels = getArgument("labels");
    labels = labels.replace(/\s/g,'');
    let blacklist = getArgument("blacklist");
    blacklist = blacklist.replace(/\s/g,'');
    if (String(curated) !== "None") {
      graphDictionary = {
        output: getArgument("output"),
        output_file: "find",
        species_id: "find",
        organism_curation: curated,
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blacklist: blacklist,
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        progress_log: progress_format,
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
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blacklist: blacklist,
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        progress_log: progress_format,
        session_data: session_format
      }
    }

    var cmd = parseCommand(graphDictionary);
    console.log("Running: " + scriptFilename + " curate " + cmd);
    execute(scriptFilename + " curate " + cmd, output => {
      console.log(output);
      update_session_info("processed", true);
    });
    return _callback;
  }
};

function displayOptions() {
  update_session_info("current_pathway", null);
  database_url = get_session_info("database_url");
  var data = JSON.parse(fs.readFileSync(database_url).toString());
  update_session_info("max_value", data.metadata.max_value);
  update_session_info("max_stat", data.metadata.max_stat);
  update_session_info("database_date", data.metadata.database_date);
  update_session_info("curation_date", data.metadata.curation_date);
  update_session_info("reactome_version", data.metadata.reactome_version);
  update_session_info("blacklist", data.metadata.blacklist);
  if (
    (transcriptomics === true) |
    (proteomics === true) |
    (metabolomics === true) |
    (get_session_info("provided_data") === true) |
    (get_session_info("provided_data") === "true") |
    (data.categories.length > 0)
  ) {
    $("#content").html(
      '<a href="motif.html"><div id="continue"><font size="3">View Motif Search</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a>'
    );
    update_session_info("provided_data", true);
  } else {
    $("#content").html(
      '<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a>'
    );
  }
}
