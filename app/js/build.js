/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
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
const exec = require("child_process").exec;
var path = require("path");
var fs = require("fs");
var $ = require("jquery");

var app = require("electron").remote.app;
var basePath = app.getAppPath();

var userDataPath = app.getPath("userData");
var session_file = userDataPath + path.sep + "session_data.json";
var progress_file = userDataPath + path.sep + "progress_log.json"

var scriptFilename;
if (navigator.appVersion.indexOf("Win") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-win.exe");
} else if (navigator.appVersion.indexOf("Mac") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-mac");
} else if (navigator.appVersion.indexOf("X11") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-mac");
} else if (navigator.appVersion.indexOf("Linux") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-linux");
} else {
  console.log("Unable to locate metaboverse-cli binary")
}

fs.copyFile(
  basePath + path.sep + "data" + path.sep + "progress_log_template.json",
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
    let blocklist = getArgument("blocklist");
    if (blocklist === "") {
      blocklist = "no_blocklist"
    } else {
      blocklist = blocklist.replace(/\s/g,'');
    }
    if (String(curated) !== "None") {
      let db_url;
      if (String(getArgument("database_url")) === "None") {
        db_url = "find";
      } else {
        db_url = getArgument("database_url");
      }
      graphDictionary = {
        output: getArgument("output"),
        output_file: db_url,
        species_id: "find",
        organism_curation: curated,
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        progress_log: "\"" + progress_file + "\"",
        session_data: "\"" + session_file + "\""
      }
    } else {
      graphDictionary = {
        output: getArgument("output"),
        output_file: getArgument("database_url"),
        species_id: getArgument("organism_id"),
        model_file:
          getArgument("output").slice(0,-1) +
          path.sep + 
          getArgument("organism_id") +
          "_metaboverse_db.pickle\"",
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        progress_log: "\"" + progress_file + "\"",
        session_data: "\"" + session_file + "\""
      }
    }

    var cmd = parseCommand(graphDictionary);
    console.log("Running: " + scriptFilename + " curate " + cmd);
    execute(scriptFilename + " curate " + cmd, output => {
      update_session_info("processed", true);
    });
    return _callback;
  }
};

function displayOptions() {
  update_session_info("current_pathway", null);
  database_url = get_session_info("database_url");
  try {
    var data = JSON.parse(fs.readFileSync(database_url).toString());
  } catch(e) {
    alert('Failed to open database URL: \n' + database_url)
    var elem = document.getElementById("progressBar");
    elem.style.width = "0%";
    elem.innerHTML = "0%";
  }
  update_session_info("max_value", data.metadata.max_value);
  update_session_info("max_stat", data.metadata.max_stat);
  update_session_info("database_date", data.metadata.database_date);
  update_session_info("curation_date", data.metadata.curation_date);
  update_session_info("reactome_version", data.metadata.reactome_version);
  update_session_info("blocklist", data.metadata.blocklist);
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

window.addEventListener("load", function(event) {
  $('#content').append('')
  runBuild()

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session()
  }

})
