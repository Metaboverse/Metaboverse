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
const exec = require("child_process").exec;
var path = require("path");
var fs = require("fs");
var $ = require("jquery");

var session_file = path.join(__dirname, "..", "data", "session_data.json");
var progress_template_file = path.join(__dirname, "..", "data", "progress_log_template.json");
var progress_file = path.join(__dirname, "..", "data", "progress_log.json");

var timer = 5000;

var scriptFilename;

console.log("Operating System information:")
console.log(navigator.appVersion)

if (navigator.appVersion.indexOf("Win") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-windows.exe");
} else if (navigator.appVersion.indexOf("Mac") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-darwin");
} else if (navigator.appVersion.indexOf("Linux") != -1) {
  scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-linux");
} else {
  console.log("Unable to locate metaboverse-cli binary")
}

fs.copyFile(
  progress_template_file,
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
    sum_values += session[j];
  }
  if (sum_values >= 100) {
    sum_values = 100;
  }
  elem.style.width = sum_values + 5 + "%";
  elem.innerHTML = sum_values + "%";

  if (sum_values >= 100) {
    setTimeout(function(){
      displayOptions();
    }, timer);
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
    if (args_dict[key] === false || args_dict[key] === "false") {
    } else if (args_dict[key] === true || args_dict[key] === "true") {
      commandString = commandString + " --" + key;
    } else {
      commandString = commandString
                    + " --" + key + " \""
                    + decodeURIComponent(
                        encodeURIComponent(
                          args_dict[key]))
                        .replace(/"/g, '')
                        .replace(/\\/g, '\\\\')
                    + "\"";
    }
  }
  return commandString;
}

function write_log(command, stdout, stderr) {

  let today = new Date().toISOString().slice(0, 10);
  let experiment_name = getArgument("experiment_name").replace(/\s+/g, '-');

  fs.writeFileSync(
    path.join(getArgument("output").replace(/\"/g, ''), "metaboverse_session_" + experiment_name + ".log"),
    "Operating System information:\n" +
    navigator.appVersion + "\n" +
    "Log date: " + today + "\n\n" +
    command + "\n" +
    stdout + "\n" +
    "########\nSTDERR:\n########\n" +
    stderr,
    function(err) {
      if (err) throw err;
      console.log("Session log written to file");
    });

}


function execute(command, callback) {
  exec(command, (error, stdout, stderr) => {
    callback(stdout);

    if (stdout.toLowerCase().includes("exception") |
      stdout.toLowerCase().includes("error") |
      stderr.toLowerCase().includes("exception") |
      stderr.toLowerCase().includes("error")) {
      
      write_log(command, stdout, stderr);
      alert("Metaboverse Build encountered an error -- please check metaboverse_session.log in your output folder for detailed information.");
      callback(stderr);
    } else {
      write_log(command, stdout, stderr);
    }
  });
}

runBuild = function(_callback) {

  if (get_session_info("processed") === true) {
    var elem = document.getElementById("progressBar");
    elem.style.width = "100%";
    elem.innerHTML = "100%";
    setTimeout(function(){
      displayOptions();
    }, timer);
  } else {
    curated = getArgument("curation_url")
    let labels = getArgument("labels");
    if (labels === "") {
      labels = "0"
    } else {
      labels = labels.replace(/\s/g, '');
    }
    let blocklist = getArgument("blocklist");
    if (blocklist === "") {
      blocklist = "no_blocklist"
    } else {
      blocklist = blocklist.replace(/\s/g, '');
    }
    if (String(curated) !== "None") {
      let db_url;
      if (String(getArgument("database_url")) === "None") {
        db_url = "\"find\"";
      } else {
        db_url = String(decodeURIComponent(getArgument("database_url")));
      }
      graphDictionary = {
        output: getArgument("output"),
        organism_id: "\"find\"",
        output_file: db_url,
        organism_curation_file: curated,
        neighbor_dictionary_file: getArgument("neighbors_url"),
        graph_template_file: getArgument("template_url"),
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        database_source: getArgument("database_source"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        force_new_curation: getArgument("forceNewCuration"),
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        broadcast_metabolites: getArgument("broadcastMetabolites"),
        collapse_threshold: parseFloat(getArgument("collapse_threshold")) / 100,
        progress_log: "\"" + progress_file + "\"",
        session_data: "\"" + session_file + "\""
      }
    } else {
      graphDictionary = {
        output: getArgument("output"),
        organism_id: getArgument("organism_id"),
        output_file: getArgument("database_url"),
        organism_curation_file: getArgument("curation_url"),
        neighbor_dictionary_file: getArgument("neighbors_url"),
        graph_template_file: getArgument("template_url"),
        transcriptomics: getArgument("transcriptomics"),
        proteomics: getArgument("proteomics"),
        metabolomics: getArgument("metabolomics"),
        database_source: getArgument("database_source"),
        //additional_reactions: getArgument("additional_reactions"),
        experiment_type: getArgument("experiment_type"),
        experiment_name: getArgument("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        force_new_curation: getArgument("forceNewCuration"),
        collapse_with_modifiers: getArgument("collapseWithModifiers"),
        broadcast_genes: getArgument("broadcastGeneExpression"),
        broadcast_metabolites: getArgument("broadcastMetabolites"),
        collapse_threshold: parseFloat(getArgument("collapse_threshold")) / 100,
        progress_log: "\"" + progress_file + "\"",
        session_data: "\"" + session_file + "\""
      }
    }
    var cmd = parseCommand(graphDictionary);
    console.log("Running: " + scriptFilename + " curate" + cmd);
    execute(scriptFilename + " curate" + cmd, output => {
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
    update_session_info("curation_url", data.metadata.organism_curation);
    update_session_info("output", data.metadata.output);
    update_session_info("output_file", data.metadata.output_file);
    update_session_info("organism_id", data.metadata.organism_id);
    update_session_info("transcriptomics", data.metadata.transcriptomics);
    update_session_info("proteomics", data.metadata.proteomics);
    update_session_info("metabolomics", data.metadata.metabolomics);
    update_session_info("experiment_name", data.metadata.experiment_name);
    update_session_info("experiment_type", data.metadata.experiment_type);
    update_session_info("max_value", data.metadata.max_value);
    update_session_info("max_stat", data.metadata.max_stat);
    update_session_info("processed", true);
    update_session_info("collapseWithModifiers",
      data.metadata.collapse_with_modifiers);
    update_session_info("broadcastGeneExpression",
      data.metadata.broadcast_genes);
    update_session_info("broadcastMetabolites",
      data.metadata.broadcast_metabolites);
    update_session_info("labels", data.metadata.labels);
    update_session_info("blocklist", data.metadata.blocklist);
    update_session_info("collapse_threshold", data.metadata.collapse_threshold);

    update_session_info("curation_url", data.metadata.network);
    update_session_info("curation_version", data.metadata.curation_version);
    update_session_info("curation_date", data.metadata.curation_date);

    update_session_info("database_version", data.metadata.database_version);

    update_session_info("neighbors_url", data.metadata.neighbors_url);
    update_session_info("neighbors_version", data.metadata.neighbors_version);
    update_session_info("neighbors_date", data.metadata.neighbors_date);

    update_session_info("template_url", data.metadata.template_url);
    update_session_info("template_version", data.metadata.template_version);
    update_session_info("template_date", data.metadata.template_date);

    update_session_info("model_version", data.metadata.model_version);
    update_session_info("model_date", data.metadata.model_date);

    if (
      (transcriptomics === true) |
      (proteomics === true) |
      (metabolomics === true) |
      (get_session_info("provided_data") === true) |
      (get_session_info("provided_data") === "true") |
      (data.categories.length > 0)
    ) {
      $("#content").html(
        '<a href="motif.html"><div id="continue"><font size="3">Pattern Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="perturbations.html"><div id="continue"><font size="3">Perturbation Networks</font></div></a>'
      );
      update_session_info("provided_data", true);
    } else {
      $("#content").html(
        '<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="perturbations.html"><div id="continue"><font size="3">Perturbation Networks</font></div></a>'
      );
    }
  } catch (e) {
    //alert('Failed to open database: \n' + database_url)
    //var elem = document.getElementById("progressBar");
    //elem.style.width = "0%";
    //elem.innerHTML = "0%";
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
