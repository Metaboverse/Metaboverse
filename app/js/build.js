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
var exec = require("child_process").exec;
var { ipcRenderer } = require("electron");
var path = require("path");
var fs = require("fs");
var $ = require("jquery");
var timer = 5000;

let database_url; 
let processed_bool; 
let provided_data_bool; 
let experiment_name;
let output;
let labels;
let blocklist;


console.log("Operating System information:")
console.log(navigator.appVersion)

var scriptFilename = get_script_name();


async function set_watch_files() {

  let paths = await ipcRenderer.invoke('get-paths');
  fs.copyFile(
    paths.progressFileTemplatePath,
    paths.progressFilePath,
    (err) => {
      if (err) throw err;
      console.log("Progress log file created");
    }
  );
  fs.watch(paths.progressFilePath, function(event, filename) {
    var elem = document.getElementById("progressBar");
    var sum_values = 0;
    var session = JSON.parse(fs.readFileSync(paths.progressFilePath).toString());
    for (j in session) {
      sum_values += session[j];
    }
    if (sum_values >= 100) {
      elem.style.width = "107%";
      elem.innerHTML = "Curation complete. Loading file...";
    } else {
      elem.style.width = sum_values + 7 + "%";
      elem.innerHTML = sum_values + "%";
    }
    
    if (sum_values >= 100) {
      setTimeout(function(){
        displayOptions();
      }, timer);
    }
  });
}
set_watch_files();

async function get_omics_bool_async() {
  var transcriptomics = false;
  var proteomics = false;
  var metabolomics = false;
  var transcriptomics_bool = await get_default_async("transcriptomics");
  var proteomics_bool = await get_default_async("proteomics");
  var metabolomics_bool = await get_default_async("metabolomics");

  if (transcriptomics_bool != null) {
    transcriptomics = true;
  }
  if (proteomics_bool != null) {
    proteomics = true;
  }
  if (metabolomics_bool != null) {
    metabolomics = true;
  }

  return {transcriptomics, proteomics, metabolomics};
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

async function write_log(command, stdout, stderr) {

  experiment_name = await get_argument_async("experiment_name");
  experiment_name = experiment_name.replace(/\s+/g, '-');
  let today = new Date().toISOString().slice(0, 10);

  output = await get_argument_async("output");
  fs.writeFileSync(
    path.join(output.replace(/\"/g, ''), "metaboverse_session_" + experiment_name + ".log"),
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

runBuild = async function(_callback) {

  console.log("Parsing build arguments...")
  var elem = document.getElementById("progressBar");
  elem.style.width = "7%";
  elem.innerHTML = "0%";

  let paths = await ipcRenderer.invoke('get-paths');
  processed_bool = await get_session_info_async("processed");
  if (processed_bool == true) {
    elem.style.width = "100%";
    elem.innerHTML = "100%";
    setTimeout(function(){
      displayOptions();
    }, timer);
  } else {
    curated = await get_argument_async("curation_url");
    labels = await get_argument_async("labels");
    blocklist = await get_argument_async("blocklist");
    database_url = await get_argument_async("database_url");

    if (labels == "" | labels == undefined) {
      labels = "0"
    } else {
      labels = labels.replace(/\s/g, '');
    }

    if (blocklist == "" | blocklist == undefined) {
      blocklist = "no_blocklist"
    } else {
      blocklist = blocklist.replace(/\s/g, '');
    }

    if (String(curated) != "None") {
      let db_url;
      if (String(database_url) == "None" | String(database_url) == "null") {
        db_url = "\"find\"";
      } else {
        db_url = String(decodeURIComponent(database_url));
      }
      graphDictionary = {
        output: await get_argument_async("output"),
        organism_id: "\"find\"",
        output_file: db_url,
        organism_curation_file: curated,
        neighbor_dictionary_file: await get_argument_async("neighbors_url"),
        graph_template_file: await get_argument_async("template_url"),
        transcriptomics: await get_argument_async("transcriptomics"),
        proteomics: await get_argument_async("proteomics"),
        metabolomics: await get_argument_async("metabolomics"),
        database_source: await get_argument_async("database_source"),
        //additional_reactions: await get_argument("additional_reactions"),
        experiment_type: await get_argument_async("experiment_type"),
        experiment_name: await get_argument_async("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        force_new_curation: await get_argument_async("forceNewCuration"),
        collapse_with_modifiers: await get_argument_async("collapseWithModifiers"),
        broadcast_genes: await get_argument_async("broadcastGeneExpression"),
        broadcast_metabolites: await get_argument_async("broadcastMetabolites"),
        collapse_threshold: parseFloat(await get_argument_async("collapse_threshold")) / 100,
        progress_log: "\"" + paths.progressFilePath + "\"",
        session_data: "\"" + paths.sessionFilePath + "\""
      }
    } else {
      graphDictionary = {
        output: await get_argument_async("output"),
        organism_id: await get_argument_async("organism_id"),
        output_file: await get_argument_async("database_url"),
        organism_curation_file: await get_argument_async("curation_url"),
        neighbor_dictionary_file: await get_argument_async("neighbors_url"),
        graph_template_file: await get_argument_async("template_url"),
        transcriptomics: await get_argument_async("transcriptomics"),
        proteomics: await get_argument_async("proteomics"),
        metabolomics: await get_argument_async("metabolomics"),
        database_source: await get_argument_async("database_source"),
        //additional_reactions: await get_argument_async("additional_reactions"),
        experiment_type: await get_argument_async("experiment_type"),
        experiment_name: await get_argument_async("experiment_name"),
        labels: labels,
        blocklist: blocklist,
        force_new_curation: await get_argument_async("forceNewCuration"),
        collapse_with_modifiers: await get_argument_async("collapseWithModifiers"),
        broadcast_genes: await get_argument_async("broadcastGeneExpression"),
        broadcast_metabolites: await get_argument_async("broadcastMetabolites"),
        collapse_threshold: parseFloat(await get_argument_async("collapse_threshold")) / 100,
        progress_log: "\"" + paths.progressFilePath + "\"",
        session_data: "\"" + paths.sessionFilePath + "\""
      }
    }
    var cmd = parseCommand(graphDictionary);
    console.log("Running command:\n" + scriptFilename + " curate" + cmd);
    execute(scriptFilename + " curate" + cmd, output => {
      update_session_info("processed", true);
    });
    return _callback;
  }
};

async function displayOptions() {

  update_session_info("current_pathway", null);
  database_url = await get_session_info_async("database_url");
  try {
    var data = JSON.parse(fs.readFileSync(database_url).toString());
    update_session_info("curation_url", data.metadata.organism_curation_file);
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

    provided_data_bool = await get_session_info_async("provided_data");
    var omics = await get_omics_bool_async();
    var transcriptomics = omics[0];
    var proteomics = omics[1];
    var metabolomics = omics[2];

    if (
      (transcriptomics === true) |
      (proteomics === true) |
      (metabolomics === true) |
      (provided_data_bool == true) |
      (provided_data_bool == "true") |
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
    var elem = document.getElementById("progressBar");
    elem.style.width = "0%";
    elem.innerHTML = "";
  } catch (e) {
    //alert('Failed to open database: \n' + database_url)
    //var elem = document.getElementById("progressBar");
    //elem.style.width = "0%";
    //elem.innerHTML = "0%";
  }
}

window.addEventListener("load", function(event) {
  $('#content').append('')
  console.log("Loading build script...")
  runBuild()

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session()
  }

})
