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
var $ = require("jquery");
var fs = require("fs");
var path = require("path");
var pixelWidth = require("string-pixel-width");
var { ipcRenderer } = require('electron');


function write_json(session_data) {
  ipcRenderer.invoke('get-paths').then((paths) => {
    fs.writeFileSync(paths.sessionFilePath, JSON.stringify(session_data), function(err) {
      if (err) throw err;
      console.log("Session data updated");
    });
  });
}

function update_session_info(key_update, value_update, abbrev_dict = null) {

  ipcRenderer.invoke('get-paths').then((paths) => {
    var session = JSON.parse(fs.readFileSync(paths.sessionFilePath).toString(), "utf8");
    session[key_update] = value_update;
    
    // Where database output location is, make this the output location
    if (key_update === "database_url") {
      file_path = value_update.substring(0, value_update.lastIndexOf(path.sep));
      session["output"] = file_path;
    } else if (key_update === "curation_url") {
      console.log(key_update + ":", value_update)
      file_path = value_update.substring(0, value_update.lastIndexOf(path.sep));
      session["output"] = file_path;
    } else {}

    if ((abbrev_dict != null) & (key_update === "organism")) {
      session["organism_id"] = abbrev_dict[value_update];
    }
    write_json(session);
    console.log("Updated session variable:", "\"" + key_update + "\"", "to ", "\"" + value_update + "\"")
  });
}

function get_session_info(key_update, callback) {
  ipcRenderer.invoke('get-paths').then((paths) => {
    var session = JSON.parse(fs.readFileSync(paths.sessionFilePath).toString());
    try {
      value = session[key_update];
      callback(null, value);
    } catch (e) {
      console.log("Could not update session variable: ", key_update)
      callback(e);
    }
  });
}

async function get_session_info_async(key_update) {
  let paths = await ipcRenderer.invoke('get-paths');
  try {
    let session = JSON.parse(fs.readFileSync(paths.sessionFilePath).toString());
    return session[key_update];
  } catch (e) {
    console.log("Could not update session variable: ", key_update);
    throw e; 
  }
}

async function get_default_async(key) {
  let paths = await ipcRenderer.invoke('get-paths');
  try {
    let session = JSON.parse(fs.readFileSync(paths.sessionFilePath).toString());
    let value = session[key];
    return value;
  } catch (e) {
    console.log("Could not parse session variable: ", key);
    throw e;  
  }
}

async function get_argument_async(key, callback) {
  let paths = await ipcRenderer.invoke('get-paths');
  try {
    let session = JSON.parse(fs.readFileSync(paths.sessionFilePath).toString());
    let value = session[key];
    if (value == null) {
      value = "None";
    }
    return value;
  } catch (e) {
    console.log("Could not parse session variable: ", key);
    throw e;  
  }
}


// http://www.alessioatzeni.com/blog/simple-tooltip-with-jquery-only-text/
$(document).ready(function() {
  // Tooltip only Text
  $(".info")
    .hover(
      function() {
        // Hover over code
        var title = $(this).attr("title");
        $(this)
          .data("tipText", title)
          .removeAttr("title");
        $('<p class="tooltip"></p>')
          .text(title)
          .appendTo("body")
          .fadeIn("slow");
      },
      function() {
        // Hover out code
        $(this).attr("title", $(this).data("tipText"));
        $(".tooltip").remove();
      }
    )
    .mousemove(function(e) {
      var mousex = e.pageX + 20; //Get X coordinates
      var mousey = e.pageY - 25; //Get Y coordinates
      $(".tooltip").css({
        top: mousey,
        left: mousex
      });
    });
});


function determineWidth(input) {
  if (pixelWidth(input) < 1662) {
    var mod_selection = input + "<br><br>";
  } else {
    var mod_selection = input;
  }

  return mod_selection;
}

// Populate dictionary to access component reactions for each pathway
function make_pathway_dictionary(data, database_key) {
  // Get pathway name and ID
  var pathways = data[database_key];
  var pathway_dict = {};
  for (var key in pathways) {
    pathway_dict[pathways[key]["name"]] = {
      id: pathways[key]["name"],
      reactome: pathways[key]["reactome"],
      reactions: pathways[key]["reactions"]
    };
  }

  return pathway_dict;
}

// Source: https://stackoverflow.com/a/48729396/9571488
function dateComponentPad(value) {
  var format = String(value);

  return format.length < 2 ? '0' + format : format;
}

// Source: https://stackoverflow.com/a/48729396/9571488
function formatDate(date) {
  var datePart = [ date.getFullYear(), date.getMonth() + 1, date.getDate() ].map(dateComponentPad);
  var timePart = [ date.getHours(), date.getMinutes(), date.getSeconds() ].map(dateComponentPad);

  return datePart.join('-') + '-' + timePart.join('-');
}


function get_script_name() {

  if (navigator.appVersion.indexOf("Win") != -1) {
    scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-windows.exe");
  } else if (navigator.appVersion.indexOf("Mac") != -1) {
    scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-darwin");
  } else if (navigator.appVersion.indexOf("Linux") != -1) {
    scriptFilename = path.join(__dirname, "..", "python", "metaboverse-cli-linux");
  } else {
    console.log("Unable to locate metaboverse-cli binary")
  }

  return scriptFilename;
}
