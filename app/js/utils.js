update_session_info/*
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
var $ = require("jquery");
var fs = require("fs");
var pixelWidth = require("string-pixel-width");

var app = require("electron").remote.app;
var userDataPath = app.getPath("userData");
var session_file = userDataPath + "/session_data.json";

console.log(session_file);

function write_json(session_data) {
  fs.writeFileSync(session_file, JSON.stringify(session_data), function(err) {
    if (err) throw err;
    console.log("Session data updated");
  });
}

function update_session_info(key_update, value_update, abbrev_dict = null) {

  try {
    var session = JSON.parse(fs.readFileSync(session_file).toString(), "utf8");
    console.log(session)
    console.log(key_update)
    console.log(value_update)
    session[key_update] = value_update;

    // Where database output location is, make this the output location
    if (key_update === "database_url") {
      path = value_update.substring(0, value_update.lastIndexOf("/"));
      session["output"] = path + "/";
    } else if (key_update === "curation_url") {
      path = value_update.substring(0, value_update.lastIndexOf("/"));
      session["output"] = path + "/";
    } else {}

    if ((abbrev_dict !== null) & (key_update === "organism")) {
      session["organism_id"] = abbrev_dict[value_update];
    }
    write_json(session);
  } catch (e) {
    console.log("Could not update session variable: ", key_update)
  }

}

function get_session_info(key_update) {

  try {
    var session = JSON.parse(fs.readFileSync(session_file).toString());
    value = session[key_update];
    return value;
  } catch (e) {
    console.log("Could not update session variable: ", key_update)
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
      $(".tooltip").css({ top: mousey, left: mousex });
    });
});

function getDefault(key) {
  var session = JSON.parse(fs.readFileSync(session_file).toString());
  value = session[key];

  return value;
}

function getArgument(key) {
  var session = JSON.parse(fs.readFileSync(session_file).toString());
  value = session[key];

  if (value === null) {
    value = "None";
  }

  return value;
}

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
