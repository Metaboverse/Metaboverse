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

var $ = require("jquery");
var reactome_api = "https://reactome.org/ContentService/data/species/all";

// Drop pre-existing metabolic network database for further analysis
window.addEventListener("load", function(event) {
  document.getElementById("database-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("database-input").value.split(".");
    console.log(inputVal[inputVal.length - 1])
    if (inputVal[inputVal.length - 1] !== "json") {
      alert(
        "Input is not a .json file. You must upload the correct file type for the analyses to work."
      );

    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f.path);
        update_session_info("database_url", f.path);

        $('#selectedFile').html('<font size="2">' + f.path + '</font>');
        var data = JSON.parse(fs.readFileSync(f.path).toString());
        if (data.categories.length > 0) {
          $("#content").html(
            '<a href="motif.html"><div id="continue"><font size="3">View Motif Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>'
          );
        } else {
          $("#content").html(
            '<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>'
          );
        }
        // Load session data values into current session
        var abbreviation_dict = {};
        $.getJSON(reactome_api, function(d) {
          d.forEach(function(datum) {
            abbreviation_dict[datum["abbreviation"]] = datum["displayName"];
          });
          update_session_info("organism",
            abbreviation_dict[data.metadata.species_id]);
        });
        update_session_info("curation_url", data.metadata.organism_curation);
        update_session_info("output", data.metadata.output);
        update_session_info("organism_id", data.metadata.species_id);
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
            data.metadata.broadcastGeneExpression);
        update_session_info("labels", data.metadata.labels);
        update_session_info("blacklist", data.metadata.blacklist);
        update_session_info("database_date", data.metadata.database_date);
        update_session_info("curation_date", data.metadata.curation_date);
        update_session_info("reactome_version", data.metadata.reactome_version);

      } catch (error) {
        console.log(error);
        alert(
          "Unable to load provided file into session info."
        );
        window.location.reload(false);
      }
    }


  };
});

ipcRenderer.on("file-input", (event, data) => {
  $("#file-input").text(data);
});
