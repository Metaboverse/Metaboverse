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

        $('#selectedFile').replaceWith('<font size="2">' + f.path + '</font>');
        var data = JSON.parse(fs.readFileSync(f.path).toString());
        if (data.categories.length > 0) {
          $("#content").replaceWith(
            '<a href="motif.html"><div id="continue"><font size="3">View Motif Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>'
          );
        } else {
          $("#content").replaceWith(
            '<a href="visualize.html"><div id="continue"><font size="3">Visualize</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="connections.html"><div id="continue"><font size="3">Connectivity</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>'
          );
        }


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
