/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) Jordan A. Berg, The University of Utah

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
var path = require("path");
var semver = require('semver');
var $ = require("jquery");
var reactome_api = "https://reactome.org/ContentService/data/species/all";

let defaultValue; 


window.addEventListener("load", function(event) {
  event.preventDefault();
  event.stopPropagation();

  get_session_info("launched", (err, value) => {
    if (err) {
      console.log(err)
    } else {
      if (value === false) {
        update_session_info("launched", true);
        $.ajax({
          url: path.join("..", "__version__.txt"),
          success: function(version) {
            $.getJSON("https://api.github.com/repos/Metaboverse/Metaboverse/tags", function(d) {
              let _v = String(version.trim().replace(/[^0-9.]/g,''))
    
              let avail_versions = [];
              let version_dict = {};
              for (_k in d) {
                let _this_version = String(d[_k].name.trim().replace(/[^0-9.]/g,''))
                avail_versions.push(_this_version);
                version_dict[_this_version] = d[_k].name;
              }
    
              // Source: https://stackoverflow.com/a/40201629
              avail_versions = avail_versions.map( a => a.split('.').map( n => +n+100000 ).join('.') ).sort()
                      .map( a => a.split('.').map( n => +n-100000 ).join('.') );
              let _c = avail_versions[avail_versions.length - 1];
    
              if (semver.gt(_c, _v)) {
                alert("A more current version of Metaboverse is available:\n\n" + version_dict[_c] + "\n\n\nPlease download this version then close this window and launch the new version.")
                window.open(
                  "https://github.com/Metaboverse/Metaboverse/releases/tag/" + version_dict[_c],
                  "_blank",
                  "top=500,left=200,frame=true,nodeIntegration=no,enableRemoteModule=no,worldSafeExecuteJavaScript=yes,contextIsolation=yes");
              }
            })
          }
        });
      }
    }
  });
  $('#content').append('<a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>');

  get_session_info("database_url", (err, value) => {
    if (err) {
      console.log(err)
    } else {
      defaultValue = "No file selected";
      if (value != "" && value != null) {
        defaultValue = value;
        $("#content").html(
          '<a href="motif.html"><div id="continue"><font size="3">Pattern Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Explore</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="perturbations.html"><div id="continue"><font size="3">Perturbation Networks</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>'
        );
      }
      $('#selectedFile').append('<font size="2">' + defaultValue + '</font>');
    }
  });
  

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session();
  }

  document.getElementById("dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    document.getElementById('database-input').click();
  }

  document.getElementById("database-input").onchange = async function(event) {
    event.preventDefault();
    event.stopPropagation();
  
    let inputVal = document.getElementById("database-input").value.split(".");
  
    if (isInvalidFileType(inputVal)) {
      alert("Input is not a .mvrs file. You must upload the correct file type for the analyses to work.");
      return;
    }
    
    try {
      let file = event.srcElement.files[0];
      await processFile(file);
    } catch (error) {
      console.error(error);
      alert("Unable to load provided file into session info.");
      window.location.reload(false);
    }
  };
  
  function isInvalidFileType(inputVal) {
    return inputVal[inputVal.length - 1] !== "mvrs";
  }
  
  async function processFile(file) {
    console.log("The file you dragged: ", file.path);
    await update_session_info("database_url", file.path);
    $('#selectedFile').html('<font size="2">' + file.path + '</font>');
    
    let data = JSON.parse(fs.readFileSync(file.path).toString());
  
    $("#content").html(generateContentHTML(data.categories.length > 0));
  
    let abbreviation_dict = {};
    await $.getJSON(reactome_api, function(d) {
      d.forEach((datum) => abbreviation_dict[datum["abbreviation"]] = datum["displayName"]);
      update_session_info("organism", abbreviation_dict[data.metadata.organism_id]);
    });
    console.log(data.metadata);
  
    let metadataKeys = Object.keys(data.metadata);
    for (let key of metadataKeys) {
      await update_session_info(key, data.metadata[key]);
      console.log(key, data.metadata[key]);
    }
    // Read session_data and print 
    print_session_info();
  }
  
  function generateContentHTML(hasCategories) {
    if (hasCategories) {
      return '<a href="motif.html"><div id="continue"><font size="3">Pattern Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="visualize.html"><div id="continue"><font size="3">Explore</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="perturbations.html"><div id="continue"><font size="3">Perturbation Networks</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>';
    } else {
      return '<a href="visualize.html"><div id="continue"><font size="3">Explore</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="perturbations.html"><div id="continue"><font size="3">Perturbation Networks</font></div></a></br></br><a href="curate.html"><div id="continue"><font size="3">Skip</font></div></a>';
    }
  }
});

ipcRenderer.on("file-input", (event, data) => {
  $("#file-input").text(data);
});
