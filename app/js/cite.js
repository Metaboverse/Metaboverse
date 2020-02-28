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
var $ = require("jquery");

// Get Metaboverse version
$.ajax({
  url: "../../__version__.txt",
  success: function(version) {
    document.getElementById("getVersion").innerHTML = version;
  }
});

// Get current Reactome version
$.ajax({
  url: "https://reactome.org/tag/release",
  success: function(result) {
    var versions = [];
    let array = [...result.matchAll("Version (.*) Released")];

    for (a in array) {
      versions.push(array[a][1]);
    }

    document.getElementById("getReactomeVersion").innerHTML = Math.max(
      ...versions
    );
  }
});
