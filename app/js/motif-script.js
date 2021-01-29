/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2021  Youjia Zhou, Jordan A. Berg
  zhou325 <at> sci <dot> utah <dot> edu
  jordan <dot> berg <at> biochem <dot> utah <dot> edu

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
var fs = require('fs');
var path = require("path");
var app = require('electron').remote.app;

var userDataPath = app.getPath('userData');
var session_file = userDataPath + path.sep + "session_data.json";
let database_url = JSON.parse(
  fs.readFileSync(session_file).toString())["database_url"];

showMotifs = function(_callback) {

  d3.json(database_url).then(data => {
    let metaGraph = new MetaGraph(data);
  });
  return _callback;
}

window.addEventListener("load", function(event) {
  showMotifs()
})
