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
var path = require('path')
var fs = require('fs')

var app = require('electron').remote.app
var basePath = app.getAppPath()

var userDataPath = app.getPath('userData');
var session_file = userDataPath + "/session_data.json"

// Copy session info template each time the app is launched

fs.copyFile(basePath + "/data/session_data_template.json" , session_file, (err) => {

  if (err) throw err;
  console.log('Session data file was copied for this session');
});