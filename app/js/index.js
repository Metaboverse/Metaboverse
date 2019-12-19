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
const {app, BrowserWindow} = require('electron')
const {ipcRenderer, remote} = require('electron')
var $ = require('jquery')

// Drop pre-existing metabolic network database for further analysis
window.addEventListener('load', function(event) {

  document.getElementById("dropDatabase").onchange = function (event) {

  event.preventDefault()
  event.stopPropagation()

  var inputVal = document.getElementById("dropDatabase").value.split(".")

  if (inputVal[inputVal.length-1] !== "json") {
    alert('Input is not a .json file. You must upload the correct file type for the analyses to work. Restarting page...')
    window.location.reload(false);
  } else {

    try {
      f = event.srcElement.files[0]
      console.log('The file you dragged: ', f)
      path = f.path;

      update_session_info("database_url", path)

      $('#content').replaceWith('<a href="../html/motif.html"><div id="continue"><font size="3">Run Motif Analysis</font></div></a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="../html/visualize.html"><div id="continue"><font size="3">Visualize</font></div></a></br></br><a href="../html/curate.html"><div id="continue"><font size="3">Skip</font></div></a>')

    } catch (error) {
      console.log(error)
      alert('Input is not a .json file. You must upload the correct file type for the analyses to work.')
    }
  }
}})

ipcRenderer.on('dropDatabase', (event, data) => {
  $("#dropDatabase").text(data)
})
