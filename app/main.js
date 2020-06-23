/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2020 Jordan A. Berg
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

// Modules to control application life and create native browser window
const { app, BrowserWindow } = require("electron");
const { ipcRenderer } = require("electron");
const path = require("path");
const fs = require("fs");
const ipcMain = require("electron").ipcMain;
const { dialog } = require("electron");

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow;

function createWindow() {
  // Create the browser window.
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 1400,
    minHeight: 900,
    webPreferences: {
      preload: path.join(__dirname, "js/preload.js"),
      nodeIntegration: true
    },
    dependencies: {
      zerorpc: "fyears/zerorpc-node"
    },
    icon: path.join(__dirname, "data/icon/metaboverse_logo.icns")
  });

  // and load the index.html of the app.
  mainWindow.loadFile("html/index.html");

  // Open the DevTools.
  //mainWindow.webContents.openDevTools()
  //mainWindow.webContents.closeDevTools()

  // Emitted when the window is closed.
  mainWindow.on("closed", function() {

    /*
    var userDataPath = app.getPath("userData");
    var session_file = userDataPath + "/session_data.json";
    var session = JSON.parse(fs.readFileSync(session_file).toString());
    var output = session["output"];

    if (output !== null) {
      fs.writeFile(
        output + "session_data.json",
        JSON.stringify(session),
        function(err) {
          if (err) throw err;
          console.log("Session data written to user output directory");
        }
      );
    } else {
      const options = {
        buttons: ["Okay"],
        title: "Alert",
        message: "Valid output directory cannot be found.",
        detail:
          "Closing will lose all session data.\n\nCopy the following information to a .json file to save your session information:\n" +
          JSON.stringify(session)
      };

      dialog.showMessageBox(null, options, (response, checkboxChecked) => {
        console.log(response);
        console.log(checkboxChecked);
      });
    }
    */

    fs.unlinkSync(session_file);
    mainWindow = null;
  });
}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on("ready", createWindow);

// Quit when all windows are closed.
app.on("window-all-closed", function() {
  // On macOS it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  if (process.platform !== "darwin") app.quit();
});

app.on("activate", function() {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) createWindow();
});

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.
let pyProc = null;
let pyPort = null;

const selectPort = () => {
  pyPort = 4242;
  return pyPort;
};

const createPyProc = () => {
  let port = "" + selectPort();
  let script = path.join(__dirname, "pycalc", "api.py");
  pyProc = require("child_process").spawn("python", [script, port]);
  if (pyProc != null) {
    console.log("child process success");
  }
};

const exitPyProc = () => {
  pyProc.kill();
  pyProc = null;
  pyPort = null;
};

app.on("ready", createPyProc);
app.on("will-quit", exitPyProc);

//var app = require("electron").remote.app;
var basePath = app.getAppPath();

var userDataPath = app.getPath("userData");
var session_file = userDataPath + "/session_data.json";

// Copy session info template each time the app is launched

fs.copyFile(
  basePath + "/data/session_data_template.json",
  session_file,
  err => {
    if (err) throw err;
    console.log("Session data file was copied for this session");
  }
);
