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

// Modules to control application life and create native browser window
const {
  app,
  BrowserWindow,
  ipcMain,
  dialog, 
} = require("electron");
const path = require("path");
const fs = require("fs");

// Get app and user paths
const appPath = app.getAppPath();
const userDataPath = app.getPath('userData');

const sessionFileTemplatePath = path.join(appPath, "data", "session_data_template.json");
const sessionFilePath = path.join(userDataPath, "session_data.json");

const progressFileTemplatePath = path.join(appPath, "data", "progress_log_template.json");
const progressFilePath = path.join(userDataPath, "progress_log.json");

ipcMain.handle('get-paths', () => {
  return {
    sessionFileTemplatePath,
    sessionFilePath,
    progressFileTemplatePath,
    progressFilePath
  };
});


// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow;

function createWindow() {
  // Create the browser window.
  mainWindow = new BrowserWindow({
    width: 1450,
    height: 1000,
    minWidth: 900,
    minHeight: 600,
    webPreferences: {
      preload: path.join(__dirname, "js", "preload.js"),
      enableRemoteModule: true,
      worldSafeExecuteJavaScript: true,
      contextIsolation: false,
      nodeIntegration: true
    },
    dependencies: {
      zerorpc: path.join("fyears", "zerorpc-node")
    },
    icon: path.join(__dirname, "data", "icon", "png", "icon_64x64.png")
  });


  // Show devtools
  //mainWindow.webContents.openDevTools();

  // and load the index.html of the app.
  mainWindow.loadFile(path.join(__dirname, "html", "index.html"));
  
  // Emitted when the window is closed.
  mainWindow.on("closed", function() {
    fs.unlinkSync(sessionFilePath);
    mainWindow = null;
  });

  mainWindow.webContents.setFrameRate(60)
}


// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on("ready", createWindow);

// Quit when all windows are closed.
app.on("window-all-closed", function() {
  // On macOS it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  app.quit();
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
  pyProc = require("child_process").spawn("python", [script, port], {shell: true});
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


// Copy session info template each time the app is launched
fs.copyFile(
  sessionFileTemplatePath,
  sessionFilePath,
  err => {
    if (err) throw err;
    console.log("Session data file was copied for this session");
  }
);


ipcMain.handle('open-file-dialog-mvrs', async (event) => {
  const result = await dialog.showOpenDialog({
    properties: ['openFile'],
    filters: [
      { name: 'Documents', extensions: ['mvrs', 'json'] },
      // You can add more types if you want
    ],
  });

  if (result.canceled) {
    return;
  } else {
    return result.filePaths;
  }
});

ipcMain.handle('open-file-dialog-mvdb', async (event) => {
  const result = await dialog.showOpenDialog({
    properties: ['openFile'],
    filters: [
      { name: 'Documents', extensions: ['mvdb', 'json'] },
      // You can add more types if you want
    ],
  });

  if (result.canceled) {
    return;
  } else {
    return result.filePaths;
  }
});

ipcMain.handle('open-file-dialog-tsv', async (event) => {
  const result = await dialog.showOpenDialog({
    properties: ['openFile'],
    filters: [
      { name: 'Documents', extensions: ['txt', 'tsv'] },
      // You can add more types if you want
    ],
  });

  if (result.canceled) {
    return;
  } else {
    return result.filePaths;
  }
});

ipcMain.handle('save-file-dialog-mvrs', async (event) => {
  const result = await dialog.showSaveDialog({
    defaultPath: 'output.mvrs'
  });

  if (result.canceled) {
    return;
  } else {
    return result.filePath;
  }
});

ipcMain.handle('save-file-dialog-svg', async (event) => {
  const result = await dialog.showSaveDialog({
    defaultPath: 'image.svg',
    filters: [
      { name: 'Documents', extensions: ['svg'] },
      // You can add more types if you want
    ],
  });

  if (result.canceled) {
    return;
  } else {
    return result.filePath;
  }
});

ipcMain.handle('save-file-dialog-table', async (event, output_data) => {
  const result = await dialog.showSaveDialog({
    defaultPath: '..' + path.sep + '..' + path.sep,
    properties: ['createDirectory'],
    filters: [
      { name: 'tab-delimited (*.tsv; *.txt)', extensions: ['txt', 'tsv'] },
    ],
  });

  if (result.canceled) {
    return;
  } else {
    // write to file
    let hasExtension = /\.[^\/\\]+$/.test(result.filePath);
    if (hasExtension === false) {
      result.filePath = `${result.filePath}.${"txt"}`;
    }

    fs.writeFileSync(result.filePath, output_data, function(err) {
      if (err) throw err;
      console.log("Table exported");
    });
    return result.filePath;
  }
});

ipcMain.handle('show-warning-dialog', async (event, options) => {
  const result = await dialog.showMessageBox({
    type: 'warning',
    title: options.title,
    message: options.message,
    buttons: ['OK'],
    defaultId: 0, // The index of the button to be selected by default
    cancelId: 0, // The index of the button to be triggered when the dialog is canceled
    noLink: true // This prevents Electron from automatically adding a link to the dialog's message
  });

  return result.response; // This will be the index of the clicked button
});

ipcMain.on('reload-window', (event, arg) => {
  mainWindow.webContents.reload();
});