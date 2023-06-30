/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) 2022 Metaboverse

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
var spawn = require("child_process").spawn;
var path = require("path");


// All of the Node.js APIs are available in the preload process.
// It has the same sandbox as a Chrome extension.
window.addEventListener("DOMContentLoaded", () => {
  const replaceText = (selector, text) => {
    const element = document.getElementById(selector);
    if (element) element.innerText = text;
  };

  for (const type of ["chrome", "node", "electron"]) {
    replaceText(`${type}-version`, process.versions[type]);
  }
});


async function showWarningPython() {
  const options = {
    title: 'Warning',
    message: 'Metaboverse cannot find a Python installation. Python is necessary for curating new data. \nhttps://www.python.org/downloads/',
  };
  const result = await ipcRenderer.invoke('show-warning-dialog', options);
  console.log(result); // This will print the index of the clicked button
}


async function showWarningDependencies() {
  const options = {
    title: 'Warning',
    message: "There was an error installing the required Python dependencies. Please check your internet connection and Python installation and try again.",  };
  const result = await ipcRenderer.invoke('show-warning-dialog', options);
  console.log(result); // This will print the index of the clicked button
}


function check_packages() {
  // Using the user's local Python, install the required packages from requirements.txt to a virtual environment 
  // Install packages from __dirname/../python/requirements.txt
  var process = spawn("python", ["-m", "pip", "install", "-r", path.join(__dirname, "..", "python", "requirements.txt")]);
  process.stdout.on("data", function (data) {
    console.log(data.toString());
    if (data.toString().toLowerCase().includes("error")) {
      console.log("Error installing packages");
      showWarningDependencies();
    } else {
      console.log("Packages installed");
    }
  })
}


function check_python() {
  var process = spawn("python", ["--version"]);
  process.stdout.on("data", function (data) {
    console.log(data.toString());
    //if (data.toString().toLowerCase().includes("python")) {
    //  console.log("Python is installed");
    //  check_packages();
    //} else{
      console.log("Python is not installed");

      showWarningPython();
    }  )
  }//);
//};

// Run check_python() and let the user continue if Python is installed
check_python()







