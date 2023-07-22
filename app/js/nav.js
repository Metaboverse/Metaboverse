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

// Hyperlinks listener
const newBrowserSettings = 'top=500,left=200,frame=true,nodeIntegration=no,enableRemoteModule=no,worldSafeExecuteJavaScript=yes,contextIsolation=yes';
window.addEventListener("load", function(event) {
  event.preventDefault();
  event.stopPropagation();

  var user_path = window.location.pathname;
  var page = user_path.split('/').pop();

  document.getElementById("issues_link").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    window.open(
      'https://github.com/Metaboverse/Metaboverse/issues',
      '_blank',
      newBrowserSettings
    )
  }

  document.getElementById("citation_link").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    window.open(
      'https://metaboverse.readthedocs.io/en/latest/content/cite.html',
      '_blank',
      newBrowserSettings
    )
  }

  document.getElementById("docs_link").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    window.open(
      'https://metaboverse.readthedocs.io/en/latest/',
      '_blank',
      newBrowserSettings
    )
  }

  if (page !== "motif.html" && page !== "visualize.html" && page !== "perturbations.html") {
    document.getElementById("usage_link").onclick = function(event) {
      event.preventDefault();
      event.stopPropagation();

      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/walkthrough.html',
        '_blank',
        newBrowserSettings
      )
    }
  }

  document.getElementById("overview_link").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    if (page === "motif.html") {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/walkthrough.html#reaction-pattern-analysis',
        '_blank',
        newBrowserSettings
      )
    } else if (page === "visualize.html") {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/walkthrough.html#general-pathway-exploration',
        '_blank',
        newBrowserSettings
      )
    } else if (page === "perturbations.html") {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/walkthrough.html#perturbation-network-modeling',
        '_blank',
        newBrowserSettings
      )
    } else {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/overview.html',
        '_blank',
        newBrowserSettings
      )
    }
  }

});