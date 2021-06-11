/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2021 Jordan A. Berg
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

// Hyperlinks listener
const newBrowserSettings = 'top=500,left=200,frame=false,nodeIntegration=no,enableRemoteModule=no,worldSafeExecuteJavaScript=yes,contextIsolation=yes';
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
        'https://metaboverse.readthedocs.io/en/latest/content/general-usage.html',
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
        'https://metaboverse.readthedocs.io/en/latest/content/general-usage.html#regulatory-hotspot-identification-pattern-analysis',
        '_blank',
        newBrowserSettings
      )
    } else if (page === "visualize.html") {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/general-usage.html#general-pathway-exploration',
        '_blank',
        newBrowserSettings
      )
    } else if (page === "perturbations.html") {
      window.open(
        'https://metaboverse.readthedocs.io/en/latest/content/general-usage.html#perturbation-network-modeling',
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