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

// Change user selection based on input
function selectPathway() {
  var selection = document.getElementById("pathwayMenu").value;
  console.log("Selection: " + selection);
  return selection;
}

function selectSuperPathway(selector) {
  var superSelection = document.getElementById("superPathwayMenu").value;
  console.log("Super Pathway: " + superSelection);
  return superSelection;
}

// Populate dictionary to access component reactions for each super-pathway
function make_superPathway_dictionary(data) {
  // Get pathway name and ID
  var superPathways = data.super_pathways;
  var superPathwayDict = {};
  for (var key in superPathways) {
    superPathwayDict[superPathways[key]["name"]] = {
      id: superPathways[key]["name"],
      pathway_id: superPathways[key]["id"],
      reactome_id: superPathways[key]["reactome"],
      reactions: superPathways[key]["reactions"]
    };
  }

  return superPathwayDict;
}

// Make Pathway menu for users to
function make_menu(
  pathway_dict,
  selector,
  opening_message,
  provide_all = false
) {
  // Get species names (keys) as list
  pathways_list = [];
  pathways_list = Object.getOwnPropertyNames(pathway_dict).map(function(k) {
    return k;
  });

  // Case-insensitive sort (Source: https://stackoverflow.com/a/26145462/9571488)
  pathways_list.sort( function(s1, s2) {
    var l=s1.toLowerCase(), m=s2.toLowerCase();
    return l===m?0:l>m?1:-1;
  });


  if (provide_all === true) {
    pathways_list.unshift("All pathways");
    pathways_list.unshift("All entities");
  }
  pathways_list.unshift(opening_message); // Add select prompt to menu bar

  // Generate drop-down menu for species select
  menu = [];
  menu = document.getElementById(selector);
  for (var i = 0; i < pathways_list.length; i++) {
    var option = document.createElement("option");
    option.innerHTML = pathways_list[i];
    option.value = pathways_list[i];
    menu.appendChild(option);
  }
}

// Remove elements from pathway menu
function emptyMenu(selectMenu) {
  var i;
  for (i = selectMenu.options.length - 1; i >= 0; i--) {
    selectMenu.remove(i);
  }
}

function parsePathways(pathway_dict, selectedReactions) {
  var parsed_pathway_dict = {};
  // adapted from https://stackoverflow.com/a/53606357
  let checker = (arr, target) => target.every(v => arr.includes(v));

  for (pathway in pathway_dict) {
    belongs = checker(selectedReactions, pathway_dict[pathway]["reactions"]);

    if (belongs === true) {
      parsed_pathway_dict[pathway] = pathway_dict[pathway];
    }
  }

  return parsed_pathway_dict;
}
