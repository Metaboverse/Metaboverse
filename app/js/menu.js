/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019 Jordan A. Berg
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
  pathways_list.sort();
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
