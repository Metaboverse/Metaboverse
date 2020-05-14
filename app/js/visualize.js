/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
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

var selection = null;
var superSelection = null;
var selector = "#graph";
var _width = window.innerWidth;
var _height = window.innerHeight - 75;

// find global motifs at beginning, save as list, check each graph for members
function gatherMotifs(data) {

  let expression_dict = {};
  for (let x in data.nodes) {
    let id = data.nodes[x]['id'];
    let expression = data.nodes[x]['values'];
    expression_dict[id] = expression;
  }

  let threshold = 1;
  let motifs_Avg = motifSearch_Avg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary)

  let motifs_MaxMax = motifSearch_MaxMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary)

  let motifs_MaxMin = motifSearch_MaxMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary)

  let motifs_Sustained = motifSearch_Sustained(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary)

  let all_motifs = motifs_Avg.concat(
    motifs_MaxMax,
    motifs_MaxMin,
    motifs_Sustained);

  let global_motifs = [];
  for (let m in all_motifs) {
    global_motifs.push(all_motifs[m].id);
  }
  return global_motifs;
}

// MAIN
database_url = get_session_info("database_url");
console.log("Database path: " + database_url);

var data = JSON.parse(fs.readFileSync(database_url).toString());
var pathway_dict = make_pathway_dictionary(
  data,
  'pathway_dictionary');
var collapsed_pathway_dict = make_pathway_dictionary(
  data,
  'collapsed_pathway_dictionary');
var superPathwayDict = make_superPathway_dictionary(data);

var global_motifs = gatherMotifs(data);

var timecourse = checkCategories(data.categories); //, data.names);

make_menu(
  superPathwayDict,
  "superPathwayMenu",
  "Select a category...",
  (provide_all = true)
);

console.log(data)

let current_pathway = get_session_info("current_pathway");
if ((current_pathway !== null) & (current_pathway !== "null")) {
  change();
} else {}

d3.select("#superPathwayMenu").on("change", changeSuper);
d3.select("#pathwayMenu").on("change", change);
