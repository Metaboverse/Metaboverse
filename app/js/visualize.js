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

var global_motifs = gatherMotifs(data, data.categories);
timecourse = checkCategories(data.categories, data.labels); //, data.names);

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
d3.select("#kNN_button").on("change", change);
d3.select("#hub_button").on("change", change);
d3.select("#stat_button").on("change", change);
