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
var value_threshold = 2.5; // ~1.5 fold
var stat_threshold = 0.01;
var position = 0;

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

var global_motifs = gatherMotifs(data);

//var timecourse = checkCategories(data.categories);

console.log(data)

// Allow absolute threshold or p-value
// Allow omics selection
graph_genes = true;
collapse_reactions = true;

//data.collapsed_reaction_dictionary;
//data.mod_collapsed_pathways;
collapsed_pathway_dict = make_pathway_dictionary(
  data,
  "collapsed_pathway_dictionary"
)

var path_mapper = data.motif_reaction_dictionary

node_dict = {};
for (let node in data.nodes) {
  node_dict[data.nodes[node]['id']] = data.nodes[node];
}

// Generate expression and stats dictionaries
let expression_dict = {};
let stats_dict = {};
for (let x in data.nodes) {
  let id = data.nodes[x]['id'];
  let expression = data.nodes[x]['values'];
  let stats = data.nodes[x]['stats'];

  if (!expression.includes(null)) {
    expression_dict[id] = expression;
  }
  if (!stats.includes(null)) {
    stats_dict[id] = stats;
  }

}

// For quicker searching with new motifs
let reaction_entity_dictionary = {};
for (rxn in data.collapsed_reaction_dictionary) {

  let r = data.collapsed_reaction_dictionary[rxn];
  let entities = [];
  for (x in r["reactants"]) {
    entities.push(r["reactants"][x]);
  }
  for (x in r["products"]) {
    entities.push(r["products"][x]);
  }
  for (x in r["modifiers"]) {
    entities.push(r["modifiers"][x][0]);
  }
  for (x in r["additional_components"]) {
    entities.push(r["additional_components"][x]);
  }

  reaction_entity_dictionary[r.id] = entities;
}

function collect_perturbations(
      reaction_entity_dictionary,
      mapping_dictionary,
      threshold,
      type) {

  let perturbed_reactions = [];
  for (rxn in reaction_entity_dictionary) {

    let perturbed_true = false;
    for (e in reaction_entity_dictionary[rxn]) {

      let entity = reaction_entity_dictionary[rxn][e];
      if (entity in mapping_dictionary) {

        if (type === "value") {
          if (Math.abs(mapping_dictionary[entity][position]) >= threshold) {
            perturbed_true = true;
          }
        } else if (type === "stat") {
          if (Math.abs(mapping_dictionary[entity][position]) <= threshold) {
            perturbed_true = true;
          }
        }
      }
    }

    if (perturbed_true === true) {
      perturbed_reactions.push(rxn);
    }
  }

  return perturbed_reactions;
}

function show_graph(data, perturbed_rxns) {

  console.log("Plotting", perturbed_rxns.length, "perturbed reactions...")
  // Parse through each reaction listed and get the component parts
  let components = [];
  var rxn = 0;
  for (rxn in perturbed_rxns) {

    var target_rxns = data.collapsed_reaction_dictionary[perturbed_rxns[rxn]];

    components.push(target_rxns.id);
    for (x in target_rxns["reactants"]) {
      components.push(target_rxns["reactants"][x]);
    }
    for (x in target_rxns["products"]) {
      components.push(target_rxns["products"][x]);
    }
    for (x in target_rxns["modifiers"]) {
      components.push(target_rxns["modifiers"][x][0]);
    }
    for (x in target_rxns["additional_components"]) {
      components.push(target_rxns["additional_components"][x]);
    }
  }

  var elements = get_nodes_links(data, components);
  var new_nodes = elements[0];
  var new_links = elements[1];

  // Initialize variables
  var node_dict = {};
  var type_dict = {};

  var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
  var node_dict = node_elements[0];
  var type_dict = node_elements[1];
  var display_analytes_dict = node_elements[2];
  var display_reactions_dict = node_elements[3];
  var entity_id_dict = node_elements[4];

  var _width = window.innerWidth;
  var _height = window.innerHeight - 75;

  make_graph(
    data,
    new_nodes,
    new_links,
    type_dict,
    node_dict,
    entity_id_dict,
    display_analytes_dict,
    display_reactions_dict,
    selector,
    _width,
    _height
  );
}

function run_value_connections() {

  let value_threshold = document.getElementById("conn_value_button").value;
  let perturbed_reactions = collect_perturbations(
    reaction_entity_dictionary,
    expression_dict,
    value_threshold,
    "value");
  show_graph(data, perturbed_reactions);
}

function run_stat_connections() {

  let stat_threshold = document.getElementById("conn_stat_button").value;
  let perturbed_reactions = collect_perturbations(
    reaction_entity_dictionary,
    stats_dict,
    stat_threshold,
    "value");
  show_graph(data, perturbed_reactions);
}

d3.select("#play_button_value").on("click", run_value_connections);
d3.select("#play_button_stat").on("click", run_stat_connections);

// Initial play
run_value_connections();
