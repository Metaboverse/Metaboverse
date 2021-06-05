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

// find global motifs at beginning, save as list, check each graph for members
function gatherMotifs(data, categories) {

  var excl_hubs = true;
  var hub_threshold = 50;

  var update_nodes = {};
  let n;
  for (n in data.nodes) {
    update_nodes[data.nodes[n].id] = data.nodes[n];
  }
  data.nodes = update_nodes;

  var update_links = {};
  let l;
  for (l in data.links) {
    let link_id = data.links[l].source + "," + data.links[l].target;
    update_links[link_id] = data.links[l];
  }
  data.links = update_links;

  let expression_dict = {};
  let stats_dict = {};
  let inferred_dict = {};
  for (let x in data.nodes) {
    let id = data.nodes[x]['id'];
    let expression = data.nodes[x]['values'];
    let stats = data.nodes[x]['stats'];
    expression_dict[id] = expression;
    stats_dict[id] = stats;
    inferred_dict[id] = data.nodes[x]['inferred']
  }

  let link_neighbors = {};
  for (let l in data.links) {
    let _source = data.links[l].source;
    let _target = data.links[l].target;
    
    if (data.nodes[_source].type === "reaction" 
    || data.nodes[_target].type === "reaction"
    || data.nodes[_source].type === "collapsed"
    || data.nodes[_target].type === "collapsed") {
    } else {
      if (!(_source in link_neighbors)) {
        link_neighbors[_source] = [];
      }
      link_neighbors[_source].push(_target);

      if (!(_target in link_neighbors)) {
        link_neighbors[_target] = [];
      }
      link_neighbors[_target].push(_source);
    }
  }

  let threshold = 1;

  let collapsed_motifs_Avg = motifSearch_Avg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_MaxMax = motifSearch_MaxMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_MinMin = motifSearch_MinMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_MaxMin = motifSearch_MaxMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_MinMax = motifSearch_MinMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_Sustained = motifSearch_Sustained(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_ModReg = modifierReg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let collapsed_motifs_ModTrans = modifierTransport(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)


  let motifs_Avg = motifSearch_Avg(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MaxMax = motifSearch_MaxMax(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MinMin = motifSearch_MinMin(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MaxMin = motifSearch_MaxMin(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MinMax = motifSearch_MinMax(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_Sustained = motifSearch_Sustained(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_ModReg = modifierReg(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_ModTrans = modifierTransport(
    threshold,
    data.reaction_dictionary,
    expression_dict,
    stats_dict,
    inferred_dict,
    link_neighbors,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)
  

  let all_collapsed_motifs = [];
  for (x in categories) {
    all_collapsed_motifs[x] = collapsed_motifs_Avg[x].concat(
      collapsed_motifs_MaxMax[x],
      collapsed_motifs_MinMin[x],
      collapsed_motifs_MaxMin[x],
      collapsed_motifs_MinMax[x],
      collapsed_motifs_Sustained[x],
      collapsed_motifs_ModReg[x],
      collapsed_motifs_ModTrans[x],
      motifs_Avg[x],
      motifs_MaxMax[x],
      motifs_MinMin[x],
      motifs_MaxMin[x],
      motifs_MinMax[x],
      motifs_Sustained[x],
      motifs_ModReg[x],
      motifs_ModTrans[x]);
  }

  let global_collapsed_motifs = [];
  for (x in categories) {
    global_collapsed_motifs[x] = [];
    for (let m in all_collapsed_motifs[x]) {
      global_collapsed_motifs[x].push(all_collapsed_motifs[x][m].id);
    }
  }

  let all_motifs = [];
  for (x in categories) {
    all_motifs[x] = motifs_Avg[x].concat(
      motifs_MaxMax[x],
      motifs_MinMin[x],
      motifs_MaxMin[x],
      motifs_MinMax[x],
      motifs_Sustained[x],
      motifs_ModReg[x],
      motifs_ModTrans[x]);
  }

  let global_motifs = [];
  for (x in categories) {
    global_motifs[x] = [];
    for (let m in all_motifs[x]) {
      global_motifs[x].push(all_motifs[x][m].id);
    }
  }


  return [global_collapsed_motifs, global_motifs];
}
