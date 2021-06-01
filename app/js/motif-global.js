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

  let expression_dict = {};
  let stats_dict = {};
  for (let x in data.nodes) {
    let id = data.nodes[x]['id'];
    let expression = data.nodes[x]['values'];
    let stats = data.nodes[x]['stats'];
    expression_dict[id] = expression;
    stats_dict[id] = stats;
  }

  let threshold = 1;
  let motifs_Avg = motifSearch_Avg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MaxMax = motifSearch_MaxMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MinMin = motifSearch_MinMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MaxMin = motifSearch_MaxMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_MinMax = motifSearch_MinMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_Sustained = motifSearch_Sustained(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_ModReg = modifierReg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

  let motifs_ModTrans = modifierTransport(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    stats_dict,
    data.motif_reaction_dictionary,
    data.degree_dictionary,
    data.blocklist,
    categories)

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

  let global_motifs_reactions = [];
  for (x in categories) {
    global_motifs_reactions[x] = [];
    for (let m in all_motifs[x]) {
      global_motifs_reactions[x].push(all_motifs[x][m].id);
    }
  }

  return global_motifs_reactions;
}
