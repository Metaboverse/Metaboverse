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

// find global motifs at beginning, save as list, check each graph for members
function gatherMotifs(data, categories) {

  var excl_hubs = false;
  var hub_threshold = 100;

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
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)

  let motifs_MaxMax = motifSearch_MaxMax(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)

  let motifs_MaxMin = motifSearch_MaxMin(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)

  let motifs_Sustained = motifSearch_Sustained(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)

  let motifs_ModReg = modifierReg(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)

  let motifs_ModTrans = modifierTransport(
    threshold,
    data.collapsed_reaction_dictionary,
    expression_dict,
    data.motif_reaction_dictionary,
    data.degree_dict,
    categories)


  let all_motifs = [];
  for (x in categories) {
    all_motifs[x] = motifs_Avg[x].concat(
      motifs_MaxMax[x],
      motifs_MaxMin[x],
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
