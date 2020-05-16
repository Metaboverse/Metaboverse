/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019  Youjia Zhou, Jordan A. Berg
  zhou325 <at> sci <dot> utah <dot> edu
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

var eval_modifiers = false;

function modifiersChecked() {
  if (eval_modifiers === false) {
    eval_modifiers = true;
  } else {
    eval_modifiers = false;
  }
  console.log("Motif evaluation includes modifiers: ", eval_modifiers)
}

function parseComponents(
    reaction,
    expression_dict,
    sample_index) {

  let reactants = reaction.reactants;
  let products = reaction.products;
  let modifiers = reaction.modifiers;

  let source_expression = [];
  let target_expression = [];

  let temp_reactants = $.extend(true, [], reactants);
  let temp_products = $.extend(true, [], products);

  if (eval_modifiers === true) {
    for (x in modifiers) {
      if (modifiers[x][1] === "catalyst") {
        temp_reactants.push(modifiers[x][0])
      } else if (modifiers[x][1] === "inhibitor") {
        temp_reactants.push(modifiers[x][0])
      } else {}
    }
  }

  temp_reactants.forEach(l=>{
    let reactant_expr = expression_dict[l][sample_index];
    if(reactant_expr !== null){
      source_expression.push(parseFloat(reactant_expr));
    }
  })

  temp_products.forEach(l=>{
    let product_expr = expression_dict[l][sample_index];
    if(product_expr !== null){
      target_expression.push(parseFloat(product_expr));
    }
  })

  let updated_source = source_expression.filter(function (value) {
      return !Number.isNaN(value);
  });
  let updated_target = target_expression.filter(function (value) {
      return !Number.isNaN(value);
  });

  return [updated_source, updated_target]
}

function computeAvg(arr){
  let arr_sum = arr[0];
  for(let i=1; i<arr.length; i++){
    arr_sum += arr[i];
  }
  let arr_avg = arr_sum / arr.length;
  return arr_avg;
}

// search for motif 1
//let threshold = d3.select("#avg_num").node().value;
function motifSearch_Avg(
      threshold,
      collapsed_reaction_dict,
      expression_dict,
      path_mapper,
      sample_indices) {
  console.log("motif search 1")
  console.log("Avg threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    for(let rxn in collapsed_reaction_dict) {
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];

      if(updated_source.length>0 && updated_target.length>0){
        let source_avg = computeAvg(updated_source);
        let target_avg = computeAvg(updated_target);

        if(Math.abs(source_avg - target_avg)>=threshold){
          sample_motifs.push(reaction);
        }
      }
    }
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  console.log(discovered_motifs);
  return discovered_motifs;
}

// MaxMax
//let threshold = d3.select("#maxmax_num").node().value;
function motifSearch_MaxMax(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    path_mapper,
    sample_indices) {
  console.log("motif search 2")
  console.log("MaxMax threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    for(let rxn in collapsed_reaction_dict){
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];

      if(updated_source.length>0 && updated_target.length>0){
        let source_max = Math.max(...updated_source);
        let target_max = Math.max(...updated_target);
        if(Math.abs(source_max - target_max)>=threshold){
          sample_motifs.push(reaction);
        }
      }
    }
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  console.log(discovered_motifs);
  return discovered_motifs;
}

//MaxMin
//let threshold = d3.select("#maxmin_num").node().value;
function motifSearch_MaxMin(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    path_mapper,
    sample_indices) {
  console.log("motif search 3")
  console.log("MaxMin threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    for(let rxn in collapsed_reaction_dict){
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];

      if(updated_source.length>0 && updated_target.length>0){
        let source_max = Math.max(...updated_source);
        let target_min = Math.min(...updated_target);
        if(Math.abs(source_max - target_min)>=threshold){
          sample_motifs.push(reaction);
        }
      }
    }
    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  console.log(discovered_motifs);
  return discovered_motifs;
}

//Sustained
//let threshold = d3.select("#sustained_num").node().value;
function motifSearch_Sustained(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    path_mapper,
    sample_indices) {
  console.log("motif search 4")
  console.log("Sustained perturbation threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    for(let rxn in collapsed_reaction_dict){
      let reaction = collapsed_reaction_dict[rxn];
      let comps = parseComponents(
        reaction,
        expression_dict,
        _idx)
      let updated_source = comps[0];
      let updated_target = comps[1];

      if(updated_source.length>0 && updated_target.length>0) {

        // Sustained up-regulation
        let up_in = false;
        let up_out = false;
        for (i in updated_source) {
          if (updated_source[i] >= threshold) {
            up_in = true;
          }
        }
        for (j in updated_target) {
          if (updated_target[j] >= threshold) {
            up_out = true;
          }
        }

        // Sustained down-regulation
        let down_in = false;
        let down_out = false;
        for (k in updated_source) {
          if (updated_source[k] <= -(threshold)) {
            down_in = true;
          }
        }
        for (l in updated_target) {
          if (updated_target[l] <= -(threshold)) {
            down_out = true;
          }
        }

        if (((down_in === true) & (down_out === true)) | ((up_in === true) & (up_out === true))) {
          sample_motifs.push(reaction);
        }
      }
    }

    for (let m in sample_motifs) {
      sample_motifs[m]['pathways'] = path_mapper[sample_motifs[m]['id']]
    }
    discovered_motifs.push(sample_motifs);
  }
  console.log(discovered_motifs);
  return discovered_motifs;
}

//Path max/min comparison
//let threshold = d3.select("#pathmax_num").node().value;
function motifSearch_PathMax(
    threshold,
    mod_collapsed_pathways,
    collapsed_reaction_dict,
    expression_dict,
    path_mapper,
    sample_indices) {
  console.log("motif search 5")
  console.log("Pathway min/max threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    // For each pathway, get reactions
    for (pathway in mod_collapsed_pathways) {

      let values = [];
      let reactions = mod_collapsed_pathways[pathway].reactions;
      for (rxn in reactions) {
        let reaction = collapsed_reaction_dict[reactions[rxn]];
        let comps = parseComponents(
          reaction,
          expression_dict,
          _idx)
        let updated_source = comps[0];
        let updated_target = comps[1];

        // Combine all values
        values = values.concat(updated_source, updated_target);
      }
      if (values.length > 0) {
        if (Math.abs(Math.max.apply(Math,values) - Math.min.apply(Math,values)) >= threshold) {
          sample_motifs.push(mod_collapsed_pathways[pathway]);
        }
      }
    }
    discovered_motifs.push(sample_motifs);
  }
  return discovered_motifs;
  // Get expression for all entities of components
  // Compare min/max
  // Return to list if true
  // Will need to reformat motif display since just showing pathways, not reactions (make dummy reactions?)
  // Make sorting index for later that is also output

}

//Path coverage comparison
//let threshold = d3.select("#pathcov_num").node().value;
function motifSearch_PathCov(
    threshold,
    min_coverage,
    mod_collapsed_pathways,
    collapsed_reaction_dict,
    stats_dict,
    path_mapper,
    sample_indices) {
  console.log("motif search 6")
  console.log("threshold set at: ", threshold)
  let discovered_motifs = [];

  for (_idx in sample_indices) {
    let sample_motifs = [];

    // For each pathway, get reactions
    for (pathway in mod_collapsed_pathways) {

      let values = 0;
      let total = 0;

      let reactions = mod_collapsed_pathways[pathway].reactions;
      for (rxn in reactions) {
        let reaction = collapsed_reaction_dict[reactions[rxn]];
        let comps = parseComponents(
          reaction,
          stats_dict,
          _idx)
        let updated_source = comps[0].filter(stat => stat < threshold);
        let updated_target = comps[1].filter(stat => stat < threshold);

        // Check that at least one component in the reaction meets thresholding
        // criteria
        if (updated_source.length + updated_target.length > 0) {
          values = values + 1;
        }
        total = total + 1;
      }

      let cov = values / total;
      if (cov >= min_coverage) {
        sample_motifs.push([
          mod_collapsed_pathways[pathway],
          cov,
          values,
          total
        ]);
      }
    }

    sample_motifs.sort(function(one, two) {
      return two[1] - one[1];
    });

    discovered_motifs[_idx] = sample_motifs;
  }
  console.log(discovered_motifs);
  return discovered_motifs;
}
