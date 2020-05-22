
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
    expression_dict) {

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
    let reactant_expr = expression_dict[l];
    if(reactant_expr !== null){
      source_expression.push(parseFloat(reactant_expr));
    }
  })

  temp_products.forEach(l=>{
    let product_expr = expression_dict[l];
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
      value_type,
      modifiers) {
  console.log("motif search 1")
  console.log("Avg threshold set at: ", threshold)
  console.log(value_type)
  console.log(expression_dict)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let comps = parseComponents(
      reaction,
      expression_dict)
    let updated_source = comps[0];
    let updated_target = comps[1];

    if(updated_source.length>0 && updated_target.length>0){
      let source_avg = computeAvg(updated_source);
      let target_avg = computeAvg(updated_target);

      if(Math.abs(source_avg - target_avg)>=threshold){
        reaction.magnitude_change = Math.abs(source_avg - target_avg);
        discovered_motifs.push(reaction);
      }
    }
  }
  for (let m in discovered_motifs) {
    discovered_motifs[m]['pathways'] = path_mapper[discovered_motifs[m]['id']]
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
    modifiers) {
  console.log("motif search 2")
  console.log("MaxMax threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let comps = parseComponents(
      reaction,
      expression_dict)
    let updated_source = comps[0];
    let updated_target = comps[1];

    if(updated_source.length>0 && updated_target.length>0){
      let source_max = Math.max(...updated_source);
      let target_max = Math.max(...updated_target);
      if(Math.abs(source_max - target_max)>=threshold){
        reaction.magnitude_change = Math.abs(source_max - target_max);
        discovered_motifs.push(reaction);
      }
    }
  }
  for (let m in discovered_motifs) {
    discovered_motifs[m]['pathways'] = path_mapper[discovered_motifs[m]['id']]
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
    modifiers) {
  console.log("motif search 3")
  console.log("MaxMin threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let comps = parseComponents(
      reaction,
      expression_dict)
    let updated_source = comps[0];
    let updated_target = comps[1];

    if(updated_source.length>0 && updated_target.length>0){
      let source_max = Math.max(...updated_source);
      let target_min = Math.min(...updated_target);
      if(Math.abs(source_max - target_min)>=threshold){
        reaction.magnitude_change = Math.abs(source_max - target_min);
        discovered_motifs.push(reaction);
      }
    }
  }
  for (let m in discovered_motifs) {
    discovered_motifs[m]['pathways'] = path_mapper[discovered_motifs[m]['id']]
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
    modifiers) {
  console.log("motif search 4")
  console.log("Sustained perturbation threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let comps = parseComponents(
      reaction,
      expression_dict)
    let updated_source = comps[0];
    let updated_target = comps[1];

    if(updated_source.length>0 && updated_target.length>0) {

      // Sustained up-regulation
      let up_in = false;
      let up_out = false;
      let magnitude_change_up;
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
      if(up_in===true && up_out===true) {
        magnitude_change_up = Math.abs(Math.max(...updated_source) - Math.max(...updated_target));
      }

      // Sustained down-regulation
      let down_in = false;
      let down_out = false;
      let magnitude_change_down;
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
      if(down_in===true && down_out===true) {
        magnitude_change_down = Math.abs(Math.min(...updated_source) - Math.min(...updated_target));
      }

      if (((down_in === true) & (down_out === true)) | ((up_in === true) & (up_out === true))) {
        let magnitude_change;
        if(magnitude_change_up && magnitude_change_down){
          magnitude_change = Math.max(magnitude_change_up, magnitude_change_down);
        } else if (magnitude_change_up){
          magnitude_change = magnitude_change_up;
        } else if(magnitude_change_down){
          magnitude_change = magnitude_change_down;
        }
        reaction.magnitude_change = magnitude_change;
        discovered_motifs.push(reaction);
      }
    }
  }

  for (let m in discovered_motifs) {
    discovered_motifs[m]['pathways'] = path_mapper[discovered_motifs[m]['id']]
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
    modifiers) {
  console.log("motif search 5")
  console.log("Pathway min/max threshold set at: ", threshold)
  let discovered_motifs = [];

  // For each pathway, get reactions
  for (pathway in mod_collapsed_pathways) {

    let values = [];

    let reactions = mod_collapsed_pathways[pathway].reactions;
    for (rxn in reactions) {
      let reaction = collapsed_reaction_dict[reactions[rxn]];
      let comps = parseComponents(
        reaction,
        expression_dict)
      let updated_source = comps[0];
      let updated_target = comps[1];

      // Combine all values
      values = values.concat(updated_source, updated_target);
    }

    if (values.length > 0) {
      if (Math.abs(Math.max.apply(Math,values) - Math.min.apply(Math,values)) >= threshold) {
        mod_collapsed_pathways[pathway].magnitude_change = Math.abs(Math.max.apply(Math,values) - Math.min.apply(Math,values));
        discovered_motifs.push(mod_collapsed_pathways[pathway]);
      }
    }
  }

  console.log(discovered_motifs);
  return discovered_motifs;
  // Get expression for all entities of components
  // Compare min/max
  // Return to list if true
  // Will need to reformat motif display since just showing pathways, not reactions (make dummy reactions?)
  // Make sorting index for later that is also output

}
