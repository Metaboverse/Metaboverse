
// search for motif 1
//let threshold = d3.select("#avg_num").node().value;
function motifSearch_Avg(
      threshold,
      collapsed_reaction_dict,
      expression_dict,
      path_mapper) {
  console.log("motif search 1")
  console.log("Avg threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let reactants = reaction.reactants;
    let products = reaction.products;
    let source_expression = [];
    let target_expression = [];

    reactants.forEach(l=>{
      let reactant_expr = expression_dict[l];
      if(reactant_expr !== null){
        source_expression.push(parseFloat(reactant_expr));
      }
    })

    products.forEach(l=>{
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

    if(updated_source.length>0 && updated_target.length>0){
      let source_avg = computeAvg(updated_source);
      let target_avg = computeAvg(updated_target);

      if(Math.abs(source_avg - target_avg)>=threshold){
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

function computeAvg(arr){
  let arr_sum = arr[0];
  for(let i=1; i<arr.length; i++){
    arr_sum += arr[i];
  }
  let arr_avg = arr_sum / arr.length;
  return arr_avg;
}

// MaxMax
//let threshold = d3.select("#maxmax_num").node().value;
function motifSearch_MaxMax(
    threshold,
    collapsed_reaction_dict,
    expression_dict,
    path_mapper) {
  console.log("motif search 2")
  console.log("MaxMax threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let reactants = reaction.reactants;
    let products = reaction.products;
    let source_expression = [];
    let target_expression = [];

    reactants.forEach(l=>{
      let reactant_expr = expression_dict[l];
      if(reactant_expr !== null){
        source_expression.push(parseFloat(reactant_expr));
      }
    })

    products.forEach(l=>{
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

    if(updated_source.length>0 && updated_target.length>0){
      let source_max = Math.max(...updated_source);
      let target_max = Math.max(...updated_target);
      if(Math.abs(source_max - target_max)>=threshold){
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
    path_mapper) {
  console.log("motif search 3")
  console.log("MaxMin threshold set at: ", threshold)
  let discovered_motifs = [];

  for(let rxn in collapsed_reaction_dict){
    let reaction = collapsed_reaction_dict[rxn];
    let reactants = reaction.reactants;
    let products = reaction.products;
    let source_expression = [];
    let target_expression = [];

    reactants.forEach(l=>{
      let reactant_expr = expression_dict[l];
      if(reactant_expr !== null){
        source_expression.push(parseFloat(reactant_expr));
      }
    })

    products.forEach(l=>{
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
    if(updated_source.length>0 && updated_target.length>0){
      let source_max = Math.max(...updated_source);
      let target_min = Math.min(...updated_target);
      if(Math.abs(source_max - target_min)>=threshold){
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
