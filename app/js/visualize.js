/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/j-berg/Metaboverse/
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
const d3 = require("d3");
var fs = require("fs");
var saveSVG = require("save-svg-as-png");

const max_nodes = 1500;
var sample = 0;
var entity = "values_js";
var graph_genes = true;
var collapse_reactions = true;
var saved_nodes = [];
var saved_links = [];
var collapsed_nodes = [];
var collapsed_links = [];

// Change user selection based on input
var selection = null;
function selectPathway() {
  var selection = document.getElementById("pathwayMenu").value;
  console.log("Selection: " + selection);
  return selection;
}

var superSelection = null;
function selectSuperPathway(selector) {
  var superSelection = document.getElementById("superPathwayMenu").value;
  console.log("Super Pathway: " + superSelection);
  return superSelection;
}

// Populate dictionary to access component reactions for each pathway
function make_pathway_dictionary(data, database_key) {
  // Get pathway name and ID
  var pathways = data[database_key];
  var pathway_dict = {};
  for (var key in pathways) {
    pathway_dict[pathways[key]["name"]] = {
      id: pathways[key]["name"],
      reactome: pathways[key]["reactome"],
      reactions: pathways[key]["reactions"]
    };
  }

  return pathway_dict;
}

// Populate dictionary to access component reactions for each super-pathway
function make_superPathway_dictionary(data) {
  // Get pathway name and ID
  var superPathways = data.super_pathways;
  var superPathwayDict = {};
  for (var key in superPathways) {
    superPathwayDict[superPathways[key]["name"]] = {
      id: superPathways[key]["name"],
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

function emptyMenu(selectMenu) {
  var i;
  for (i = selectMenu.options.length - 1; i >= 0; i--) {
    selectMenu.remove(i);
  }
}

function parsePathways(pathway_dict, selectedReactions) {
  var parsed_pathway_dict = {};
  // from https://stackoverflow.com/a/53606357
  let checker = (arr, target) => target.every(v => arr.includes(v));

  for (pathway in pathway_dict) {
    belongs = checker(selectedReactions, pathway_dict[pathway]["reactions"]);

    if (belongs === true) {
      parsed_pathway_dict[pathway] = pathway_dict[pathway];
    }
  }

  return parsed_pathway_dict;
}

function initialize_nodes(nodes, node_dict, type_dict) {
  var expression_dict = {};
  var stats_dict = {};
  var display_analytes_dict = {};
  var display_reactions_dict = {};
  var entity_id_dict = {};

  // Make dictionary of node color values and types
  nodes.forEach(function(node) {
    node_dict[node["name"]] = node["values_js"];
    type_dict[node["name"]] = node["type"];

    try {
      expression_dict[node["name"]] = node["values"][0];
      stats_dict[node["name"]] = node["stats"][0];
    }
    catch (err) {
      console.log(node);
    }

    entity_id_dict[node["name"]] = node["id"];
    entity_id_dict[node["id"]] = node["name"];

    if (node["type"] === "reaction") {
      display_analytes_dict[node["name"]] = "none";
      display_reactions_dict[node["name"]] = "inline";
    } else {
      display_reactions_dict[node["name"]] = "none";
      display_analytes_dict[node["name"]] = "inline";
    }
  });

  return [
    node_dict,
    type_dict,
    expression_dict,
    stats_dict,
    display_analytes_dict,
    display_reactions_dict,
    entity_id_dict
  ];
}

function linkArc(d) {
  var dx = d.target.x - d.source.x;
  var dy = d.target.y - d.source.y;
  var dr = Math.sqrt(dx * dx + dy * dy);

  return (
    "M" +
    d.source.x +
    "," +
    d.source.y +
    "A" +
    dr +
    "," +
    dr +
    " 0 0,1 " +
    d.target.x +
    "," +
    d.target.y
  );
}

function transform(d) {
  return "translate(" + d.x + "," + d.y + ")";
}

function parse_pathway(data, reactions) {
  var reaction_dictionary = data.reaction_dictionary;

  // Parse through each reaction listed and get the component parts
  var components = [];
  for (rxn in reactions) {
    var target_rxns = reaction_dictionary[reactions[rxn]];
    components.push(reactions[rxn]);
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
  var new_elements = get_nodes_links(data, components);

  return new_elements;
}

function get_nodes_links(data, components) {
  var nodes = data["nodes"];
  var links = data["links"];

  // Parse the nodes of interest
  var add_nodes = [];
  nodes.forEach(function(node) {
    if (components.includes(node["id"])) {
      var node_copy = $.extend(true, {}, node);
      add_nodes.push(node_copy);
    }
  });

  var add_links = [];
  // Parse out links of interest
  links.forEach(function(link) {
    if (
      components.includes(link.source) &&
      components.includes(link.target) &&
      link.source !== link.target
    ) {
      var link_copy = $.extend(true, {}, link);
      add_links.push(link_copy);
    }
  });

  return [add_nodes, add_links];
}

function nearest_neighbors(data, entity_id) {
  // Get current nearest neighbors value
  var kNN = document.getElementById("kNN_button").value;
  document.getElementById("pathwayMenu").value = "";

  // Curate kNN to node of interest
  var nn_elements = parse_kNN_pathway(data, entity_id, kNN);
  var new_nodes = nn_elements[0];
  var new_links = nn_elements[1];

  // Save old graph in burner variable

  // Recurate the graph
  // Initialize variables
  var node_dict = {};
  var type_dict = {};

  var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
  var node_dict = node_elements[0];
  var type_dict = node_elements[1];
  var expression_dict = node_elements[2];
  var stats_dict = node_elements[3];
  var display_analytes_dict = node_elements[4];
  var display_reactions_dict = node_elements[5];
  var entity_id_dict = node_elements[6];

  // Remove old plot and plot this one
  make_graph(
    data,
    new_nodes,
    new_links,
    type_dict,
    node_dict,
    entity_id_dict,
    expression_dict,
    stats_dict,
    display_analytes_dict,
    display_reactions_dict
  );
}

function parse_kNN_pathway(data, entity_id, kNN) {
  // Reset warning messages
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

  var reaction_dictionary = data.reaction_dictionary;

  // Parse through each reaction where entity is a component
  nn_components = [entity_id];
  for (reaction in reaction_dictionary) {
    //Return all reactions that contain the entity
    if (
      reaction_dictionary[reaction]["reactants"].includes(entity_id) ||
      reaction_dictionary[reaction]["products"].includes(entity_id) ||
      reaction_dictionary[reaction]["modifiers"].includes(entity_id) ||
      reaction_dictionary[reaction]["additional_components"].includes(entity_id)
    ) {
      nn_components.push(reaction_dictionary[reaction]["id"]);
      for (x in reaction_dictionary[reaction]["reactants"]) {
        nn_components.push(reaction_dictionary[reaction]["reactants"][x]);
      }
      for (x in reaction_dictionary[reaction]["products"]) {
        nn_components.push(reaction_dictionary[reaction]["products"][x]);
      }
      for (x in reaction_dictionary[reaction]["modifiers"]) {
        nn_components.push(reaction_dictionary[reaction]["modifiers"][x][0]);
      }
      for (x in reaction_dictionary[reaction]["additional_components"]) {
        nn_components.push(
          reaction_dictionary[reaction]["additional_components"][x]
        );
      }
    }
  }

  var nodes = data["nodes"];
  nodes.forEach(function(node) {
    if (node.degree > 500) {
      console.log(
        "filtering " +
          node.name +
          " (" +
          node.degree +
          ") -- may cause edge loss"
      );
      nn_components = nn_components.filter(x => x !== node.id);
    }
  });

  // If too many nodes for first neighborhood, stop plotting
  var escape = nn_components.length;
  if (escape > max_nodes) {
    document.getElementById("warning_line_1").innerHTML =
      '<i style="color: red;">Too many entities to plot</i><br><i style="color: red;">Will not plot</i>';
    document.getElementById("warning_line_2").innerHTML = "<br>";
    kNN = 0;
    nn_components = [];
  }

  if (kNN > 1) {
    document.getElementById("warning_line_1").innerHTML =
      '<i style="color: red;">Please wait<br><br></i>';
    document.getElementById("warning_line_2").innerHTML = "";

    // Filter out any components that are above hub threshold for kNN
    var hub_limit = document.getElementById("hub_button").value;
    var hub_exclusion = new Set();
    if (hub_limit > 0) {
      nodes.forEach(function(node) {
        if (node.degree > hub_limit) {
          hub_exclusion.add(node.id);
        }
      });
    }
    console.log(hub_exclusion);

    n = 1;
    while (n < kNN) {
      var components = [];

      for (element in nn_components) {
        for (reaction in reaction_dictionary) {
          //Return all reactions that contain the entity
          if (
            (reaction_dictionary[reaction]["reactants"].includes(
              nn_components[element]
            ) ||
              reaction_dictionary[reaction]["products"].includes(
                nn_components[element]
              ) ||
              reaction_dictionary[reaction]["modifiers"].includes(
                nn_components[element]
              ) ||
              reaction_dictionary[reaction]["additional_components"].includes(
                nn_components[element]
              )) &&
            !hub_exclusion.has(reaction_dictionary[reaction]["id"])
          ) {
            components.push(reaction_dictionary[reaction]["id"]);

            for (x in reaction_dictionary[reaction]["reactants"]) {
              if (
                !hub_exclusion.has(
                  reaction_dictionary[reaction]["reactants"][x]
                )
              ) {
                components.push(reaction_dictionary[reaction]["reactants"][x]);
              }
            }
            for (x in reaction_dictionary[reaction]["products"]) {
              if (
                !hub_exclusion.has(reaction_dictionary[reaction]["products"][x])
              ) {
                components.push(reaction_dictionary[reaction]["products"][x]);
              }
            }
            for (x in reaction_dictionary[reaction]["modifiers"]) {
              if (
                !hub_exclusion.has(
                  reaction_dictionary[reaction]["modifiers"][x][0]
                )
              ) {
                components.push(
                  reaction_dictionary[reaction]["modifiers"][x][0]
                );
              }
            }
            for (x in reaction_dictionary[reaction]["additional_components"]) {
              if (
                !hub_exclusion.has(
                  reaction_dictionary[reaction]["additional_components"][x]
                )
              ) {
                components.push(
                  reaction_dictionary[reaction]["additional_components"][x]
                );
              }
            }
          }
        }
      }

      // Only combine the lists if they pass the threshold
      escape = escape + components.length;

      if (escape > max_nodes) {
        console.log(escape);
        document.getElementById("warning_line_1").innerHTML =
          '<i style="color: red;">Too many entities to plot. Will only plot first ' +
          n +
          " neighborhood(s)</i>";
        document.getElementById("warning_line_2").innerHTML = "";
        n = kNN + 2;
      } else {
        nn_components = nn_components.concat(components);
        document.getElementById("warning_line_1").innerHTML = "<br>";
        document.getElementById("warning_line_2").innerHTML = "<br><br>";
      }
      n = n + 1;
    }
  }

  var new_elements = get_nodes_links(data, nn_components);
  return new_elements;
}

function kNN_input(d) {
  var knn_value = document.getElementById("kNN_button").value;
  console.log("k-NN parameter now set to: " + knn_value);
}

function hub_input(d) {
  var hub_value = document.getElementById("hub_button").value;
  console.log("Hub threshold now set to: " + hub_value);
}

function get_link(d) {
  if (d.type === "complex_component") {
    if (
      d.sub_type === "metabolite_component" ||
      d.sub_type === "protein_component" ||
      d.sub_type === "gene_component"
    ) {
      return d.sub_type;
    } else {
      return d.type;
    }
  } else {
    return d.type;
  }
}

function make_graph(
  data,
  new_nodes,
  new_links,
  type_dict,
  node_dict,
  entity_id_dict,
  expression_dict,
  stats_dict,
  display_analytes_dict,
  display_reactions_dict
) {
  console.log(new_nodes);
  console.log(new_links);

  // Allow flexible window dimensions based on initial window size when opened
  var width = window.innerWidth;
  var height = window.innerHeight - 75;

  // Restart graph
  d3.select("svg").remove();

  // Initialize force graph object
  var svg = d3
    .select("#graph")
    .append("svg")
    .attr("width", width - 5)
    .attr("height", height - 80)
    .call(
      d3.zoom().on("zoom", function() {
        svg.attr("transform", d3.event.transform);
      })
    )
    .on("dblclick.zoom", null)
    .append("g");

  const forceX = d3.forceX(width / 2).strength(0.015);
  const forceY = d3.forceY(height / 2).strength(0.015);

  const simulation = d3
    .forceSimulation(new_nodes)
    .force(
      "link",
      d3
        .forceLink(new_links)
        .id(d => d.id)
        .distance(40)
        .strength(1)
    )
    .force("charge", d3.forceManyBody().strength(-1000))
    .force("center", d3.forceCenter(width / 2, height / 2));

  simulation
    .alphaTarget(0.01)
    .alphaMin(0.1)
    .velocityDecay(0.7);

  // Generate edges with style attributes
  svg
    .append("defs")
    .selectAll("marker")
    .data([
      "collapsed",
      "reaction",
      "reactant",
      "product",
      "inhibitor",
      "catalyst",
      "gene_component",
      "complex_component",
      "mirna_component",
      "other",
      "protein_component",
      "metabolite_component"
    ])
    .enter()
    .append("marker")
    .attr("id", function(d) {
      return d;
    })
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 15.7)
    .attr("refY", -0.18)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0, -5L10, 0L0, 5");

  var link = svg
    .append("g")
    .selectAll("path")
    .data(new_links)
    .enter()
    .append("path")
    .attr("class", function(d) {
      return "link " + get_link(d);
    })
    .attr("marker-end", function(d) {
      return "url(#" + get_link(d) + ")";
    });

  var node = svg
    .selectAll(".node")
    .data(new_nodes)
    .enter()
    .append("g")
    .attr("class", "node")
    .style("--node_color", function(d) {
      return "rgba(" + d[entity].toString() + ")";
    })
    .style("r", function(d) {
      return 6;
    })
    .style("stroke-dasharray", function(d) {
      if (d.inferred === "true" || d.type === "collapsed") {
        return "2,2";
      } else {
        return "none";
      }
    })
    .call(
      d3
        .drag()
        .subject(dragsubject)
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended)
    );

  var circle = node.append(function(d) {
    if (d.type === "gene_component" || d.sub_type === "gene") {
      return document.createElementNS("http://www.w3.org/2000/svg", "ellipse");
    } else if (d.complex === "true" || d.sub_type === "protein_component") {
      return document.createElementNS("http://www.w3.org/2000/svg", "rect");
    } else {
      return document.createElementNS("http://www.w3.org/2000/svg", "circle");
    }
  });

  circle.on("dblclick", function(d) {
    document.getElementById("reaction_notes").innerHTML = "";

    if (
      type_dict[d["name"]] === "reaction" ||
      type_dict[d["name"]] === "collapsed"
    ) {
      console.log("Selected a reaction, will not perform kNN graphing");
      document.getElementById("reaction_notes").innerHTML =
        "<b><i>" + d.name + "</i></b>: " + d.notes;
    } else {
      document.getElementById("type_selection_type").innerHTML =
        "Nearest Neighbor";

      var mod_selection = determineWidth(d["name"]);
      document.getElementById("type_selection").innerHTML = mod_selection;

      graph_genes = true;
      nearest_neighbors(data, entity_id_dict[d["name"]]);
    }
  });

  var text = node.append("text").html(function(d) {
    if (type_dict[d.name] === "reaction" ||
      type_dict[d.name] === "collapsed"
    ) {
      // Label other nodes with expression value in parentheses
      return (
        "<tspan dx='16' y='.31em' style='font-weight: bold;'>" +
        d.name +
        "</tspan>"
      );
    } else {
      return (
        "<tspan dx='16' y='-.5em' style='font-weight: bold;'>" +
        d.name +
        "</tspan>" +
        "<tspan x='16' y='.7em'>Value: " +
        parseFloat(expression_dict[d.name]).toFixed(2) +
        "</tspan>" +
        "<tspan x='16' y='1.7em'>Statistic: " +
        parseFloat(stats_dict[d.name]).toFixed(2) +
        "</tspan>"
      );
    }
  });

  simulation.on("tick", tick);

  // Refresh current graph
  d3.select("#restartButton").on("click", function() {
    simulation
      .alphaTarget(0.01)
      .alphaMin(0.1)
      .velocityDecay(0.7)
      .restart();
  });

  d3.select("#saveGraph").on("click", function() {
    saveSVG.saveSvgAsPng(
      document.getElementsByTagName("svg")[0],
      "plot.png",
      (encoderOptions = 1),
      (scale = 5),
      (encoderType = "image/png")
    );
  });

  d3.select("#openPathway").on("click", function() {
    pathway = selectPathway();

    if (pathway !== "") {
      pathway_id = pathway_dict[pathway].reactome;
      reactome_string = "https://reactome.org/PathwayBrowser/#/" + pathway_id;
      window.open(reactome_string, "window name", "window settings");
    }
  });

  toggle_e = true;
  d3.select("#toggleExpression").on("click", function() {
    if (toggle_e === false) {
      toggle_e = true;
      text.html(function(d) {
        if (type_dict[d.name] === "reaction") {
          // If reaction node, do not display expression value
          return (
            "<tspan dx='16' y='.31em' style='font-weight: bold;'>" +
            d.name +
            "</tspan>"
          );
        } else {
          // Label other nodes with expression value in parentheses
          return (
            "<tspan dx='16' y='-.5em' style='font-weight: bold;'>" +
            d.name +
            "</tspan>" +
            "<tspan x='16' y='.7em'>Value: " +
            parseFloat(expression_dict[d.name]).toFixed(2) +
            "</tspan>" +
            "<tspan x='16' y='1.7em'>Statistic: " +
            parseFloat(stats_dict[d.name]).toFixed(2) +
            "</tspan>"
          );
        }
      });
    } else {
      toggle_e = false;
      text.html(function(d) {
        return (
          "<tspan dx='16' y='.31em' style='font-weight: bold;'>" +
          d.name +
          "</tspan>"
        );
      });
    }
  });

  toggle_a = false;
  toggle_r = false;
  text.style("--node_display", function(d) {
    return "none";
  });

  d3.select("#toggleAnalytes").on("click", function() {
    if (toggle_a === false) {
      toggle_a = true;
      determine_displays(toggle_a, toggle_r);
    } else {
      toggle_a = false;
      determine_displays(toggle_a, toggle_r);
    }
  });

  d3.select("#toggleReactions").on("click", function() {
    if (toggle_r === false) {
      toggle_r = true;
      determine_displays(toggle_a, toggle_r);
    } else {
      toggle_r = false;
      determine_displays(toggle_a, toggle_r);
    }
  });

  $("#toggleColors").click(function(e) {
    if (toggle_c === false) {
      entity = "stats_js";
      node.style("--node_color", function(d) {
        return "rgba(" + d[entity].toString() + ")";
      });
    } else {
      entity = "values_js";
      node.style("--node_color", function(d) {
        return "rgba(" + d[entity].toString() + ")";
      });
    }
  });

  function determine_displays(toggle_a, toggle_r) {
    if (toggle_a === false && toggle_r === false) {
      var display_labels = "none";
    } else if (toggle_a === false && toggle_r === true) {
      var display_labels = "reactions";
    } else if (toggle_a === true && toggle_r === false) {
      var display_labels = "analytes";
    } else if (toggle_a === true && toggle_r === true) {
      var display_labels = "all";
    } else {
      var display_labels = "none";
    }

    // options:

    // none -> all labels hidden until hovered
    if (display_labels === "none") {
      text.style("--node_display", function(d) {
        return "none";
      });
    }
    // analytes -> analytes shown, reactions hovered
    else if (display_labels === "analytes") {
      text.style("--node_display", function(d) {
        return display_analytes_dict[d.name];
      });
    }
    // reactions -> reactions shown, analytes hovered
    else if (display_labels === "reactions") {
      text.style("--node_display", function(d) {
        return display_reactions_dict[d.name];
      });
    }
    // analytes + reactions -> all -> all shown
    else if (display_labels === "all") {
      text.style("--node_display", function(d) {
        return "inline";
      });
    } else {
      text.style("--node_display", function(d) {
        return "none";
      });
    }
  }

  d3.select("#toggleGenes").on("click", function() {
    if (graph_genes === false) {
      graph_genes = true;
      new_nodes = saved_nodes;
      new_links = saved_links;

      make_graph(
        data,
        new_nodes,
        new_links,
        type_dict,
        node_dict,
        entity_id_dict,
        expression_dict,
        stats_dict,
        display_analytes_dict,
        display_reactions_dict
      );
    } else {
      graph_genes = false;
      saved_nodes = new_nodes;
      saved_links = new_links;

      var newer_components = [];
      for (x in new_nodes) {
        if (
          new_nodes[x]["type"] !== "gene_component" &&
          new_nodes[x]["sub_type"] !== "gene_component"
        ) {
          newer_components.push(new_nodes[x]["id"]);
        }
      }

      var newer_elements = get_nodes_links(data, newer_components);
      var newer_nodes = newer_elements[0];
      var newer_links = newer_elements[1];

      // Initialize variables
      var new_node_dict = {};
      var new_type_dict = {};

      var new_node_elements = initialize_nodes(
        newer_nodes,
        new_node_dict,
        new_type_dict
      );
      var new_node_dict = new_node_elements[0];
      var new_type_dict = new_node_elements[1];
      var new_expression_dict = new_node_elements[2];
      var new_stats_dict = new_node_elements[3];
      var new_display_analytes_dict = new_node_elements[4];
      var new_display_reactions_dict = new_node_elements[5];
      var new_entity_id_dict = new_node_elements[6];

      make_graph(
        data,
        newer_nodes,
        newer_links,
        new_type_dict,
        new_node_dict,
        new_entity_id_dict,
        new_expression_dict,
        new_stats_dict,
        new_display_analytes_dict,
        new_display_reactions_dict
      );
    }
  });

  d3.select("#collapseNodes").on("click", function() {
    if (collapse_reactions === false) {
      collapse_reactions = true;
      new_nodes = collapsed_nodes;
      new_links = collapsed_links;

      make_graph(
        data,
        new_nodes,
        new_links,
        type_dict,
        node_dict,
        entity_id_dict,
        expression_dict,
        stats_dict,
        display_analytes_dict,
        display_reactions_dict
      );
    } else {
      collapse_reactions = false;
      collapsed_nodes = new_nodes;
      collapsed_links = new_links;

      var selection = document.getElementById("pathwayMenu").value;
      var reactions = collapsed_pathway_dict[selection]["reactions"];
      var collapsed_reaction_dictionary = data.collapsed_reaction_dictionary;

      // Parse through each reaction listed and get the component parts
      var components = [];
      for (rxn in reactions) {
        var target_rxns = collapsed_reaction_dictionary[reactions[rxn]];
        components.push(reactions[rxn]);
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
      var newer_elements = get_nodes_links(data, components);
      var newer_nodes = newer_elements[0];
      var newer_links = newer_elements[1];

      // Initialize variables
      var new_node_dict = {};
      var new_type_dict = {};

      var new_node_elements = initialize_nodes(
        newer_nodes,
        new_node_dict,
        new_type_dict
      );
      var new_node_dict = new_node_elements[0];
      var new_type_dict = new_node_elements[1];
      var new_expression_dict = new_node_elements[2];
      var new_stats_dict = new_node_elements[3];
      var new_display_analytes_dict = new_node_elements[4];
      var new_display_reactions_dict = new_node_elements[5];
      var new_entity_id_dict = new_node_elements[6];

      make_graph(
        data,
        newer_nodes,
        newer_links,
        new_type_dict,
        new_node_dict,
        new_entity_id_dict,
        new_expression_dict,
        new_stats_dict,
        new_display_analytes_dict,
        new_display_reactions_dict
      );
    }
  });


  var cell = node.append("path").attr("class", "cell");

  // Draw curved edges
  function tick() {
    link.attr("d", linkArc);
    circle.attr("transform", transform);
    text.attr("transform", transform);
  }

  function dragsubject() {
    return simulation.find(d3.event.x, d3.event.y);
  }

  function dragstarted() {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    d3.event.subject.fx = d3.event.subject.x;
    d3.event.subject.fy = d3.event.subject.y;
  }

  function dragged() {
    d3.event.subject.fx = d3.event.x;
    d3.event.subject.fy = d3.event.y;
  }

  function dragended() {
    if (!d3.event.active)
      simulation
        .alphaTarget(0.05)
        .alphaMin(0.06)
        .velocityDecay(0.7);
    d3.event.subject.fx = null;
    d3.event.subject.fy = null;
  }

  function drawLink(d) {
    context.moveTo(d.source.x, d.source.y);
    context.lineTo(d.target.x, d.target.y);
  }

  function drawNode(d) {
    context.moveTo(d.x + 3, d.y);
    context.arc(d.x, d.y, 3, 0, 2 * Math.PI);
  }
}

function changeSuper() {
  var superSelection = document.getElementById("superPathwayMenu").value;
  emptyMenu(document.getElementById("pathwayMenu"));

  // Limit pathways to super-pathway
  if (superSelection === "All pathways") {
    var parsed_pathway_dict = pathway_dict;
  } else if (superSelection === "All entities") {
    var entity_dictionary = parseEntities(data.nodes);
  } else {
    var selectedReactions = superPathwayDict[superSelection]["reactions"];
    var parsed_pathway_dict = parsePathways(pathway_dict, selectedReactions);
  }

  if (superSelection !== "All entities") {
    make_menu(parsed_pathway_dict, "pathwayMenu", "Select a pathway...");
  } else {
    make_menu(entity_dictionary, "pathwayMenu", "Select an entity...");
  }
}

// Create dictionary of all entity names and their corresponding entity IDs
function parseEntities(nodes) {
  entity_dictionary = {};
  for (node in nodes) {
    if (nodes[node].type !== "reaction") {
      entity_dictionary[nodes[node].name] = nodes[node].id;
    }
  }

  return entity_dictionary;
}

// Graphing
function change() {
  graph_genes = true;
  collapse_reactions = true;
  var selection = document.getElementById("pathwayMenu").value;
  var superSelection = document.getElementById("superPathwayMenu").value;

  document.getElementById("reaction_notes").innerHTML = "";

  var mod_selection = determineWidth(selection);
  document.getElementById("type_selection_type").innerHTML = "Pathway";
  document.getElementById("type_selection").innerHTML = mod_selection;
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

  if (superSelection !== "All entities") {
    // Run normal first plot
    var reactions = pathway_dict[selection]["reactions"];
    var elements = parse_pathway(data, reactions);
    var new_nodes = elements[0];
    var new_links = elements[1];

    // Initialize variables
    var node_dict = {};
    var type_dict = {};

    var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
    var node_dict = node_elements[0];
    var type_dict = node_elements[1];
    var expression_dict = node_elements[2];
    var stats_dict = node_elements[3];
    var display_analytes_dict = node_elements[4];
    var display_reactions_dict = node_elements[5];
    var entity_id_dict = node_elements[6];

    make_graph(
      data,
      new_nodes,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      expression_dict,
      stats_dict,
      display_analytes_dict,
      display_reactions_dict
    );
  } else {
    // Get kNN of entity selected
    var mod_selection = determineWidth(selection);
    document.getElementById("type_selection").innerHTML = mod_selection;
    document.getElementById("type_selection_type").innerHTML =
      "Nearest Neighbor";

    graph_genes = true;
    var entity_dictionary = parseEntities(data.nodes);
    nearest_neighbors(data, entity_dictionary[selection]);
  }
}

// Check number of categories
function checkCategories(categories) {
  //change to > 1 after testing
  if (data.categories.length > 1) {
    timecourse = true;
    timecourse_fill =
      '<div id="play-button" align="center">Pause<div id="bar" align="center"></div></div>';
    document.getElementById("slider").innerHTML = timecourse_fill;

    //buildSlider();
  } else {
    timecourse = false;
  }
  return timecourse;
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

var timecourse = checkCategories(data.categories);

make_menu(
  superPathwayDict,
  "superPathwayMenu",
  "Select a category...",
  (provide_all = true)
);

d3.select("#superPathwayMenu").on("change", changeSuper);
d3.select("#pathwayMenu").on("change", change);
