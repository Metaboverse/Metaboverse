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


Portions of the force graphing below based on or adapted from code from Mike Bostock
The original code is under the GNU General Public License v3.0, allowing for modification
and distribution
License and copyright notice: GNU General Public License v3.0
Changes:
  - Heavily modified and added to the style CSS for more flexibility in plotting
  - Adapted general D3 plotting functions and commands to work with input data and accept flexibility
  - Modified plotting functions to allow for the differential shading of nodes
  - All other components are original
Source:
http://bl.ocks.org/mbostock/1153292
https://bl.ocks.org/mbostock/1212215
*/

// Change user selection based on input
var selection = null;
console.log(selection)
function selectPathway() {

  var selection = document.getElementById("pathwayMenu").value;
  console.log(selection)
  return selection;
};

// Populate dictionary to access component reactions for each pathway
function make_pathway_dictionary(data) {

  // Get pathway name and ID
  var master = data[0].master_reference;
  var pathways = data[0].pathway_dictionary;
  var pathway_dict = {}
  for (var key in pathways) {

    if (key in master) {
      var name = master[key];
      pathway_dict[name] = {
        'id': key,
        'reactions': pathways[key]
      };

    } else {
      //var name = key; //eventually need to figure out mapping for remaining pathways that dont have a name, maybe they should be R-ALL
    };

  };

  return pathway_dict;
};

// Make Pathway menu for users to
function make_menu(pathway_dict) {

  // Get species names (keys) as list
  pathways_list = Object.getOwnPropertyNames(
    pathway_dict
  ).map(function(k) {
    return k;
  });
  pathways_list.sort();
  pathways_list.unshift("Select a pathway..."); // Add select prompt to menu bar

  // Generate drop-down menu for species select
  var menu = document.getElementById("pathwayMenu");
  for (var i = 0; i < pathways_list.length; i++) {

    var option = document.createElement("option");
    option.innerHTML = pathways_list[i];
    option.value = pathways_list[i];
    menu.appendChild(option);

  };

};

function initialize_nodes(nodes, node_dict, type_dict) {

  var expression_dict = {};
  var display_analytes_dict = {};
  var display_reactions_dict = {};
  var entity_id_dict = {};

  // Make dictionary of node color values and types
  nodes.forEach(function(node) {

    node_dict[node['name']] = node['rgba_js']
    type_dict[node['name']] = node['type']
    expression_dict[node['name']] = node['expression'][0]
    entity_id_dict[node['name']] = node['entity_id']
    entity_id_dict[node['entity_id']] = node['name']

    if (node['type'] === 'reaction') {
      display_analytes_dict[node['name']] = 'none'
      display_reactions_dict[node['name']] = 'inline'
    } else {
      display_reactions_dict[node['name']] = 'none'
      display_analytes_dict[node['name']] = 'inline'
    };

  });

  return [node_dict, type_dict, expression_dict, display_analytes_dict, display_reactions_dict, entity_id_dict];

};

function initialize_links(links, nodex, node_dict) {

  console.log(">>>>>>>>>")
  console.log(links)
  // Populate node object with relevant information and data
  links.forEach(function(link) {

    // Add pertinent node information
    console.log(link)
    console.log(nodex)
    console.log("+++++")

    //link.source = nodex[link.source] || (nodex[link.source] = {name: link.source});
    //link.target = nodex[link.target] || (nodex[link.target] = {name: link.target});

    link.source = nodex[link.source] || (nodex[link.source] = {name: link.source});
    link.target = nodex[link.target] || (nodex[link.target] = {name: link.target});

    console.log(link)
    console.log(nodex)
    console.log("=====")

    nodex[link.source["name"]]["color"] = node_dict[link.source["name"]];
    nodex[link.target["name"]]["color"] = node_dict[link.target["name"]];

  });

  return nodex;

};

function linkArc(d) {

  var dx = d.target.x - d.source.x;
  var dy = d.target.y - d.source.y;
  var dr = Math.sqrt(dx * dx + dy * dy);

  return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;

};

function transform(d) {

  return "translate(" + d.x + "," + d.y + ")";

};

function parse_pathway(data, reactions) {

  var master = data[0].master_reference;
  var reactions_dictionary = data[0].reactions_dictionary;

  // Parse through each reaction listed and get the component parts
  var components = [];
  for (rxn in reactions) {

    var target_rxns = reactions_dictionary[reactions[rxn]];
    for (x in target_rxns) {

      components.push(target_rxns[x]);
      components.push(master[target_rxns[x]]);

    };

  };

  var new_elements = get_nodes_links(data, components);

  return new_elements;

};

function get_nodes_links(data, components) {

  var nodes = data[0].nodes;
  var links = data[0].links;

  // Parse the nodes of interest
  var new_nodes = [];
  var node_components = [];
  nodes.forEach(function(node) {

    if (components.includes(node['entity_id']) || components.includes(node['name'])) {

      new_nodes.push(node)
      node_components.push(node["name"])
      node_components.push(node["entity_id"])

    };

  });

  // Parse out links of interest
  var new_links = [];
  links.forEach(function(link) {

    if ((node_components.includes(link.source) && node_components.includes(link.target)) ||
        (node_components.includes(link.source["name"]) && node_components.includes(link.target["name"]))) {
      console.log("relinking")
      new_links.push(link);

    };

  });

  console.log(new_links)

  return [new_nodes, new_links];

};

function nearest_neighbors(data, entity_id) {

  var master = data[0].master_reference;

  // Get current nearest neighbors value
  var kNN = document.getElementById("kNN_button").value;

  // Curate kNN to node of interest
  var nn_elements = parse_kNN_pathway(data, entity_id, kNN);
  var new_nodes = nn_elements[0];
  var new_links = nn_elements[1];

  // Save old graph in burner variable

  // Recurate the graph
  // Initialize variables
  var nodex = {};
  var node_dict = {};
  var type_dict = {};

  var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
  var node_dict = node_elements[0];
  var type_dict = node_elements[1];
  var expression_dict = node_elements[2];
  var display_analytes_dict = node_elements[3];
  var display_reactions_dict = node_elements[4];
  var entity_id_dict = node_elements[5];

  console.log("?????????????????????????????")
  console.log(new_links)
  console.log("?????????????????????????????")

  var nodex = initialize_links(new_links, nodex, node_dict);

  make_graph(
      data,
      nodex,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      expression_dict,
      display_analytes_dict,
      display_reactions_dict)



  // Remove old plot and plot this one

  // Provide a go back button to get back to the graph they made previously
  // Run by dblclick on node of interest
  // Make node of interest bigger than others

};

function parse_kNN_pathway(data, entity_id, kNN) {

  var master = data[0].master_reference;
  var reactions_dictionary = data[0].reactions_dictionary;

  // Parse through each reaction where entity is a component
  var components = [];

  //Return all reactions that contain the entity
  nn_components = [entity_id];
  for (reaction in reactions_dictionary) {

    if (reactions_dictionary[reaction].includes(entity_id)) {

      for (entity in reactions_dictionary[reaction]) {

        nn_components.push(reactions_dictionary[reaction][entity]);

      };

    };

  };

  var new_elements = get_nodes_links(data, nn_components);

  return new_elements;

};

function kNN_input(d) {

  var knn_value = document.getElementById("kNN_button").value;
  console.log(knn_value)
};





function make_graph(
    data,
    nodex,
    new_links,
    type_dict,
    node_dict,
    entity_id_dict,
    expression_dict,
    display_analytes_dict,
    display_reactions_dict) {

  // Allow flexible window dimensions based on initial window size when opened
  var width = window.innerWidth;
  var height = window.innerHeight;

  // Restart graph
  d3.select("svg").remove();
  d3.select("force").remove();
  d3.select("g_nodes").remove();

  // Initialize force graph object
  var svg = d3
    .select("body")
    .append("svg")
      .attr("width", width)
      .attr("height", height)
    .call(d3.behavior.zoom()
    .on("zoom", activate_zoom));

  var force = d3.layout.force()
    .gravity(0.1)
    .linkDistance(60)
    .charge(-500)
    .on("tick", tick)
    .size([width, height]);

  var g_nodes = d3.values(nodex);

  // Build graph
  force
    .nodes(g_nodes)
    .links(new_links)
    .start();

  // Generate edges with style attributes
  svg.append("defs").selectAll("marker")
    .data([
      "reaction",
      "reactant",
      "product",
      "inhibitor",
      "catalyst",
      "complex_component"])
    .enter()
    .append("marker")
      .attr("id", function(d) { return d; })
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 15)
      .attr("refY", -1.5)
      .attr("markerWidth", 6)
      .attr("markerHeight", 6)
      .attr("orient", "auto")
    .append("path")
      .attr("d", "M0, -5L10, 0L0, 5");

  var path = svg.append("g").selectAll("path")
    .data(force.links())
    .enter().append("path")
      .attr("class", function(d) { return "link " + d.type; })
      .attr("marker-end", function(d) { return "url(#" + d.type + ")"; });

  var node = svg.selectAll(".node")
    .data(g_nodes)
    .enter().append("g")
      .attr("class", "node")
    .style("--node_color", function(d) {
      try {
        return "rgba(" + d.color[0].toString() + ")";
      } catch (e) {
        return "rgba(" + d.name.color[0].toString() + ")";
      }})
    .style("--node_radius", function(d) { return 6; })
    .call(force.drag);

  var circle = node
    .append("circle")
      .attr("r", 6)
    .on("dblclick", function(d) {

      if ((type_dict[entity_id_dict[d["name"]]] === "reaction") || (type_dict[d["name"]] === "reaction")) {

        console.log("Selected a reaction, will not perform kNN graphing")

      } else {

        console.log(entity_id_dict[d["name"]])
        nearest_neighbors(data, entity_id_dict[d["name"]]);


      };

    });

  var text = node
    .append("text")
      .attr("x", 16)
      .attr("y", ".31em")
    .text(function(d) {

      if (type_dict[d.name] === "reaction") {

        // If reaction node, do not display expression value
        return d.name;

      } else {

        // Label other nodes with expression value in parentheses
        return d.name + ' (' + parseFloat(expression_dict[d.name]).toFixed(2) + ')';

      }
    });

  // Not working right now
  toggle_e = true;
  d3.select("#toggleExpression")
    .on("click", function() {

      if (toggle_e === false) {

        toggle_e = true;
        text.text(function(d) {

          if (type_dict[d.name] === "reaction") {

            // If reaction node, do not display expression value
            return d.name;

          } else {

            // Label other nodes with expression value in parentheses
            return d.name + ' (' + parseFloat(expression_dict[d.name]).toFixed(2) + ')';

          }

        });

      } else {
        toggle_e = false;
        text.text(function(d) { return d.name });
      }

   });

  toggle_a = false;
  toggle_r = false;
  text.style("--node_display", function(d) { return "none"; });

  d3.select("#toggleAnalytes")
    .on("click", function() {

      if (toggle_a === false) {

        toggle_a = true;
        determine_displays(toggle_a, toggle_r);

      } else {

        toggle_a = false;
        determine_displays(toggle_a, toggle_r);
      }

    });

  d3.select("#toggleReactions")
    .on("click", function() {

      if (toggle_r === false) {

        toggle_r = true;
        determine_displays(toggle_a, toggle_r);

      } else {

        toggle_r = false;
        determine_displays(toggle_a, toggle_r);

      }

    });

  function determine_displays(toggle_a, toggle_r) {

    if ((toggle_a === false) && (toggle_r === false)) {
      var display_labels = "none";
    }
    else if ((toggle_a === false) && (toggle_r === true)) {
      var display_labels = "reactions";
    }
    else if ((toggle_a === true) && (toggle_r === false)) {
      var display_labels = "analytes";
    }
    else if ((toggle_a === true) && (toggle_r === true)) {
      var display_labels = "all";
    }
    else {
      var display_labels = "none";
    }

    // options:

    // none -> all labels hidden until hovered
    if (display_labels === "none") {
      text.style("--node_display", function(d) { return "none"; })
    }
    // analytes -> analytes shown, reactions hovered
    else if (display_labels === "analytes") {
      text.style("--node_display", function(d) { return display_analytes_dict[d.name]; })
    }
    // reactions -> reactions shown, analytes hovered
    else if (display_labels === "reactions") {
      text.style("--node_display", function(d) { return display_reactions_dict[d.name]; })
    }
    // analytes + reactions -> all -> all shown
    else if (display_labels === "all") {
      text.style("--node_display", function(d) { return "inline"; })
    }
    else {
      text.style("--node_display", function(d) { return "none"; })
    }

  };

  var cell = node
    .append("path")
      .attr("class", "cell");

  d3.select("body")
    .on("keydown", function () {

      toggle_zoom = d3.event.altKey;

    });

  d3.select("body")
    .on("keyup", function () {

      toggle_zoom = false;

    });

  // Draw curved edges
  function tick() {

    path.attr("d", linkArc);
    circle.attr("transform", transform);
    text.attr("transform", transform);

  };

  // Toggle zoom and pan by pressing the Alt key
  // This section adapted from Pedro Tabacof, https://stackoverflow.com/a/34815469/9571488
  var toggle_zoom = false;
  function activate_zoom() {

      if (toggle_zoom === true) {

        svg.attr(
          "transform",
          "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")"

        );
      }
  };
};

// MAIN
d3.json("data/HSA_global_reactions.json", function(data) {

  var pathway_dict = make_pathway_dictionary(data);
  make_menu(pathway_dict);
  d3.select("#pathwayMenu").on("change", change);

  // Graphing
  function change() {

    var selection = document.getElementById("pathwayMenu").value;

    var reactions = pathway_dict[selection]['reactions'];
    var elements = parse_pathway(data, reactions);
    var new_nodes = elements[0];
    var new_links = elements[1];

    // Initialize variables
    var nodex = {};
    var node_dict = {};
    var type_dict = {};

    var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
    var node_dict = node_elements[0];
    var type_dict = node_elements[1];
    var expression_dict = node_elements[2];
    var display_analytes_dict = node_elements[3];
    var display_reactions_dict = node_elements[4];
    var entity_id_dict = node_elements[5];

    console.log("?????????????????????????????")
    console.log(new_links)
    console.log("?????????????????????????????")

    var nodex = initialize_links(new_links, nodex, node_dict);

    make_graph(
        data,
        nodex,
        new_links,
        type_dict,
        node_dict,
        entity_id_dict,
        expression_dict,
        display_analytes_dict,
        display_reactions_dict)

    };

  });
