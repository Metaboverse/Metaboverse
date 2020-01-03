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
const d3 = require('d3')

const max_nodes = 1500;

// Change user selection based on input
var selection = null;
function selectPathway() {

  var selection = document.getElementById("pathwayMenu").value;
  console.log("Pathway: " + selection)
  return selection;
};

// Populate dictionary to access component reactions for each pathway
function make_pathway_dictionary(data) {

  // Get pathway name and ID
  var master = data.master_reference;
  var pathways = data.pathway_dictionary;
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



  var master = data.master_reference;
  var reactions_dictionary = data.reactions_dictionary;

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

  var nodes = data["nodes"];
  var links = data["links"];

  // Parse the nodes of interest
  var new_nodes = [];
  var node_components = [];
  nodes.forEach(function(node) {

    if (components.includes(node['entity_id']) || components.includes(node['name'])) {

      var node_copy = $.extend(true,{},node)

      node_copy["rgba"] = node_copy["rgba"][0]
      node_copy["rgba_js"] = node_copy["rgba_js"][0]

      new_nodes.push(node_copy)
      node_components.push(node_copy["name"])
      node_components.push(node_copy["entity_id"])

    };

  });

  var new_links = []
  // Parse out links of interest
  links.forEach( function (link) {

    if (node_components.includes(link.source) && node_components.includes(link.target)) {
      var link_copy = $.extend(true,{},link)
      new_links.push(link_copy)
    };

  });

  return [new_nodes, new_links];

};

function nearest_neighbors(data, entity_id) {

  var master = data.master_reference;

  // Get current nearest neighbors value
  var kNN = document.getElementById("kNN_button").value;

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
  var display_analytes_dict = node_elements[3];
  var display_reactions_dict = node_elements[4];
  var entity_id_dict = node_elements[5];

  make_graph(
      data,
      new_nodes,
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

  // Reset warning messages
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

  var master = data.master_reference;
  var reactions_dictionary = data.reactions_dictionary;

  // Parse through each reaction where entity is a component
  nn_components = [entity_id];
  for (reaction in reactions_dictionary) {

    //Return all reactions that contain the entity
    if (reactions_dictionary[reaction].includes(entity_id)) {

      for (entity in reactions_dictionary[reaction]) {

        nn_components.push(reactions_dictionary[reaction][entity]);

      };
    };
  };

  var escape = nn_components.length
  console.log(escape)
  // If too many nodes for first neighborhood, stop plotting
  if (escape > max_nodes) {

    document.getElementById("warning_line_1").innerHTML = "<i style=\"color: red;\">Too many entities to plot</i>";
    document.getElementById("warning_line_2").innerHTML = "<i style=\"color: red;\">Will not plot<br></i>";
    kNN = 0
    nn_components = [];
  }

  if (kNN > 1) {

    document.getElementById("warning_line_1").innerHTML = "<i style=\"color: red;\">Please wait<br><br></i>";
    document.getElementById("warning_line_2").innerHTML = "";

    n = 1
    while (n < kNN) {

      var components = [];

      for (element in nn_components) {

        for (reaction in reactions_dictionary) {

          if (reactions_dictionary[reaction].includes(nn_components[element])) {

            for (el in reactions_dictionary[reaction]) {

              components.push(reactions_dictionary[reaction][el]);

            };
          };
        };
      }
      n = n + 1

      // Only combine the lists if they pass the threshold
      escape = escape + components.length
      console.log(escape)
      if (escape > max_nodes) {
        n = kNN + 2
        document.getElementById("warning_line_1").innerHTML = "<i style=\"color: red;\">Too many entities to plot. Will only plot first neighborhood</i>";
        document.getElementById("warning_line_2").innerHTML = "";

      } else {
        nn_components = nn_components.concat(components)
        document.getElementById("warning_line_1").innerHTML = "<br>";
        document.getElementById("warning_line_2").innerHTML = "<br><br>";
      }

    }
  }
  console.log("++++++")
  var new_elements = get_nodes_links(data, nn_components);

  return new_elements;

};

function kNN_input(d) {

  var knn_value = document.getElementById("kNN_button").value;
  console.log("k-NN parameter now set to: " + knn_value)
};

function make_graph(
    data,
    new_nodes,
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
    .select("#graph")
    .append("svg")
      .attr("width", width - 5)
      .attr("height", height - 80)
    .call(d3.zoom().on("zoom", function () {
      svg.attr("transform", d3.event.transform)
      }))
    .append("g")

  const forceX = d3.forceX(width / 2).strength(0.015)
  const forceY = d3.forceY(height / 2).strength(0.015)

  const simulation = d3.forceSimulation(new_nodes)
      .force("link", d3.forceLink(new_links).id(d => d.name)
        .distance(40)
        .strength(1))
      .force("charge", d3.forceManyBody().strength(-1000))
      .force("center", d3.forceCenter(width / 2, height / 2))

  simulation.alphaTarget(0.01).alphaMin(0.1).velocityDecay(0.70)

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

  var link = svg.append("g").selectAll("path")
    .data(new_links)
    .enter().append("path")
      .attr("class", function(d) { return "link " + d.type; })
      .attr("marker-end", function(d) { return "url(#" + d.type + ")"; });

  var node = svg.selectAll(".node")
    .data(new_nodes)
    .enter().append("g")
      .attr("class", "node")
    .style("--node_color", function(d) {
      return "rgba(" + d.rgba_js.toString() + ")";
      }
    )
    .style("--node_radius", function(d) { return 6; })
    .call(d3.drag()
          .subject(dragsubject)
          .on("start", dragstarted)
          .on("drag", dragged)
          .on("end", dragended));

  var circle = node
    .append("circle")
      .attr("r", 6)
    .on("dblclick", function(d) {

      if ((type_dict[entity_id_dict[d["name"]]] === "reaction") || (type_dict[d["name"]] === "reaction")) {

        console.log("Selected a reaction, will not perform kNN graphing")

      } else {
        document.getElementById("type_selection_type").innerHTML = "Nearest Neighbor";
        document.getElementById("type_selection").innerHTML = d["name"];
        nearest_neighbors(data, entity_id_dict[d["name"]]);

      };

    });

  var text = node
    .append("text")
      .attr("x", 14)
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

  simulation.on("tick", tick);

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

  // Draw curved edges
  function tick() {

    link.attr("d", linkArc);
    circle.attr("transform", transform);
    text.attr("transform", transform);

  };

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
    if (!d3.event.active) simulation.alphaTarget(0.05).alphaMin(0.06).velocityDecay(0.70);
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

};

// MAIN
database_url = get_session_info("database_url")
console.log("Database path: " + database_url)

var data = JSON.parse(fs.readFileSync(database_url).toString());

var pathway_dict = make_pathway_dictionary(data);
make_menu(pathway_dict);
d3.select("#pathwayMenu").on("change", change);

// Graphing
function change() {

  var selection = document.getElementById("pathwayMenu").value;
  document.getElementById("type_selection_type").innerHTML = "Pathway";
  document.getElementById("type_selection").innerHTML = selection;
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

  var reactions = pathway_dict[selection]['reactions'];
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
  var display_analytes_dict = node_elements[3];
  var display_reactions_dict = node_elements[4];
  var entity_id_dict = node_elements[5];

  make_graph(
      data,
      new_nodes,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      expression_dict,
      display_analytes_dict,
      display_reactions_dict)

  };
