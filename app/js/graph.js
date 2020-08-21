/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2020 Jordan A. Berg
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

var d3 = require("d3");
var fs = require("fs");
var saveSVG = require("save-svg-as-png");

const hullPadding = 60;
const max_nodes = 1500;
var sample = 0;
var entity = "values_js";
var knn_value = 1;
var hub_value = 1000000;
var stat_value = 0.05;
var graph_genes = true;
var collapse_reactions = true;
var saved_nodes = [];
var saved_links = [];
var collapsed_nodes = [];
var collapsed_links = [];

// Define the div for the tooltips
var div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

d3.select("button#options_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "56px");
      div
        .html("Click on the following buttons to toggle the display of the respective feature. Hover over the button for more details.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#knn_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "97px");
      div
        .html("Change value to expand number of nearest neighbors to plot for the selected analyte node. To expand neighbors, double-click on a non-reaction node. (Currently limited to 2 neighborhoods due to rapid graph expansion).")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#hub_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "97px");
      div
        .html("Provide a value to threshold displayed entities to those that have no more than the number of connections to other entities you provide. By default, graphing will display all entities no matter the number of connections.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#stat_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "40px");
      div
        .html("Provide a value to threshold statistical value where node borders are bolded.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#notes_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "54px");
      div
        .html("Display reaction details by double-clicking on a reaction node. Display metabolite synonyms by single-clicking on metabolite node.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#notes_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "54px");
      div
        .html("Display reaction details by double-clicking on a reaction node. Display metabolite synonyms by single-clicking on metabolite node.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);
d3.select("button#legend_info")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "43px");
      div
        .html("Click and drag the background to pan, or use the mouse wheel to zoom.")
      }
    )
  .on("mouseout", function(d) {
  div.style("opacity", 0);
  div.html("")
  }
);

d3.select("button#shape_legend")
  .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "435px")
        .style("width", "200px");
      div
        .html("<div style='margin-left:15px;margin-top:15px;'><font size='3'><b><u>Relationships</u></b></font></br><div class='arrow'><div class='line grey-arrow'></div><div class='point grey-arrow'></div></div>&nbsp;&nbsp;&nbsp;Core interaction</br><div class='arrow'><div class='line green-arrow'></div><div class='point green-arrow'></div></div>&nbsp;&nbsp;&nbsp;Catalyst</br><div class='arrow'><div class='line red-arrow'></div><div class='point red-arrow'></div></div>&nbsp;&nbsp;&nbsp;Inhibitor</br><div class='arrow'><div class='line2 blue-arrow'></div></div>&nbsp;&nbsp;&nbsp;Metabolite Component</br><div class='arrow'><div class='line2 orange-arrow'></div></div>&nbsp;&nbsp;&nbsp;Protein Component</br><div class='arrow'><div class='line2 purple-arrow'></div></div>&nbsp;&nbsp;&nbsp;Gene Component</br></br><font size='3'><b><u>Entities</u></b></font></br></br><span class='fas fa-star fa-lg grey-shader'></span>&nbsp;&nbsp;&nbsp;&nbsp;Reaction</br></br><span class='fas fa-star fa-lg purple-shader'></span>&nbsp;&nbsp;&nbsp;&nbsp;Regulated reaction</br></br>&nbsp;<span class='dot white-dot'></span>&nbsp;&nbsp;&nbsp;&nbsp;Metabolite</br></br>&nbsp;<span class='square white-dot'></span>&nbsp;&nbsp;&nbsp;&nbsp;Complex</br></br>&nbsp;<span class='diamond white-dot'></span>&nbsp;&nbsp;&nbsp;&nbsp;Protein</br></br>&nbsp;<span class='black-triangle'><span class='white-triangle'></span></span>&nbsp;&nbsp;&nbsp;&nbsp;Gene</br></br>* Bolded shapes meet statistical threshold</div>")
      }
    )
  .on("mouseout", function(d) {
    div.style("opacity", 0);
    div.html("")
  }
);

function checkReaction(
    reaction,
    element) {
  // Check that element is included in reaction components

  if (reaction["reactants"].includes(element)
      || reaction["products"].includes(element)
      || reaction["modifiers"].includes(element)
      || reaction["additional_components"].includes(element)) {
    return true;
  } else {
    return false;
  }
}

function checkPlotting(
    data,
    reaction) {
  // Check whether a neighboring reaction should be included

  let these_degrees = [];
  let rxn_comps = [];
  for (x in reaction["reactants"]) {
    these_degrees.push(data.nodes[reaction["reactants"][x]].degree);
    rxn_comps.push(data.nodes[reaction["reactants"][x]].name);
  }
  for (x in reaction["products"]) {
    these_degrees.push(data.nodes[reaction["products"][x]].degree);
    rxn_comps.push(data.nodes[reaction["products"][x]].name);
  }
  for (x in reaction["modifiers"]) {
    these_degrees.push(data.nodes[reaction["modifiers"][x][0]].degree);
    rxn_comps.push(data.nodes[reaction["modifiers"][x][0]].name);
  }
  for (x in reaction["additional_components"]) {
    these_degrees.push(data.nodes[reaction["additional_components"][x]].degree);
    rxn_comps.push(data.nodes[reaction["additional_components"][x]].name);
  }
  let _max = Math.max(...these_degrees);

  let blocklisted_node = false;
  for (r in rxn_comps) {
    if (data.blocklist.includes(rxn_comps[r])) {
      blocklisted_node = true;
    }
  }

  return [_max, blocklisted_node];
}

function initialize_nodes(
    nodes,
    node_dict,
    type_dict) {
  var display_analytes_dict = {};
  var display_reactions_dict = {};
  var entity_id_dict = {};

  // Make dictionary of node color values and types
  nodes.forEach(function(node) {
    node_dict[node["name"]] = node["values_js"];
    type_dict[node["name"]] = node["type"];

    entity_id_dict[node["name"]] = node["id"];
    entity_id_dict[node["id"]] = node["name"];

    if ((node["type"] === "reaction") || (node["type"] === "collapsed")) {
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

function parse_pathway(
    data,
    reactions,
    reaction_dictionary,
    degree_dictionary) {

  // Parse through each reaction listed and get the component parts
  var components = [];
  for (rxn in reactions) {
    var target_rxns = reaction_dictionary[reactions[rxn]];

    if (target_rxns !== undefined) {
      components.push(reactions[rxn]);
      for (x in target_rxns["reactants"]) {
        if (degree_dictionary[target_rxns["reactants"][x]] <= hub_value) {
          components.push(target_rxns["reactants"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["reactants"][x] +
              " (" +
              degree_dictionary[target_rxns["reactants"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in target_rxns["products"]) {
        if (degree_dictionary[target_rxns["products"][x]] <= hub_value) {
          components.push(target_rxns["products"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["products"][x] +
              " (" +
              degree_dictionary[target_rxns["products"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in target_rxns["modifiers"]) {
        if (degree_dictionary[target_rxns["modifiers"][x][0]] <= hub_value) {
          components.push(target_rxns["modifiers"][x][0]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["modifiers"][x][0] +
              " (" +
              degree_dictionary[target_rxns["modifiers"][x][0]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in target_rxns["additional_components"]) {
        if (degree_dictionary[target_rxns["additional_components"][x]] <= hub_value) {
          components.push(target_rxns["additional_components"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["additional_components"][x] +
              " (" +
              degree_dictionary[target_rxns["additional_components"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
    } else {
      console.log("Could not access the following:")
      console.log(rxn)
      console.log(reactions[rxn])
      console.log(reaction_dictionary[reactions[rxn]])
      alert("Could not access the following id: ", reactions[rxn], ". This could be related to a failure to updated de-duplicated reactions.")
    }

  }
  var new_elements = get_nodes_links(data, components);

  return new_elements;
}

function get_nodes_links(data, components) {

  var nodes = data.nodes;
  var links = data.links;
  components = [...new Set(components)];

  // Parse the nodes of interest
  var add_nodes = [];
  for (c in components) {
    if (components[c] in nodes) {
      var node_copy = $.extend(true, {}, nodes[components[c]]);
      add_nodes.push(node_copy);
    }
  }

  var add_links = [];
  for (c1 in components) {
    for (c2 in components) {
      _key = components[c1] + "," + components[c2];
      if (_key in links && components[c1] !== components[c2]) {
        var link_copy = $.extend(true, {}, links[_key]);
        add_links.push(link_copy);
      }
    }
  }

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
  var display_analytes_dict = node_elements[2];
  var display_reactions_dict = node_elements[3];
  var entity_id_dict = node_elements[4];

  // Remove old plot and plot this one
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
    _height,
    global_motifs)
}

function parse_kNN_pathway(data, entity_id, kNN) {
  // Reset warning messages
  try {
    document.getElementById("warning_line_1").innerHTML = "<br>";
    document.getElementById("warning_line_2").innerHTML = "<br><br>";
  } catch(e) {}

  if (collapse_reactions !== true) {
    var reaction_dictionary = data.reaction_dictionary;
  } else {
    var reaction_dictionary = data.collapsed_reaction_dictionary;
  }
  var degree_dictionary = data.degree_dictionary;

  // Parse through each reaction where entity is a component
  nn_components = [entity_id];
  for (reaction in reaction_dictionary) {
    //Return all reactions that contain the entity
    let target_rxns = reaction_dictionary[reaction];
    if (checkReaction(
        target_rxns,
        entity_id) === true) {
      nn_components.push(target_rxns["id"]);
      for (x in target_rxns["reactants"]) {
        if (degree_dictionary[target_rxns["reactants"][x]] <= hub_value) {
          nn_components.push(target_rxns["reactants"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["reactants"][x] +
              " (" +
              degree_dictionary[target_rxns["reactants"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in reaction_dictionary[reaction]["products"]) {
        if (degree_dictionary[target_rxns["products"][x]] <= hub_value) {
          nn_components.push(target_rxns["products"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["products"][x] +
              " (" +
              degree_dictionary[target_rxns["products"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in reaction_dictionary[reaction]["modifiers"]) {
        if (degree_dictionary[target_rxns["modifiers"][x][0]] <= hub_value) {
          nn_components.push(target_rxns["modifiers"][x][0]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["modifiers"][x][0] +
              " (" +
              degree_dictionary[target_rxns["modifiers"][x][0]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
      for (x in reaction_dictionary[reaction]["additional_components"]) {
        if (degree_dictionary[target_rxns["additional_components"][x]] <= hub_value) {
          nn_components.push(target_rxns["additional_components"][x]);
        } else {
          console.log(
            "Filtering " +
              target_rxns["additional_components"][x] +
              " (" +
              degree_dictionary[target_rxns["additional_components"][x]] +
              " degrees) -- may cause edge loss"
          );
        }
      }
    }
  }

  // If too many nodes for first neighborhood, stop plotting
  var escape = nn_components.length;
  if (escape > max_nodes) {
    try {
      document.getElementById("warning_line_1").innerHTML =
        '<i class="red-text">Too many entities to plot</i><br><i class="red-text">Will not plot</i>';
      document.getElementById("warning_line_2").innerHTML = "<br>";
      alert("Too many entities to plot")
    } catch(e) {}

    kNN = 0;
    nn_components = [];
  }

  if (kNN > 1) {
    try {
      document.getElementById("warning_line_1").innerHTML =
        '<i class="red-text">Please wait<br><br></i>';
      document.getElementById("warning_line_2").innerHTML = "";
    } catch(e) {}
    // Filter out any components that are above hub threshold for kNN
    n = 1;
    nn_components = [...new Set(nn_components)];
    while (n < kNN) {
      var components = [];

      for (element in nn_components) {
        if (data.nodes[nn_components[element]].degree <= hub_value) {

          for (reaction in reaction_dictionary) {
            //Return all reactions that contain the entity
            if (checkReaction(
                reaction_dictionary[reaction],
                nn_components[element]) === true) {

              let outputs = checkPlotting(
                data,
                reaction_dictionary[reaction])
              let _max = outputs[0];
              let blocklisted_node = outputs[1];

              if (_max > hub_value) {
                console.log("Not plotting reaction: " + reaction_dictionary[reaction]["id"] + "; entity with too many connections.")
              } else if (blocklisted_node === true) {
                console.log("Not plotting reaction: " + reaction_dictionary[reaction]["id"] + "; contains blocklisted entity.")
              } else {
                components.push(reaction_dictionary[reaction]["id"]);

                for (x in reaction_dictionary[reaction]["reactants"]) {
                  if (
                    data.nodes[reaction_dictionary[reaction]["reactants"][x]].degree <= hub_value
                  ) {
                    components.push(reaction_dictionary[reaction]["reactants"][x]);
                  }
                }
                for (x in reaction_dictionary[reaction]["products"]) {
                  data.nodes[reaction_dictionary[reaction]["products"][x]]
                  if (
                    data.nodes[reaction_dictionary[reaction]["products"][x]].degree <= hub_value
                  ) {
                    components.push(reaction_dictionary[reaction]["products"][x]);
                  }
                }
                for (x in reaction_dictionary[reaction]["modifiers"]) {
                  if (
                    data.nodes[reaction_dictionary[reaction]["modifiers"][x][0]].degree <= hub_value
                  ) {
                    components.push(
                      reaction_dictionary[reaction]["modifiers"][x][0]
                    );
                  }
                }
                for (x in reaction_dictionary[reaction]["additional_components"]) {
                  if (
                    data.nodes[reaction_dictionary[reaction]["additional_components"][x]].degree <= hub_value
                  ) {
                    components.push(
                      reaction_dictionary[reaction]["additional_components"][x]
                    );
                  }
                }
              }
            }
          }
        }
      }

      // Only combine the lists if they pass the threshold
      escape = escape + components.length;
      if (escape > max_nodes) {
        console.log(escape);
        try {
          document.getElementById("warning_line_1").innerHTML =
            '<i class="red-text">Too many entities to plot. Will only plot first ' +
            n +
            " neighborhood(s)</i>";
          document.getElementById("warning_line_2").innerHTML = "";
          alert("Too many entities to plot. Will only plot first " + n + " neighborhood(s)")
        } catch(e) {}
        n = kNN + 2;
      } else {
        nn_components = nn_components.concat(components);
        try {
          document.getElementById("warning_line_1").innerHTML = "<br>";
          document.getElementById("warning_line_2").innerHTML = "<br><br>";
        } catch(e) {}
      }
      n = n + 1;
    }
  }

  var new_elements = get_nodes_links(data, nn_components);
  return new_elements;
}

function kNN_input(d) {
  knn_value = document.getElementById("kNN_button").value;
  console.log("k-NN parameter now set to: " + knn_value);
  change();
}

function stat_input(d) {
  stat_value = document.getElementById("stat_button").value;
  if ((stat_value === null) || (stat_value === "")) {
    stat_value = 1;
  }
  console.log("Stat threshold now set to: " + stat_value);
  change();
}

function hub_input(d) {
  hub_value = document.getElementById("hub_button").value;
  if ((hub_value === null) || (hub_value === "")) {
    hub_value = 1000000;
  }
  console.log("Hub threshold now set to: " + hub_value);
  change();
}

function get_link(d) {
  if (d.type === "complex_component") {
    if (
      d.sub_type === "metabolite_component" ||
      d.sub_type === "protein_component"
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
    display_analytes_dict,
    display_reactions_dict,
    selector,
    _width,
    _height,
    global_motifs) {

  console.log(new_nodes)

  // Final update to prevent plotting of blocklisted nodes
  id_blocklist = [];
  if (graph_genes === false) {
    for (n in new_nodes) {
      if (new_nodes[n].type === "gene_component") {
        id_blocklist.push(new_nodes[n].id)
      }
    }
  }

  node_keep = [];
  for (n in new_nodes) {
    if (data.metadata.blocklist.split(",").includes(new_nodes[n].name)) {
      id_blocklist.push(new_nodes[n].id);
    } else if (id_blocklist.includes(new_nodes[n].id)) {

    } else {
      node_keep.push(new_nodes[n]);
    }
  }

  link_keep = [];
  for (l in new_links) {

    if (id_blocklist.includes(new_links[l].target) || id_blocklist.includes(new_links[l].source)) {
    } else if (id_blocklist.includes(new_links[l].target.id) || id_blocklist.includes(new_links[l].source.id)) {
    } else {
      link_keep.push(new_links[l]);
    }
  }

  graph_nodes = node_keep;
  graph_links = link_keep;

  // Restart graph
  d3.selectAll("#svg_viewer_id").remove();

  if (timecourse === true) {
    var sample = d3.select("circle#dot").attr("x");
    if (sample === null) {
      sample = 0;
    }
  } else {
    var sample = 0;
  }

  console.log("Building graph for sample: ", sample);
  console.log("Hub threshold set at: ", hub_value);

  // Initialize force graph object
  var svg_viewer = d3
    .select(selector)
    .append("svg")
    .attr("id", "svg_viewer_id")
    .attr("width", _width - 5)
    .attr("height", _height - 80)
    .call(
      d3.zoom().on("zoom", function() {
        svg_viewer.attr("transform", d3.event.transform);
      })
    )
    .on("dblclick.zoom", null)
    .append("g");

  const forceX = d3.forceX(_width / 2).strength(0.015);
  const forceY = d3.forceY(_height / 2).strength(0.015);

  const simulation = d3
    .forceSimulation(graph_nodes)
    .force(
      "link",
      d3
        .forceLink(graph_links)
        .id(d => d.id)
        .distance(40)
        .strength(1)
    )
    .force("charge", d3.forceManyBody().strength(-1000))
    .force("center", d3.forceCenter(_width / 2, _height / 2));

  simulation
    .alphaTarget(0.01)
    .alphaMin(0.1)
    .velocityDecay(0.7);

  // Generate edges with style attributes
  svg_viewer
    .append("defs")
    .selectAll("marker")
    .data([
      "collapsed",
      "collapsed_catalyst",
      "collapsed_inhibitor",
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
    .attr("markerHeight", 10)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0, -5L10, 0L0, 5");

  var offset = 30;
  function getGroup(d, links) {

    if (d.type === "gene_component") {
      return "none"
    } else if (d.compartment === "none") {
      for (l in links) {
        if (d.id === links[l].source.id && links[l].target.compartment !== "none") {
          return links[l].target.compartment;
        } else if (d.id === links[l].target.id && links[l].source.compartment !== "none") {
          return links[l].source.compartment;
        } else {}
      }
      return "none";

    } else {
      return d.compartment;
    }
  }

  var curve = d3.line()
    .curve(d3.curveBasis);

  function drawCluster(d) {
    return curve(d.path); // 0.8
  }

  var hullg = svg_viewer.append("g")

  var link = svg_viewer
    .append("g")
    .selectAll("path")
    .data(graph_links)
    .enter()
    .append("path")
    .attr("class", function(d) {
      return "link " + get_link(d);
    })
    .attr("marker-end", function(d) {
      return "url(#" + get_link(d) + ")";
    });

  var node = svg_viewer
    .selectAll(".node")
    .data(graph_nodes)
    .enter()
    .append("g")
    .attr("class", "node")
    .attr("id", function(d) {return d.id})
    .call(
      d3
        .drag()
        .subject(dragsubject)
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended)
    );

  var circle = node
    .append("path")
    .style("fill", function(d) {
      return "rgba(" + d[entity][sample].toString() + ")";
     })
    .style("stroke", "black")
    .style("stroke-width", function(d) {
      if ((d['stats'][sample] === undefined) || (d['stats'][sample] === null)) {
        return 1;
      } else if (d['stats'][sample] < stat_value) {
        return 2;
      } else {
        return 1;
      }
    })
    .style("stroke-dasharray", function(d) {
      if (d.inferred === "true" || d.type === "collapsed") {
        return "2,2";
      } else {
        return "none";
      }
    })
    .attr("d", d3.symbol()
      .size(function(d) {
        return 175;
      })
      .type(function(d) {
        if (d.type === "gene_component" || d.sub_type === "gene") {
          return d3.symbolTriangle;
        } else if (d.complex === "true") {
          return d3.symbolSquare;
        } else if (d.sub_type === "protein_component") {
          return d3.symbolDiamond;
        }  else if (d.type === "metabolite_component" || d.sub_type === "metabolite_component") {
          return d3.symbolCircle;
        } else if (d.type === "reaction" || d.sub_type === "reaction") {
          return d3.symbolStar;
        } else {
          return d3.symbolStar;
        }
      }))
      .attr("id", function(d) {return d.id})

  // Motif page will only send the current time-point's motifs
  var page_path = window.location.pathname;
  var page_name = page_path.substring(page_path.lastIndexOf('/') + 1);
  if (page_name === "motif.html") {
    if (global_motifs !== undefined) {
      let motif_ids = [];
      for (key in global_motifs) {
        motif_ids.push(global_motifs[key].id)
      }
      if (motif_ids.length > 0) {
        graph_nodes.forEach(node=>{
          let rxn_id = node.id;
          if (motif_ids.includes(rxn_id)) {
            d3.selectAll("path#" + rxn_id)
            .style("stroke", "purple")
            .style("stroke-width", 3)
            .attr("d", d3.symbol()
              .size(function(d) {
                return 400;
              })
              .type(d3.symbolStar))
          }
        })
      }
    }
  } else {
    if (global_motifs[sample] !== undefined) {
      if (global_motifs[sample].length > 0) {
        graph_nodes.forEach(node=>{
          let rxn_id = node.id;
          if (global_motifs[sample].includes(rxn_id)) {
            d3.selectAll("path#" + rxn_id)
              .style("stroke", "purple")
              .style("stroke-width", 3)
              .attr("d", d3.symbol()
                .size(function(d) {
                  return 400;
                })
                .type(d3.symbolStar))
          }
        })
      }
    }
  }

  function getCategories(nodes) {

    var categories = new Set();
    for (n in nodes) {
      if (nodes[n].compartment === "none" || nodes[n].compartment === null || nodes[n].type === "gene_component" || nodes[n].compartment === "undefined" || nodes[n].compartment === undefined) {
        categories.add("none");
      } else {
        categories.add(nodes[n].compartment);
      }
    }

    var category_key = new Object();
    var counter = 0;
    categories.forEach( s => {
      category_key[s] = counter;
      counter = counter + 1;
    })

    return category_key
  }

  var categories = getCategories(graph_nodes);
  var fill = d3.schemeCategory10;
  var fill2 = d3.schemeTableau10;

  hullg.selectAll("path.hull").remove();
  hull = hullg
    .selectAll("path.hull")
      .data(convexHulls(graph_nodes, graph_links, getGroup, offset))
    .enter().append("path")
      .attr("class", "hull")
      .attr("d", drawCluster)
      .style("fill", function(d) {
        if (d.group !== "undefined" && d.group !== "none") {
          if (categories[d.group] > 9) {
            return fill2[categories[d.group] - 10];
          } else {
            return fill[categories[d.group]];
          }
        }
      }
    )

    let compartment_dictionary = {};
    for (let _n in graph_nodes) {
      if (graph_nodes[_n]['compartment'] !== undefined) {
        compartment_dictionary[graph_nodes[_n]['compartment']] = graph_nodes[_n]['compartment_display']
      }
    }
    d3.select("button#compartment_legend")
      .on("mouseover", function(d) {
        let category_number = 0;
        let make_string = "";
        for (let s in categories) {
          if (s !== null && s !== "none" && s !== undefined && s !== "undefined" && compartment_dictionary[s] !== undefined) {
            make_string = make_string + "&nbsp;";
            make_string = make_string + "<span class='ellipse' style='--dot_color:" + hull._groups[0][categories[s]].style.fill + ";'></span>";
            make_string = make_string + "&nbsp;&nbsp;&nbsp;&nbsp;";
            make_string = make_string + compartment_dictionary[s];
            make_string = make_string + "</br>";
            category_number = category_number + 1;
          }
        }
        // add one for formatting
        if (category_number > 15) {
          category_number = category_number + 1;
        }

        div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", (60 + (15 * category_number * 1.2)).toString() + "px")
          .style("width", "275px");
        div
          .html("<div style='margin-left:15px;margin-top:15px;'><font size='3'><b><u>Compartments</u></b></font></br></br>" + make_string)
        }
      )
      .on("mouseout", function(d) {
        div.style("opacity", 0);
        div.html("")
        category_number = 0;
      });

  var timer = 0;
  var delay = 200;
  var prevent = false;

  circle
    .on("click", function(d) {
      timer = setTimeout(function() {
        if (!prevent) {
          document.getElementById("reaction_notes").innerHTML = "";

          if (d.synonyms.length > 0) {

            let synonym_string = "";
            for (let s in d.synonyms) {
              synonym_string = synonym_string + "<br /> - " + d.synonyms[s]
            }

            document.getElementById("reaction_notes").innerHTML =
              "<b><i>" + d.name + " Synonyms</i></b>: " + synonym_string;
          }
        }
        prevent = false;
      }, delay);
    })
    .on("dblclick", function(d) {
      clearTimeout(timer);
      prevent = true;
      document.getElementById("reaction_notes").innerHTML = "";

      if (
        type_dict[d["name"]] === "reaction" ||
        type_dict[d["name"]] === "collapsed"
      ) {
        document.getElementById("reaction_notes").innerHTML =
          "<b><i>" + d.name + "</i></b>: " + d.notes;
      } else {
        document.getElementById("type_selection_type").innerHTML =
          "Nearest Neighbor";

        var mod_selection = determineWidth(d["name"]);
        document.getElementById("type_selection").innerHTML = mod_selection;

        if (data.metadata.transcriptomics !== null) {
          graph_genes = true;
        } else {
          graph_genes = false;
        }
        nearest_neighbors(data, entity_id_dict[d["name"]]);
      }
    });

  var text = node
    .append("text")
    .attr("id", function(d) {return d.id})
    .html(function(d) {
      if (type_dict[d.name] === "reaction") {
        // Label other nodes with expression value in parentheses
        return (
          "<tspan dx='16' y='.31em' class='bold-text'>"
          + d.name
          + "</tspan>"
          + "<tspan x='16' y='1.7em'>Compartment: "
          + d.compartment_display
          + "</tspan>"
        );
      } else if (type_dict[d.name] === "collapsed") {
        return (
          "<tspan dx='16' y='.31em' class='bold-text'>"
          + d.name
          + "</tspan>"
        );
      } else {
        if (d.values[sample] === null
        && d.stats[sample] === null) {
          return (
            "<tspan dx='16' y='0em' class='bold-text'>"
            + d.name
            + "</tspan>"
          );
        } else {
          let display_stat;
          if (parseFloat(d.stats[sample]) < 0.01) {
            display_stat = "< 0.01"
          } else {
            display_stat = parseFloat(d.stats[sample]).toFixed(2)
          }
          return (
            "<tspan dx='16' y='-.5em' class='bold-text'>"
            + d.name
            + "</tspan>"
            + "<tspan x='16' y='.7em'>Value: "
            + parseFloat(d.values[sample]).toFixed(2)
            + "</tspan>"
            + "<tspan x='16' y='1.7em'>Statistic: "
            + display_stat
            + "</tspan>"
          );
        }
      }
    });

  simulation.on("tick", tick);

  // Toggle compartment view
  toggle_comp = true;
  d3.select("#toggleCompartments").on("click", function() {
    if (toggle_comp === false) {
      toggle_comp = true;
      hull = hullg
        .selectAll("path.hull")
          .data(convexHulls(graph_nodes, graph_links, getGroup, offset))
        .enter().append("path")
          .attr("class", "hull")
          .attr("d", drawCluster)
          .style("fill", function(d) {
            if (d.group !== "undefined" && d.group !== "none") {
              if (categories[d.group] > 9) {
                return fill2[categories[d.group]];
              } else {
                return fill[categories[d.group]];
              }
            }
          })

    } else {
      toggle_comp = false;
      hullg.selectAll("path.hull").remove();
  }
});

  d3.select("#saveGraph").on("click", function() {

    saveSVG.saveSvgAsPng(
      d3.select("#svg_viewer_id")._groups[0][0],
      "plot.png",
      {
        encoderOptions: 1,
        scale: 10,
        encoderType: "image/png"
      }
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
            "<tspan dx='16' y='.31em' class='bold-text'>" +
            d.name +
            "</tspan>"
            + "<tspan x='16' y='1.7em'>Compartment: "
            + d.compartment_display
            + "</tspan>"
          );
        } else {
          // Label other nodes with expression value in parentheses
          if (d.values[sample] === null
          && d.stats[sample] === null) {
            return (
              "<tspan dx='16' y='0em' class='bold-text'>"
              + d.name
              + "</tspan>"
            );
          } else {
            let display_stat;
            if (parseFloat(d.stats[sample]) < 0.01) {
              display_stat = "< 0.01"
            } else {
              display_stat = parseFloat(d.stats[sample]).toFixed(2)
            }
            return (
              "<tspan dx='16' y='-.5em' class='bold-text'>"
              + d.name
              + "</tspan>"
              + "<tspan x='16' y='.7em'>Value: "
              + parseFloat(d.values[sample]).toFixed(2)
              + "</tspan>"
              + "<tspan x='16' y='1.7em'>Statistic: "
              + display_stat
              + "</tspan>"
            );
          }
        }
      });
    } else {
      toggle_e = false;
      text.html(function(d) {
        return (
          "<tspan dx='16' y='.31em' class='bold-text'>" +
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
    } else {
      graph_genes = false;
    }
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
      _height,
      global_motifs)
  });

  d3.select("#collapseNodes").on("click", function() {

    if (collapse_reactions === false) {
      new_nodes = collapsed_nodes;
      new_links = collapsed_links;

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
        _height,
        global_motifs)
      collapse_reactions = true;
    } else {
      collapsed_nodes = new_nodes;
      collapsed_links = new_links;

      var selection = document.getElementById("pathwayMenu").value;

      for (x in data.pathway_dictionary) {
        if (data.pathway_dictionary[x]['name'] === selection) {
          let reactions = data.pathway_dictionary[x]["reactions"];

          var newer_elements = parse_pathway(
            data,
            reactions,
            data.reaction_dictionary,
            data.degree_dictionary);

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
          var new_display_analytes_dict = new_node_elements[2];
          var new_display_reactions_dict = new_node_elements[3];
          var new_entity_id_dict = new_node_elements[4];

          make_graph(
            data,
            newer_nodes,
            newer_links,
            new_type_dict,
            new_node_dict,
            new_entity_id_dict,
            new_display_analytes_dict,
            new_display_reactions_dict,
            selector,
            _width,
            _height,
            global_motifs)
          collapse_reactions = false;
        } else {
        }
      }
    }
  });

  var cell = node.append("path").attr("class", "cell");

  function convexHulls(nodes, links, index, offset) {
    // Function adapted from: http://bl.ocks.org/GerHobbelt/3071239

    var hulls = {};

    // create point sets
    for (var k = 0; k < nodes.length; ++k) {
      var n = nodes[k];
      if (n.size) continue;
      var i = index(n, links),
          l = hulls[i] || (hulls[i] = []);
      l.push([n.x-offset, n.y-offset]);
      l.push([n.x-offset, n.y+offset]);
      l.push([n.x+offset, n.y-offset]);
      l.push([n.x+offset, n.y+offset]);
    }

    // create convex hulls
    var hullset = [];
    for (i in hulls) {

      var hull_space = d3.polygonHull(hulls[i]);
      if (i !== undefined) {
        hullset.push({
          group: i,
          path: hull_space,
          centroid: d3.polygonCentroid(hull_space)
        });
      }
    }

    return hullset;
  }

  // Draw curved edges
  function tick() {
    link
      .attr("d", linkArc);
    circle
      .attr("transform", transform);
    text
      .attr("transform", transform);
    hull
      .data(convexHulls(graph_nodes, graph_links, getGroup, offset))
      .attr("d", drawCluster);
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

  if (data.metadata.transcriptomics !== "None") {
    graph_genes = true;
  } else {
    graph_genes = false;
  }
  collapse_reactions = true;

  let current_pathway = get_session_info("current_pathway");
  if ((current_pathway !== null) && (current_pathway !== "null")) {
    var selection = data.mod_collapsed_pathways[current_pathway].name;
    var superSelection = "All pathways";
  } else {
    var selection = document.getElementById("pathwayMenu").value;
    var superSelection = document.getElementById("superPathwayMenu").value;
  }

  document.getElementById("reaction_notes").innerHTML = "";

  var mod_selection = determineWidth(selection);
  document.getElementById("type_selection_type").innerHTML = "Pathway";
  document.getElementById("type_selection").innerHTML = mod_selection;
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

  if (superSelection !== "All entities") {
    // Run normal first plot
    var reactions = collapsed_pathway_dict[selection]["reactions"];
    var elements = parse_pathway(
      data,
      reactions,
      data.collapsed_reaction_dictionary,
      data.degree_dictionary);

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
      _height,
      global_motifs);
  } else {
    // Get kNN of entity selected
    var mod_selection = determineWidth(selection);
    document.getElementById("type_selection").innerHTML = mod_selection;
    document.getElementById("type_selection_type").innerHTML =
      "Nearest Neighbor";

    var entity_dictionary = parseEntities(data.nodes);
    nearest_neighbors(data, entity_dictionary[selection]);
  }
}

function test() {
  var assert = require('assert')

  let test_data = {
    'nodes': {
      'N1':{'id':'n1'},
      'N2':{'id':'n2'},
      'N3':{'id':'n3'},
      'N4':{'id':'n4'},
      'N5':{'id':'n5'},
      'N6':{'id':'n6'},
      'N7':{'id':'n7'},
      'N8':{'id':'n8'},
      'N9':{'id':'n9'},
      'N10':{'id':'n10'},
      'N11':{'id':'n11'},
    },
    'links': {
      'N1,N2':{'id':'l1'},
      'N2,N3':{'id':'l2'},
      'N3,N4':{'id':'l3'},
      'N4,N5':{'id':'l4'},
      'N5,N6':{'id':'l5'},
      'N6,N1':{'id':'l6'},
      'N6,N3':{'id':'l7'},
      'N3,N7':{'id':'l8'},
      'N7,N8':{'id':'l9'},
      'N8,N9':{'id':'l10'},
      'N9,N10':{'id':'l11'},
      'N10,N11':{'id':'l12'},
    }
  };
  let test_components = ['N1','N4','N5'];

  describe('graph.js', function () {
    // get_nodes_links()
    describe('get_nodes_links()', function () {
      it('should return 3 nodes and 1 link', function () {

        let test_items = get_nodes_links(
            test_data,
            test_components
        )
        let test_nodes = test_items[0];
        let test_links = test_items[1];
        assert(test_nodes[0].id === 'n1')
        assert(test_nodes[1].id === 'n4')
        assert(test_nodes[2].id === 'n5')
        assert(test_links[0].id === 'l4')

      })
    })
    // parse_kNN_pathway()
    describe('parse_kNN_pathway()', function () {
      it('should return ...', function () {

        let test_data_kNN = {
          'nodes': {
            'R1':{
              'degree':5,
              'name':'r1',
              'id':'R1'},
            'R2':{
              'degree':5,
              'name':'r2',
              'id':'R2'},
            'R3':{
              'degree':5,
              'name':'r3',
              'id':'R3'},
            'R4':{
              'degree':5,
              'name':'r4',
              'id':'R4'},
            'R5':{
              'degree':5,
              'name':'r5',
              'id':'R5'},
            'R6':{
              'degree':5,
              'name':'r6',
              'id':'R6'},
            'N1':{
              'degree':1000,
              'name':'n1',
              'id':'N1'},
            'N2':{
              'degree':10,
              'name':'n2',
              'id':'N2'},
            'N3':{
              'degree':50,
              'name':'n3',
              'id':'N3'},
            'N4':{
              'degree':9,
              'name':'n4',
              'id':'N4'},
            'N5':{
              'degree':10,
              'name':'n5',
              'id':'N5'},
            'N6':{
              'degree':11,
              'name':'n6',
              'id':'N6'},
            'N7':{
              'degree':10,
              'name':'n7',
              'id':'N7'},
            'N8':{
              'degree':5,
              'name':'n8',
              'id':'N8'},
            'N9':{
              'degree':4,
              'name':'n9',
              'id':'N9'},
            'N10':{
              'degree':3,
              'name':'n10',
              'id':'N10'},
            'N11':{
              'degree':7,
              'name':'n11',
              'id':'N11'}
          },
          'links': {
            'N7,R1':{'id':'l1'},
            'R1,N8':{'id':'l2'},
            'N3,R1':{'id':'l3'},
            'N5,R2':{'id':'l4'},
            'N6,R2':{'id':'l5'},
            'R2,N7':{'id':'l6'},
            'N11,R2':{'id':'l7'},
            'N4,R3':{'id':'l8'},
            'R3,N5':{'id':'l9'},
            'R3,N6':{'id':'l10'},
            'N11,R3':{'id':'l11'},
            'N8,R4':{'id':'l12'},
            'N9,R4':{'id':'l13'},
            'R4,N1':{'id':'l14'},
            'R4,N2':{'id':'l15'},
            'N2,R5':{'id':'l16'},
            'R5,N3':{'id':'l17'},
            'N3,R6':{'id':'l18'},
            'R4,R6':{'id':'l19'},
            'R6,N5':{'id':'l20'},
            'N11,R6':{'id':'l21'},
            'N1,R6':{'id':'l22'},
          },
          'reaction_dictionary': {
            'R1': {
              'id': 'R1',
              'reactants': ['N7'],
              'products': ['N8'],
              'modifiers': [['N3','catalyst']],
              'additional_components':[]
            },
            'R2': {
              'id': 'R2',
              'reactants': ['N5','N6'],
              'products': ['N7'],
              'modifiers': [['N11','catalyst']],
              'additional_components':[]
            },
            'R3': {
              'id': 'R3',
              'reactants': ['N4'],
              'products': ['N5','N6'],
              'modifiers': [['N11','inhibitor']],
              'additional_components':[]
            },
            'R4': {
              'id': 'R4',
              'reactants': ['N8','N9'],
              'products': ['N1','N2'],
              'modifiers': [],
              'additional_components':[]
            },
            'R5': {
              'id': 'R5',
              'reactants': ['N2'],
              'products': ['N3'],
              'modifiers': [],
              'additional_components':[]
            },
            'R6': {
              'id': 'R6',
              'reactants': ['N3','N4'],
              'products': ['N5'],
              'modifiers': [['N11','catalyst'],['N1','inhibitor']],
              'additional_components':[]
            }
          },
          'collapsed_reaction_dictionary': {
            'R1': {
              'id': 'R1',
              'reactants': ['N7'],
              'products': ['N8'],
              'modifiers': [['N3','catalyst']],
              'additional_components':[]
            },
            'R2': {
              'id': 'R2',
              'reactants': ['N5','N6'],
              'products': ['N7'],
              'modifiers': [['N11','catalyst']],
              'additional_components':[]
            },
            'R3': {
              'id': 'R3',
              'reactants': ['N4'],
              'products': ['N5','N6'],
              'modifiers': [['N11','inhibitor']],
              'additional_components':[]
            },
            'R4': {
              'id': 'R4',
              'reactants': ['N8','N9'],
              'products': ['N1','N2'],
              'modifiers': [],
              'additional_components':[]
            },
            'R5': {
              'id': 'R5',
              'reactants': ['N2'],
              'products': ['N3'],
              'modifiers': [],
              'additional_components':[]
            },
            'R6': {
              'id': 'R6',
              'reactants': ['N3','N4'],
              'products': ['N5'],
              'modifiers': [['N11','catalyst'],['N1','inhibitor']],
              'additional_components':[]
            }
          },
          'degree_dictionary': {
            'N1':1000,
            'N2':10,
            'N3':50,
            'N4':9,
            'N5':10,
            'N6':11,
            'N7':10,
            'N8':5,
            'N9':4,
            'N10':3,
            'N11':7,
            'R1':5,
            'R2':5,
            'R3':5,
            'R4':5,
            'R5':5,
            'R6':5,
          },
          'blocklist': ['n2'],

        }
        let test_entity_id = "N7";
        let test_kNN = 1;
        let test_new_elements = parse_kNN_pathway(
          test_data_kNN,
          test_entity_id,
          test_kNN)
        let el1 = test_new_elements[0];
        let el2 = test_new_elements[1];
        for (el in el1) {
          if (el1[el].name === 'r3'
              || el1[el].name === 'r4'
              || el1[el].name === 'r5'
              || el1[el].name === 'r6') {
            assert(false)
          }
        }
        test_kNN = 2;
        let test_new_elements2 = parse_kNN_pathway(
          test_data_kNN,
          test_entity_id,
          test_kNN)
        let el1_2 = test_new_elements2[0];
        let el2_2 = test_new_elements2[1];
        for (el in el1_2) {
          if (el1_2[el].name === 'r4'
              || el1_2[el].name === 'r5') {
            assert(false)
          }
        }
      })
    })
    // checkReaction()
    describe('checkReaction()', function () {
      it('should return true and false, respectively', function () {
        let test_reaction_1 = {
          'reactants':['N1'],
          'products':['N2','N3','N4'],
          'modifiers':[['N5','catalyst'],['N6','inhibitor']],
          'additional_components':['G1']
        }
        let test_output;
        test_output_reaction_1 = checkReaction(
          test_reaction_1,
          'N1'
        )
        assert(test_output_reaction_1 === true)
        let test_output_reaction_2 = checkReaction(
          test_reaction_1,
          'N6'
        )
        assert(test_output_reaction_2 === false)
      })
    })
    // get_link()
    describe('get_link()', function () {
      it('should return metabolite_component, not_complex_component, and protein_component, respectively', function () {
        let test_link1 = {
          'type':'complex_component',
          'sub_type':'metabolite_component'
        }
        let test_out1 = get_link(test_link1)
        assert(test_out1 === 'metabolite_component')
        let test_link2 = {
          'type':'not_complex_component',
          'sub_type':'ppp_component'
        }
        let test_out2 = get_link(test_link2)
        assert(test_out2 === 'not_complex_component')
        let test_link3 = {
          'type':'complex_component',
          'sub_type':'protein_component'
        }
        let test_out3 = get_link(test_link3)
        assert(test_out3 === 'protein_component')

      })
    })
  })
}
module.exports = test
