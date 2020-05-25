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

function initialize_nodes(nodes, node_dict, type_dict) {
  var display_analytes_dict = {};
  var display_reactions_dict = {};
  var entity_id_dict = {};

  // Make dictionary of node color values and types
  nodes.forEach(function(node) {
    node_dict[node["name"]] = node["values_js"];
    type_dict[node["name"]] = node["type"];

    entity_id_dict[node["name"]] = node["id"];
    entity_id_dict[node["id"]] = node["name"];

    if ((node["type"] === "reaction") | (node["type"] === "collapsed")) {
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
  document.getElementById("warning_line_1").innerHTML = "<br>";
  document.getElementById("warning_line_2").innerHTML = "<br><br>";

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
    if (
      target_rxns["reactants"].includes(entity_id) ||
      target_rxns["products"].includes(entity_id) ||
      target_rxns["modifiers"].includes(entity_id) ||
      target_rxns["additional_components"].includes(entity_id)
    ) {
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
    var hub_exclusion = new Set();
    if (hub_value > 0) {
      data.nodes.forEach(function(node) {
        if (node.degree > hub_value) {
          hub_exclusion.add(node.id);
          console.log(
            "Filtering " +
              node.id +
              " (" +
              node.degree +
              " degrees) -- may cause edge loss"
          );
        }
      });
    }

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
  knn_value = document.getElementById("kNN_button").value;
  console.log("k-NN parameter now set to: " + knn_value);
}

function stat_input(d) {
  stat_value = document.getElementById("stat_button").value;

  if ((stat_value === null) | (stat_value === "")) {
    stat_value = 1;
  }

  console.log("Stat threshold now set to: " + stat_value);
}

function hub_input(d) {
  hub_value = document.getElementById("hub_button").value;

  if ((hub_value === null) | (hub_value === "")) {
    hub_value = 1000000;
  }

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
    display_analytes_dict,
    display_reactions_dict,
    selector,
    _width,
    _height,
    global_motifs) {

  // Final update to prevent plotting of blacklisted nodes
  node_keep = [];
  id_blacklist = [];
  for (n in new_nodes) {
    if (data.metadata.blacklist.includes(new_nodes[n].name)) {
      id_blacklist.push(new_nodes[n].id);
    } else {
      node_keep.push(new_nodes[n]);
    }
  }

  link_keep = [];
  for (l in new_links) {
    if (id_blacklist.includes(new_links[l].target) | id_blacklist.includes(new_links[l].source)) {
    } else {
      link_keep.push(new_links[l]);
    }
  }

  new_nodes = node_keep;
  new_links = link_keep;

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

  console.log(new_nodes)

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

    if (d.compartment === "none") {
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
    .data(new_links)
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
    .data(new_nodes)
    .enter()
    .append("g")
    .attr("class", "node")
    .attr("id", function(d) {return d.id})
    .style("--node_color", function(d) {
      return "rgba(" + d[entity][sample].toString() + ")";
    })
    .style("--node_border", function(d) {
      if ((d['stats'][sample] === undefined) | (d['stats'][sample] === null)) {
        return 1;
      } else if (d['stats'][sample] < stat_value) {
        return 2;
      } else {
        return 1;
      }
    })
    .style("r", function() {
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
  })
  .attr("id", function(d) {return d.id});

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
        new_nodes.forEach(node=>{
          let rxn_id = node.id;
          if (motif_ids.includes(rxn_id)) {
            d3.selectAll("circle#" + rxn_id)
              .style("r", "10px")
              .style("stroke", "purple")
              .style("--node_border", 5)
          }
        })
      }
    }
  } else {
    if (global_motifs[sample] !== undefined) {
      if (global_motifs[sample].length > 0) {
        new_nodes.forEach(node=>{
          let rxn_id = node.id;
          if (global_motifs[sample].includes(rxn_id)) {
            d3.selectAll("circle#" + rxn_id)
              .style("r", "10px")
              .style("stroke", "purple")
              .style("--node_border", 5)
          }
        })
      }
    }
  }

  function getCategories(nodes) {

    var categories = new Set();
    for (n in nodes) {
      if (nodes[n].compartment == "none") {
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

  var categories = getCategories(new_nodes);
  var fill = d3.schemeCategory10;

  hullg.selectAll("path.hull").remove();
  hull = hullg
    .selectAll("path.hull")
      .data(convexHulls(new_nodes, new_links, getGroup, offset))
    .enter().append("path")
      .attr("class", "hull")
      .attr("d", drawCluster)
      .style("fill", function(d) {
        if (d.group !== "undefined") {
          return fill[categories[d.group]];
        }
      })

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

        graph_genes = true;
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
          "<tspan dx='16' y='.31em' style='font-weight: bold;'>"
          + d.name
          + "</tspan>"
          + "<tspan x='16' y='1.7em'>Compartment: "
          + d.compartment_display
          + "</tspan>"
        );
      } else if (type_dict[d.name] === "collapsed") {
        return (
          "<tspan dx='16' y='.31em' style='font-weight: bold;'>"
          + d.name
          + "</tspan>"
        );
      } else {
        return (
          "<tspan dx='16' y='-.5em' style='font-weight: bold;'>"
          + d.name
          + "</tspan>"
          + "<tspan x='16' y='.7em'>Value: "
          + parseFloat(d.values[sample]).toFixed(2)
          + "</tspan>"
          + "<tspan x='16' y='1.7em'>Statistic: "
          + parseFloat(d.stats[sample]).toFixed(2)
          + "</tspan>"
        );
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
          .data(convexHulls(new_nodes, new_links, getGroup, offset))
        .enter().append("path")
          .attr("class", "hull")
          .attr("d", drawCluster)
          .style("fill", function(d) {
            if (d.group !== "undefined") {
              return fill[categories[d.group]];
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
            parseFloat(d.values[sample]).toFixed(2) +
            "</tspan>" +
            "<tspan x='16' y='1.7em'>Statistic: " +
            parseFloat(d.stats[sample]).toFixed(2) +
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
    if (entity === "values_js") {
      entity = "stats_js";
      node.style("--node_color", function(d) {
        return "rgba(" + d[entity][sample].toString() + ")";
      });
    } else {
      entity = "values_js";
      node.style("--node_color", function(d) {
        return "rgba(" + d[entity][sample].toString() + ")";
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
        display_analytes_dict,
        display_reactions_dict,
        selector,
        _width,
        _height,
        global_motifs)
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
          if (data.degree_dictionary[new_nodes[x]["id"]] <= hub_value) {
            newer_components.push(new_nodes[x]["id"]);
          }
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
    }
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
      .data(convexHulls(new_nodes, new_links, getGroup, offset))
      .attr("d", drawCluster);

    /*Revisit labeling hull later
    hull_text
      .attr("transform", function(d) {
        return "translate(" + d.centroid + ")";
      });
    */

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

  let current_pathway = get_session_info("current_pathway");
  if ((current_pathway !== null) & (current_pathway !== "null")) {
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

    graph_genes = true;
    var entity_dictionary = parseEntities(data.nodes);
    nearest_neighbors(data, entity_dictionary[selection]);
  }
}
