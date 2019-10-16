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

d3.json("test_data.json", function(data) {
  
  // Extract links links from data for downstream use
  var links = data[0].links;

  // Set nodes variables up and extract node data
  var nodex = {}
  var nodes = data[0].nodes;

  // Make dictionary of node color values and types
  node_dict = {}
  type_dict = {}
  nodes.forEach(function(node) {

    node_dict[node['name']] = node['rgba_js']
    type_dict[node['name']] = node['type']

  });

  // Populate node object with relevant information and data
  links.forEach(function(link) {

    // Add pertinent node information
    link.source = nodex[link.source] || (nodex[link.source] = {name: link.source});
    link.target = nodex[link.target] || (nodex[link.target] = {name: link.target});
    nodex[link.source["name"]]["color"] = node_dict[link.source["name"]];
    nodex[link.target["name"]]["color"] = node_dict[link.target["name"]];

    // Prioritize reaction nodes in layout
    if (type_dict[link.source["name"]] == "reaction") {
      nodex[link.source["name"]]["weight"] = 100;
    } else {
      nodex[link.source["name"]]["weight"] = 1;
    };

    if (type_dict[link.target["name"]] == "reaction") {
      nodex[link.target["name"]]["weight"] = 100;
    } else {
      nodex[link.target["name"]]["weight"] = 1;
    };

  });

  // Allow flexible window dimensions based on initial window size when opened
  var width = window.innerWidth;
  var height = window.innerHeight;

  var svg = d3.select("body").append("svg")
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

  d3.select("#gravity").on("input", function() {

    console.log(this.value);

  });

  var voronoi = d3.geom.voronoi()
      .x(function(d) { return d.x; })
      .y(function(d) { return d.y; })
      .clipExtent([[0, 0], [width, height]]);

  // Initialize force graph object
  var g_nodes = d3.values(nodex);

  // Build graph
  force
      .nodes(g_nodes)
      .links(links)
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
      .attr("d", "M0,-5L10,0L0,5");

  var path = svg.append("g").selectAll("path")
      .data(force.links())
    .enter().append("path")
      .attr("class", function(d) { return "link " + d.type; })
      .attr("marker-end", function(d) { return "url(#" + d.type + ")"; });

  var node = svg.selectAll(".node")
      .data(g_nodes)
    .enter().append("g")
      .attr("class", "node")
      .style("--node_color", function(d) { return "rgba(" + d.color + ")"; })
      .call(force.drag);

  var circle = node.append("circle")
      .attr("r", 6);

  var text = node.append("text")
      .attr("x", 16)
      .attr("y", ".31em")
      .text(function(d) { return d.name; });

  var cell = node.append("path")
      .attr("class", "cell");

  // Draw curved edges
  function tick() {
    path.attr("d", linkArc);
    circle.attr("transform", transform);
    text.attr("transform", transform);
  }

  function linkArc(d) {

    var dx = d.target.x - d.source.x,
        dy = d.target.y - d.source.y,
        dr = Math.sqrt(dx * dx + dy * dy);

    return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;

  }

  function transform(d) {

    return "translate(" + d.x + "," + d.y + ")";

  }

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

  d3.select("body")
    .on("keydown", function () {

      toggle_zoom = d3.event.altKey;

    });

  d3.select("body")
    .on("keyup", function () {

      toggle_zoom = false;

    });
  // End of adaption of code by Pedro Tabacof

  toggle_a = false
  d3.select("#toggleAnalytes").on("click", function() {
    if (toggle_a === false) {

      toggle_a = true;
      console.log(text)
      var text = svg.append("g").selectAll("text")
        .data(force.nodes())
        .attr("class", "update")
        .enter().append("text")
          .attr("class", "enter")
          .attr("x", function(d, i) { return i * 32; })
          .attr("dy", ".35em")
        .merge(text)
          .text(function(d) { return d.name; });

    } else {

      toggle_a = false;

    }

    console.log(toggle_a)

  });

  toggle_r = false
  d3.select("#toggleReactions").on("click", function() {

    if (toggle_r === false) {

      toggle_r = true;

    } else {

      toggle_r = false;

    }

    console.log('HELLO')

  });

});
