/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2021 Jordan A. Berg
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

const data = {
  nodes: [{
      id: "node1",
      color: [191, 191, 191, 1]
    },
    {
      id: "node2",
      color: [191, 191, 191, 1]
    },
    {
      id: "node3",
      color: [191, 191, 191, 1]
    },
    {
      id: "node4",
      color: [191, 191, 191, 1]
    },
    {
      id: "node5",
      color: [191, 191, 191, 1]
    },
    {
      id: "node6",
      color: [191, 191, 191, 1]
    },

    {
      id: "node7",
      color: [255, 255, 255, 1]
    },
    {
      id: "node8",
      color: [255, 255, 255, 1]
    },
    {
      id: "node9",
      color: [255, 255, 255, 1]
    },
    {
      id: "node10",
      color: [255, 255, 255, 1]
    },
    {
      id: "node11",
      color: [255, 255, 255, 1]
    },
    {
      id: "node12",
      color: [255, 255, 255, 1]
    }
  ],

  links: [{
      source: "node1",
      target: "node2",
      type: "reaction"
    },
    {
      source: "node2",
      target: "node3",
      type: "reaction"
    },
    {
      source: "node3",
      target: "node4",
      type: "reaction"
    },
    {
      source: "node4",
      target: "node5",
      type: "reaction"
    },
    {
      source: "node5",
      target: "node6",
      type: "reaction"
    },
    {
      source: "node6",
      target: "node1",
      type: "reaction"
    },

    {
      source: "node7",
      target: "node1",
      type: "up"
    },
    {
      source: "node8",
      target: "node2",
      type: "down"
    },
    {
      source: "node9",
      target: "node3",
      type: "up"
    },
    {
      source: "node10",
      target: "node4",
      type: "down"
    },
    {
      source: "node11",
      target: "node5",
      type: "up"
    },
    {
      source: "node12",
      target: "node6",
      type: "gene"
    }
  ]
};

height = 600;
width = 850;

const links = data.links;
const nodes = data.nodes;

const forceX = d3.forceX(width / 2).strength(0.015);
const forceY = d3.forceY(height / 2).strength(0.015);

const simulation = d3
  .forceSimulation(nodes)
  .force(
    "link",
    d3
    .forceLink(links)
    .id(d => d.id)
    .distance(40)
    .strength(1)
  )
  .force("charge", d3.forceManyBody().strength(-1000))
  .force("center", d3.forceCenter(width / 2, height / 2));

simulation
  .alphaTarget(0.3)
  .alphaMin(0.1)
  .velocityDecay(0.7);

const svg = d3
  .select("#graph")
  .append("svg")
  .attr("viewBox", [0, 0, width, height]);

svg
  .append("defs")
  .selectAll("marker")
  .data(["reaction", "up", "down", "gene"])
  .enter()
  .append("marker")
  .attr("id", function(d) {
    return d;
  })
  .attr("viewBox", "0 -5 10 10")
  .attr("refX", 18.7)
  .attr("refY", -2)
  .attr("markerWidth", 6)
  .attr("markerHeight", 6)
  .attr("orient", "auto")
  .append("path")
  .attr("d", "M0, -5L10, 0L0, 5");

var link = svg
  .append("g")
  .selectAll("path")
  .data(links)
  .enter()
  .append("path")
  .attr("class", function(d) {
    return "link " + d.type;
  })
  .attr("marker-end", function(d) {
    return "url(#" + d.type + ")";
  });

const node = svg
  .append("g")
  .attr("stroke", "#fff")
  .attr("stroke-width", 1.5)
  .selectAll("circle")
  .data(nodes)
  .join("circle")
  .attr("r", 9)
  .style("--node_color", function(d) {
    return "rgba(" + d.color + ")";
  });
//.call(d3.drag()
//    .on("start", dragstarted)
//    .on("drag", dragged)
//    .on("end", dragended));

node.append("title").text(d => d.id);

simulation.on("tick", tick);

// Draw curved edges
function tick() {
  link.attr("d", linkArc);
  node.attr("transform", transform);
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
