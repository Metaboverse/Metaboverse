var data = [
    {
        "directed": true,
        "multigraph": false,
        "graph": {},
        "nodes": [
            {
                "name": "node1",
                "id": "node1",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "node2",
                "id": "node2",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "node3",
                "id": "node3",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "node4",
                "id": "node4",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "node5",
                "id": "node5",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "node6",
                "id": "node6",
                "type": "reaction",
                "color": "grey",
                "rgba_js": [
                    191,
                    191,
                    191,
                    1
                ],
            },
            {
                "name": "mod1",
                "id": "mod1",
                "type": "down",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
            {
                "name": "mod2",
                "id": "mod2",
                "type": "up",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
            {
                "name": "mod3",
                "id": "mod3",
                "type": "down",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
            {
                "name": "mod4",
                "id": "mod4",
                "type": "up",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
            {
                "name": "mod5",
                "id": "mod5",
                "type": "down",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
            {
                "name": "mod6",
                "id": "mod6",
                "type": "up",
                "color": "white",
                "rgba_js": [
                    255,
                    255,
                    255,
                    1
                ],
            },
          ],
        "links": [
            {
                "type": "reaction",
                "source": "node1",
                "target": "node2",
                "color": "grey",
            },
            {
                "type": "reaction",
                "source": "node2",
                "target": "node3",
                "color": "grey",
            },
            {
                "type": "reaction",
                "source": "node3",
                "target": "node4",
                "color": "grey",
            },
            {
                "type": "reaction",
                "source": "node4",
                "target": "node5",
                "color": "grey",
            },
            {
                "type": "reaction",
                "source": "node5",
                "target": "node6",
                "color": "grey",
            },
            {
                "type": "reaction",
                "source": "node6",
                "target": "node1",
                "color": "grey",
            },
            {
                "type": "down",
                "source": "mod1",
                "target": "node1",
                "color": "red",
            },
            {
                "type": "up",
                "source": "mod2",
                "target": "node2",
                "color": "green",
            },
            {
                "type": "down",
                "source": "mod3",
                "target": "node3",
                "color": "red",
            },
            {
                "type": "up",
                "source": "mod4",
                "target": "node4",
                "color": "purple",
            },
            {
                "type": "gene",
                "source": "mod5",
                "target": "node5",
                "color": "purple",
            },
            {
                "type": "up",
                "source": "mod6",
                "target": "node6",
                "color": "green",
            },

          ]
        }
      ];

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

});

// Allow flexible window dimensions based on initial window size when opened
var width = window.innerWidth;
var height = window.innerHeight / 2;

// Initialize force graph object
var force = d3.layout.force()
    .nodes(d3.values(nodex))
    .links(links)
    .size([width, height])
    .linkDistance(40)
    .charge(-1000)
    .on("tick", tick)
    .start();

setInterval(function(){force.alpha(0.075);},1);

// Initialize plot area
var svg = d3.select("#graph").append("svg")
    .attr("width", width)
    .attr("height", height)
  .append("g");

// Generate edges with style attributes
svg.append("defs").selectAll("marker")
    .data([
      "reaction",
      "up",
      "down",
      "gene"])
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

// Generate nodes with style attributes and colors based on expression from JSON graph file
var circle = svg.append("g").selectAll("circle")
    .data(force.nodes())
    .enter().append("circle")
    .style("--node_color", function(d) { return "rgba(" + d.color + ")"; })
    .attr("r", 6)
    .call(force.drag)
    .on('mouseover', function(d, i) {

      d3.select(this).attr({
          r: 12
        })
      })
    .on('mouseout', function(d, i) {

      d3.select(this).attr({
          r: 6
        })
      })

// Draw curved edges
function tick() {
  path.attr("d", linkArc);
  circle.attr("transform", transform);
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
