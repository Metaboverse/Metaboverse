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
var d3 = require("d3");

var width = window.innerWidth * 0.8;
var height = window.innerHeight * 0.8;
var timeTransition = 8000;
var numberConditions = 1; // determine based on user data
var start = 0; // data min index
var end = 1; // data max index
var slider_index = 0;
var timecourse = false;

function buildSlider(categories, names) {

  end = Math.max(...categories);
  numberConditions = categories.length;

  var x = d3
    .scaleLinear()
    .domain([start, end])
    .range([90, width - 200])
    .clamp(true);

  var slider = d3
    .select("#slider")
    .append("svg")
    .attr("id", "slide")
    .attr("width", width)
    .attr("height", 50)
    .attr("viewBox", [0, -15, width, 50])
    .attr("overflow", "visible");

  slider
    .append("g")
    .attr("class", "slider")
    .attr("transform", "translate(" + 50 + "," + height / 2 + ")");

  slider
    .append("line")
    .attr("class", "track")
    .attr("id", "track")
    .attr("x1", x.range()[0])
    .attr("x2", x.range()[1])
    .select(function() {
      return this.parentNode.appendChild(this.cloneNode(true));
    })
    .attr("class", "track-inset")
    .select(function() {
      return this.parentNode.appendChild(this.cloneNode(true));
    })
    .attr("class", "track-overlay")
    .call(
      d3
        .drag()
        .on("start.interrupt", function() {
          slider.interrupt();
        })
        .on("start drag", function() {
          hue(x.invert(d3.event.x));
        })
    );

  slider
    .insert("g", ".track-overlay")
    .attr("class", "ticks")
    .attr("transform", "translate(0," + 24 + ")")
    .selectAll("text")
    .data(x.ticks(numberConditions))
    .enter()
    .append("text")
    .attr("id", "tick")
    .attr("x", x)
    .attr("text-anchor", "middle")
    .text(function(d, i) {
      return names[i];
    });

  var handle = slider
    .insert("circle", ".track-overlay")
    .attr("id", "dot")
    .attr("class", "handle")
    .attr("r", 12);

  slider.transition() // Gratuitous intro!
    .duration(1)
    .tween("hue", function() {
      var i = d3.interpolate(end, start);
      return function(t) { hue(i(t)); };
    });

  function hue(h) {

    var slider_index = Math.floor(h + 0.5);
    handle
      .attr("cx", x(slider_index))
      .attr("x", slider_index);

    let node = d3.selectAll("circle")
    node.each(function(d) {

      if (d !== undefined) {

        // if rectangle, ellipse, circle
        // change fill and text
        if (d.sub_type === "metabolite_component") {
          try {
            d3.select("circle#" + d.id)
            .style("--node_color", function(d) {
              return "rgba(" + d["values_js"][slider_index].toString() + ")";
            })
            .style("--node_border", function(d) {
              if ((d.stats[slider_index] === undefined) | (d.stats[slider_index] === null)) {
                return 1;
              } else if (d.stats[slider_index] < stat_value) {
                return 2;
              } else {
                return 1;
              }
            })
            d3.select("text#" + d.id)
              .html(function(d) {
                return (
                  "<tspan dx='16' y='-.5em' style='font-weight: bold;'>"
                  + d.name
                  + "</tspan>"
                  + "<tspan x='16' y='.7em'>Value: "
                  + parseFloat(d.values[slider_index]).toFixed(2)
                  + "</tspan>"
                  + "<tspan x='16' y='1.7em'>Statistic: "
                  + parseFloat(d.stats[slider_index]).toFixed(2)
                  + "</tspan>"
                );
              })
          } catch(err) {}
        }
        if (d.sub_type === "gene") {
          try {
            d3.select("ellipse#" + d.id)
              .style("--node_color", function(d) {
                return "rgba(" + d["values_js"][slider_index].toString() + ")";
              })
              .style("--node_border", function(d) {
                if ((d.stats[slider_index] === undefined) | (d.stats[slider_index] === null)) {
                  return 1;
                } else if (d.stats[slider_index] < stat_value) {
                  return 2;
                } else {
                  return 1;
                }
              })
            d3.select("text#" + d.id)
              .html(function(d) {
                return (
                  "<tspan dx='16' y='-.5em' style='font-weight: bold;'>"
                  + d.name
                  + "</tspan>"
                  + "<tspan x='16' y='.7em'>Value: "
                  + parseFloat(d.values[slider_index]).toFixed(2)
                  + "</tspan>"
                  + "<tspan x='16' y='1.7em'>Statistic: "
                  + parseFloat(d.stats[slider_index]).toFixed(2)
                  + "</tspan>"
                );
              })
          } catch(err) {}
        }
        if (d.sub_type === "protein_component") {
          try {
            d3.select("rect#" + d.id)
              .style("--node_color", function(d) {
                return "rgba(" + d["values_js"][slider_index].toString() + ")";
              })
              .style("--node_border", function(d) {
                if ((d.stats[slider_index] === undefined) | (d.stats[slider_index] === null)) {
                  return 1;
                } else if (d.stats[slider_index] < stat_value) {
                  return 2;
                } else {
                  return 1;
                }
              })
            d3.select("text#" + d.id)
              .html(function(d) {
                return (
                  "<tspan dx='16' y='-.5em' style='font-weight: bold;'>"
                  + d.name
                  + "</tspan>"
                  + "<tspan x='16' y='.7em'>Value: "
                  + parseFloat(d.values[slider_index]).toFixed(2)
                  + "</tspan>"
                  + "<tspan x='16' y='1.7em'>Statistic: "
                  + parseFloat(d.stats[slider_index]).toFixed(2)
                  + "</tspan>"
                );
              })
          } catch(err) {}
        }

        // if reaction and in current motif set, enlarge, if not, reset
        try {
          if (global_motifs !== undefined) {
            if (global_motifs[slider_index] !== undefined) {
              if (global_motifs[slider_index].length > 0) {
                if (global_motifs[slider_index].includes(d.id)) {
                  d3.selectAll("circle#" + d.id)
                    .style("r", "10px")
                    .style("stroke", "purple")
                    .style("--node_border", 5)
                } else {
                  d3.selectAll("circle#" + d.id)
                  .style("stroke", "black")
                  .style("--node_color", function(d) {
                    return "rgba(" + d[entity][slider_index].toString() + ")";
                  })
                  .style("--node_border", function(d) {
                    if ((d.stats[slider_index] === undefined) | (d.stats[slider_index] === null)) {
                      return 1;
                    } else if (d.stats[slider_index] < stat_value) {
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
                }
              }
            }
          }
        } catch(err) {}
      }
    })
  }
}

// Check number of categories
function checkCategories(categories, labels) { //, names) {

  if (categories.length > 1) {
    timecourse = true;
    let names = labels.split(',');
    names = names.map(function (n) {
      return n.trim();
    });
    buildSlider(categories, names);
  } else {
    timecourse = false;
  }
  return timecourse;
}
