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
    .data(function() {
      let _range = x.ticks(numberConditions);
      let _iter = Math.min(..._range);
      let _max = Math.max(..._range);
      let _step = (_max - _iter) / (numberConditions - 1);
      let _ticks = [];
      while (_iter <= _max) {
        _ticks.push(_iter);
        _iter = _iter + _step;
      }
      return _ticks;
    })
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
      return function(t) {
        hue(i(t));
      };
    });

  function hue(h) {

    var slider_index = Math.floor(h + 0.5);
    handle
      .attr("cx", x(slider_index))
      .attr("x", slider_index);

    let node = d3.selectAll("path")
    node.each(function(d) {

      // if rectangle, ellipse, circle
      // change fill and text
      try {
        if (d !== undefined) {
          if (d.type === "reaction" || d.type === "collapsed") {
            // if reaction and in current motif set, enlarge, if not, reset
            if (global_motifs !== undefined) {
              if (global_motifs[slider_index] !== undefined) {
                if (global_motifs[slider_index].length > 0) {
                  let rxn_id = d.id;
                  if (global_motifs[slider_index].includes(rxn_id)) {
                    d3.select("path#" + rxn_id)
                      .style("stroke", "purple")
                      .style("stroke-width", 3)
                      .attr("d", d3.symbol()
                        .size(function(d) {
                          return 400;
                        })
                        .type(d3.symbolStar))
                  } else {
                    d3.select("path#" + rxn_id)
                      .style("stroke", "black")
                      .style("stroke-width", 1)
                      .attr("d", d3.symbol()
                        .size(function(d) {
                          return 175;
                        })
                        .type(d3.symbolStar))
                  }
                }
              }
            }
          } else {
            d3.select("path#" + d.id)
              .style("fill", function(d) {
                return "rgba(" + d["values_js"][slider_index].toString() + ")";
              })
              .style("stroke", "black")
              .style("stroke-width", function(d) {
                if ((d['stats'][slider_index] === undefined) || (d['stats'][slider_index] === null)) {
                  return 1;
                } else if (d['stats'][slider_index] < stat_value) {
                  return 2;
                } else {
                  return 1;
                }
              })
            d3.select("text#" + d.id)
              .html(function(d) {
                console.log(slider_index)
                if (d.values[slider_index] === null &&
                  d.stats[slider_index] === null) {
                  return (
                    "<tspan dx='16' y='0em' class='bold-text'>" +
                    d.name +
                    "</tspan>"
                  );
                } else {
                  let display_stat;
                  if (parseFloat(d.stats[slider_index]) < 0.01) {
                    display_stat = "< 0.01"
                  } else {
                    display_stat = parseFloat(d.stats[slider_index]).toFixed(2)
                  }
                  return (
                    "<tspan dx='16' y='-.5em' class='bold-text'>" +
                    d.name +
                    "</tspan>" +
                    "<tspan x='16' y='.7em'>Value: " +
                    parseFloat(d.values[slider_index]).toFixed(2) +
                    "</tspan>" +
                    "<tspan x='16' y='1.7em'>Statistic: " +
                    display_stat +
                    "</tspan>"
                  );
                }
              })
          }
        }
      } catch (err) {}
    })
  }
}


// Check number of categories
function checkCategories(categories, names) {

  if (categories.length > 1) {
    timecourse = true;
    names = names.map(function(n) {
      return n.trim();
    });
    buildSlider(categories, names);
  } else {
    timecourse = false;
  }
  return timecourse;
}

function populateExclusions(categories, names) {
  var select = document.getElementById("exclude_type");
  for (let i = categories.length - 1; i >= 0; --i) {
    var option = document.createElement('option');
    option.text = option.value = names[i];
    select.add(option, 0);
  }
  var option = document.createElement('option');
  option.text = option.value = "No exclusion";
  select.add(option, 0);
  $("#exclude_type").val("No exclusion");
}
