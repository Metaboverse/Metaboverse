/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) Jordan A. Berg, The University of Utah

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

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
        var these_motifs;
        if (collapse_reactions === true) {
          these_motifs = collapsed_global_motifs;
        } else {
          these_motifs = global_motifs;
        }
        if (d !== undefined) {
          if (d.type === "reaction" || d.type === "collapsed") {
            // if reaction and in current motif set, enlarge, if not, reset
            if (these_motifs !== undefined) {
              if (these_motifs[slider_index] !== undefined) {
                if (these_motifs[slider_index].length > 0) {
                  let rxn_id = d.id;
                  if (these_motifs[slider_index].includes(rxn_id)) {
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
                return set_significance_weight(d, slider_index, stat_type, stat_value);
              })
            d3.select("text#" + d.id)
              .html(function(d) {
                let this_name;
                if (d.user_label !== undefined && d.sub_type === "metabolite_component") {
                  this_name = d.user_label;
                } else {
                  this_name = d.name;
                }
                if (d.values[slider_index] === null &&
                  d.stats[slider_index] === null) {
                  return (
                    "<tspan dx='16' y='0em' class='bold-text'>" +
                    this_name +
                    "</tspan>"
                  );
                } else {
                  let display_stat;
                  if (parseFloat(d.stats[slider_index]) < 0.01) {
                    display_stat = "< 0.01"
                  } else {
                    display_stat = parseFloat(d.stats[slider_index]).toFixed(2)
                  }
                  let output_stat_string = ("<tspan dx='16' y='-.5em' class='bold-text'>" +
                    this_name +
                    "</tspan>" +
                    "<tspan x='16' y='.7em'>Value: " +
                    parseFloat(d.values[slider_index]).toFixed(2) +
                    "</tspan>");
                  if (stat_type !== "array") {
                    output_stat_string = (output_stat_string + 
                      "<tspan x='16' y='1.7em'>Statistic: " +
                      display_stat +
                      "</tspan>"
                      );
                  }
                  return output_stat_string;
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
