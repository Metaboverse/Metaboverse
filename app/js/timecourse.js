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
var width = window.innerWidth * 0.7;
var height = window.innerHeight * 0.8;
var timeTransition = 8000;
var numberConditions = 5; // determine based on user data
var start = 90; // data min
var end = 160; // data max

var x = d3
  .scaleLinear()
  .domain([start, end])
  .range([90, width - 200])
  .clamp(true);

var slider = d3
  .select("#bar")
  .append("svg")
  .attr("width", width - 250)
  .attr("height", 50)
  .attr("overflow", "visible");

slider
  .append("g")
  .attr("class", "slider")
  .attr("transform", "translate(" + 50 + "," + height / 2 + ")");

slider
  .append("line")
  .attr("class", "track")
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
  .attr("transform", "translate(0," + 18 + ")")
  .selectAll("text")
  .data(x.ticks(numberConditions))
  .enter()
  .append("text")
  .attr("x", x)
  .attr("text-anchor", "middle")
  .text(function(d) {
    return d;
  });

var handle = slider
  .insert("circle", ".track-overlay")
  .attr("class", "handle")
  .attr("r", 12);

/*
slider.transition() // Gratuitous intro!
    .duration(timeTransition)
    .tween("hue", function() {
      var i = d3.interpolate(start, end);
      return function(t) { hue(i(t)); };
    });
*/

var playButton = d3.select("#play-button");
/*
playButton.on("click", function() {
  var button = d3.select(this);
  if (button.text() == "Pause") {

    slider.transition() // Gratuitous intro!
        .duration(0)
        .tween("hue", function() {
          var i = d3.interpolate(x.invert(handle.attr("cx")), x.invert(handle.attr("cx")));
          return function(t) { hue(i(t)); };
        });

    button.text("Play");

  } else {

    // If button at end, start it back from the beginning
    if (x.invert(handle.attr("cx")) > end - 1) {
      starter = start
    } else {
      starter = x.invert(handle.attr("cx"))
    }

    slider.transition() // Gratuitous intro!
        .duration(timeTransition - (1000 * end / x.invert(handle.attr("cx"))))
        .tween("hue", function() {
          var i = d3.interpolate(starter, end);
          return function(t) { hue(i(t)); };
        });

    button.text("Pause");
  }

})

function hue(h) {
  handle
    .attr("cx", x(h));
  label
    .text(Math.round(h))
    .attr("class", "label")
    .attr("text-anchor", "middle")
    .attr("transform", "translate(" + x(h) + ",7)")

  node.each(function(d) {

      if (d.time < h) {
          d3.select(this)
            .select("circle")
            .style("fill", "rgba(214, 69, 65, 1)")
            .style("stroke", "black")
      } else {
        d3.select(this)
          .select("circle")
          .style("fill", "white")
          .style("stroke", "black")
      }

    })

    // If slide takes it to the end, reset button
    if (h > end - 0.01) {

      d3.select("#play-button").text("Play");
    }
  }
*/

// Check number of categories
function checkCategories(categories) {
  //change to > 1 after testing
  console.log(data)
  if (data.categories.length > 1) {
    timecourse = true;
    timecourse_fill =
      '<div id="play-button" align="center">Pause<div id="bar" align="center"></div></div>';
    document.getElementById("slider").innerHTML = timecourse_fill;

    //buildSlider();
  } else {
    timecourse = false;
  }
  return timecourse;
}
