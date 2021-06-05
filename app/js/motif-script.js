/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019-2021  Youjia Zhou, Jordan A. Berg
  zhou325 <at> sci <dot> utah <dot> edu
  jordan <dot> berg <at> biochem <dot> utah <dot> edu

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
var fs = require('fs');
var path = require("path");
var app = require('electron').remote.app;

var userDataPath = app.getPath('userData');
var session_file = userDataPath + path.sep + "session_data.json";
let database_url = JSON.parse(
  fs.readFileSync(session_file).toString())["database_url"];

showMotifs = function(_callback) {

  d3.json(database_url).then(data => {
    let metaGraph = new MetaGraph(data);
  });
  return _callback;
}

window.addEventListener("load", function(event) {

  set_tooltips();
  showMotifs()
})

function set_tooltips() {

  let app_path = window.location.pathname.split('/').slice(0, -2).join('/');

  var div = d3.select("body").append("div")
  .attr("class", "tooltip")
  .style("opacity", 0);
  
  d3.select("div#motif1.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px")
        .style("transition-delay", "1s");
      div
        .html(
          'Compare the average values between reactants and products that pass the provided threshold.\n' +
          '<img class="motif-example" src="' + app_path + '/data/examples/average_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif2.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "140px")
        .style("transition-delay", "1s");
      div
        .html(
          'Find instances of sustained perturbation along a reaction.\n' + 
          '<img class="motif-example" src="' + app_path + '/data/examples/sustained_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif3.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px")
        .style("transition-delay", "1s");
      div
        .html(
          'Find instances with one regulated modifier and one core component in the reaction.' +
          '<img class="motif-example bump-up" src="' + app_path + '/data/examples/modreg_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif4.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "200px")
        .style("transition-delay", "1s");
      div
        .html(
          'Find instances where a component is the same for input and output, is regulated, along with a modifier being regulated.' +
          '<img class="motif-example bump-up" src="' + app_path + '/data/examples/transreg_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif5.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "145px")
        .style("transition-delay", "1s");
      div
        .html(
          'Identify two neighboring enzymes with activity.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/enzyme_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif6.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "120px")
        .style("transition-delay", "1s");
      div
        .html(
          'Find instances of several neighboring metabolites with continued activity.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/metabolite_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif7_1.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "180px")
        .style("transition-delay", "1s");
      div
        .html(
          'Compare the absolute minimum values between reactants and products that pass the provided threshold.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/maxmax_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif7_2.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "180px")
        .style("transition-delay", "1s");
      div
        .html(
          'Compare the absolute maximum values between reactants and products that pass the provided threshold.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/minmin_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif8_1.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "190px")
        .style("transition-delay", "1s");
      div
        .html(
          'Compare the absolute maximum and minimum values between reactants and products, respectively, that pass the provided threshold.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/maxmin_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

  d3.select("div#motif8_2.motif_button")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px")
        .style("transition-delay", "1s");
      div
        .html(
          'Compare the absolute minimum and maximum values between reactants and products, respectively, that pass the provided threshold.' +
          '<img class="motif-example" src="' + app_path + '/data/examples/minmax_example.svg" />'
        )
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });


}