/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) Youjia Zhou, Jordan A. Berg

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
    
    // BEGIN: Initialize stat threshold button and functions
    var pathway_dict = make_pathway_dictionary(
      data,
      'pathway_dictionary');
    var collapsed_pathway_dict = make_pathway_dictionary(
      data,
      'collapsed_pathway_dictionary');
    d3.select("#stat_button").on("change", function() {
      stat_input(data, collapsed_pathway_dict, pathway_dict)
    });
    var div = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0);
    d3.select("button#stat_info")
    .on("mouseover", function(d) {
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "200px")
        .style("height", "210px");
      div
        .html(`
          Provide a value to threshold statistical value where node borders are bolded.
          <br><br>
          - <b>Statistical values</b>: Node borders will be bolded if node's statistical value is less than the specified threshold.
          <br><br>
          - <b>Confidence Intervals</b>: Node borders will be bolded if the selected confidence intervals ranges between samples do not overlap.
          `)
    })
    .on("mouseout", function(d) {
      div.style("opacity", 0);
      div.html("")
    });

    var non_reaction_nodes = new Set();
    for (let n in data.nodes) {
      if (data.nodes[n].type !== "reaction" && data.nodes[n].type !== "collapsed") {
        non_reaction_nodes.add(data.nodes[n].name);
      }
    }
    make_menu(Object.assign(...Array.from(non_reaction_nodes, v => ({[v]:''}))), "pathwayMenu-motif", "Filter pattern results by an entity...");
    // END: Initialize stat threshold button and functions
  });
  return _callback;
}

window.addEventListener("load", function(event) {
  set_tooltips();
  showMotifs();
})

function set_tooltips() {

  let app_path = window.location.pathname.split('/').slice(0, -2).join('/');

  var div = d3.select("body").append("div")
  .attr("class", "tooltip")
  .style("opacity", 0);
  
  d3.select("div#motif1.motif_button")
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "140px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "200px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "145px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "120px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "180px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "180px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "190px");
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
    .on("contextmenu", function(d) {
      d3.event.preventDefault();
      div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("width", "220px")
        .style("height", "160px");
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