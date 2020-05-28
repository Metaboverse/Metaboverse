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
var plotly = require("'plotly.js-dist'");
var fs = require("fs");
var saveSVG = require("save-svg-as-png");

/*
var trace1 = {
  x: [1, 2, 3],
  y: [4, 5, 6],
  type: 'scatter'
};

var trace2 = {
  x: [20, 30, 40],
  y: [50, 60, 70],
  xaxis: 'x2',
  yaxis: 'y2',
  type: 'scatter'
};

var data = [trace1, trace2];

var layout = {
  grid: {rows: 1, columns: 2, pattern: 'independent'},
};

Plotly.newPlot('myDiv', data, layout);
*/

// For each omics data-type provided, plot timecourse trends (if only one
// timepoint or condition)
function volcano(
    rna_data,
    protein_data,
    metabolite_data) {

    /*
    var trace1 = {
      x: [1, 2, 3, 4, 5],
      y: [1, 6, 3, 6, 1],
      mode: 'markers',
      type: 'scatter',
      name: 'Team A',
      text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
      marker: { size: 12 }
    };

    var trace2 = {
      x: [1.5, 2.5, 3.5, 4.5, 5.5],
      y: [4, 1, 7, 1, 4],
      mode: 'markers',
      type: 'scatter',
      name: 'Team B',
      text: ['B-a', 'B-b', 'B-c', 'B-d', 'B-e'],
      marker: { size: 12 }
    };

    var data = [ trace1, trace2 ];

    var layout = {
      xaxis: {
        range: [ 0.75, 5.25 ]
      },
      yaxis: {
        range: [0, 8]
      },
      title:'Data Labels Hover'
    };

    Plotly.newPlot('myDiv', data, layout);

    // Allow user's to select list of elements to highlight
    // Allow threshold viz
    */


}

// For each omics data-type provided, plot timecourse trends (if more than one
// timepoint). Will only work for timecourse data
function lineTimecourse(
    rna_data,
    protein_data,
    metabolite_data) {

    /*
    var trace1 = {
      x: [1, 2, 3, 4],
      y: [10, 15, 13, 17],
      type: 'scatter'
    };

    var trace2 = {
      x: [1, 2, 3, 4],
      y: [16, 5, 11, 9],
      type: 'scatter'
    };

    var data = [trace1, trace2];

    Plotly.newPlot('myDiv', data);

    //Add label on hover
    // Allow user's to select list of elements to highlight
    */


}

// For each omics data-type provided, plot the statistical distribution
function pvalDist(
    rna_stats,
    protein_stats,
    metabolite_stats) {

    /*
    var data = [
      {
        x: ['giraffes', 'orangutans', 'monkeys'],
        y: [20, 14, 23],
        type: 'bar'
      }
    ];

    Plotly.newPlot('myDiv', data);
    */


}



// Add note that this is only considering values integrated into the global network
