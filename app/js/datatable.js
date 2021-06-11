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

var { ipcRenderer, ipcMain, remote } = require("electron");
var { dialog } = require("electron").remote;
var $ = require("jquery");
var dt = require("datatables.net")();
var fs = require("fs");
var path = require("path");
var d3 = require("d3");

var { jStat } = require("jstat");

var app = require("electron").remote.app;
var userDataPath = app.getPath("userData");
var session_file = userDataPath + path.sep + "session_data.json";

var info_string = `
	<b><u>Example Data Format:</u></b>
	<br><br>
	&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Sample1&emsp;&emsp;Sample2&emsp;&emsp;Sample3&emsp;&emsp;...
	<br>
	metabolite name 1&emsp;&emsp;243&emsp;&emsp;&emsp;&emsp;414&emsp;&emsp;&emsp;&emsp;&ensp;14&emsp;&emsp;&emsp;&emsp;&ensp;...
	<br>
	metabolite name 2&emsp;&emsp;172&emsp;&emsp;&emsp;&emsp;23&emsp;&emsp;&emsp;&emsp;&emsp;41&emsp;&emsp;&emsp;&emsp;&ensp;...
	<br>
	metabolite name 3&emsp;&emsp;49&emsp;&emsp;&emsp;&emsp;&ensp;532&emsp;&emsp;&emsp;&emsp;&ensp;77&emsp;&emsp;&emsp;&emsp;&ensp;...
	<br>
	metabolite name 4&emsp;&emsp;23&emsp;&ensp;&emsp;&emsp;&emsp;1423&emsp;&emsp;&emsp;&emsp;86&emsp;&emsp;&emsp;&emsp;&ensp;...
	<br>
	...&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;...&emsp;&emsp;&emsp;&emsp;&emsp;...&emsp;&emsp;&emsp;&emsp;&ensp;...&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;...
	<br><br>
	The statistical procedures used below assume data are <b><i>normally distributed</i></b>. For sequencing data, which follows a negative binomial distribution, an appropriate tool, such as DESeq2 or Limma, should be used.
	<br><br>
	Table should be tab-delimited (.txt; .tsv) or comma-delimited (.csv). <b><i>However, we recommend tab-delimited as some gene or metabolite names may contain commas which can complicate formatting.</i></b> 
	This can be created by saving the table in Microsoft Excel using the option \"Save as type\": \"Text (Tab delimited)\"; \"CSV (Comma delimited)\")
`

var table_string = `
	<table id="loaded-table-el" class="display" width="90%"></table>
`

var processed_string = `
	<table id="processed-table-el" class="display" width="90%"></table>
`

var datatypes_string = `
	Select data type:
	<br>
	<div id="button-2condition">
		2-Condition
	</div>
	<div id="button-timecourse">
		Time-course
	</div>
	<div id="button-multicondition">
		Multi-condition
	</div>

`

var groups_string = `
	Select column headers and assign to the following groups:
	<br>
	<br>
	<u>Comparison Name</u>:
	<div id="text-group-name" class="bump-right-10">
		<input type="text" id="fname" name="fname">
	</div>
	<br>
	<div id="button-group-control">
		Control
	</div>
	<div id="text-group-control" class="bump-right-10"></div>
	<br>
	<div id="button-group-experiment">
		Experiment
	</div>
	<div id="text-group-experiment" class="bump-right-10"></div>
	<br>
	<div id="button-group-add">
		Add Comparison
	</div>
	<br>
`

var export_string = `
	<div id="button-group-export">
		Export Table
	</div>
`

window.addEventListener("load", function(event) {
    event.preventDefault();
    event.stopPropagation();

	var data_div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

	d3.select("button#data_info_flat")
		.on("mouseover", function(d) {
			data_div
				.style("opacity", 0.95)
				.style("left", (d3.event.pageX + 20) + "px")
				.style("top", (d3.event.pageY - 10) + "px")
				.style("height", "280px")
				.style("width", "360px");
			data_div.html(info_string)
		})
		.on("mouseout", function(d) {
			data_div.style("opacity", 0);
			data_div.html("")
		});

    // Format page
    document.getElementById("format-datatable-input").onclick = function(event) {
        event.preventDefault();
        event.stopPropagation();
        document.getElementById('format-input').click();
    }

    document.getElementById("format-input").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();
    
        var inputVal = document
          .getElementById("format-input")
          .value.split(".");
    
        if (
            (inputVal[inputVal.length - 1] !== "csv") &
            (inputVal[inputVal.length - 1] !== "tsv") &
            (inputVal[inputVal.length - 1] !== "txt")
        ) {
          alert(
            "Input is not a .csv, .txt, or .tsv file. You must upload the correct file type for the analyses to work."
          );
        } else {
          	try {
				f = event.srcElement.files[0];
				console.log("The file you dragged: ", f);
				$('#selectedFormat').html('<font size="2">' + f.path + '</font>');
				update_session_info("formatDatatable", f.path);
          	} catch (error) {
				console.log(error);
				alert(
					"Input is not a .csv, .txt, or .tsv file. You must upload the correct file type for the analyses to work."
				);
          	}

			// Read in datatable and show next options
			var table;
			var datatable = read_table(f.path);
			$('#loaded-table').html(table_string);
			let _columns = parse_columns(datatable);

			$(document).ready(function() {
	
				table = $('#loaded-table-el').DataTable( {
					data: datatable.slice(1),
					columns: _columns,
					select: {
						style:    'api',
            			selector: 'td:first-child'
					},
					scrollX: true,
					autoWidth: false,
					"aoColumnDefs": [
						{ "bSortable": true, "aTargets": [ 0 ] },
						{ "bSortable": false, "aTargets": Array.from(_columns.keys()).slice(1) }
					]
				} );

				// from https://www.gyrocode.com/articles/jquery-datatables-how-to-select-columns/
				$('#loaded-table').on('click', 'thead th:not(:first-child)', function(e) {
					if (table.column(this, { selected: true }).length) {
						table.column(this).deselect();
						$( table.column( this ).nodes() ).removeClass( 'highlight' );
					} else {
						table.column(this).select();    
						$( table.column( this ).nodes() ).addClass( 'highlight' );  
					}
				} );
				$('#loaded-table-el').on( 'click', 'tbody td', function () {
					table.cell( this ) //.edit();
				} );
				
				$('#select-datatype').html(datatypes_string);
				d3.select("#button-2condition")
					.on("click", function(d) {
						d3.select("#button-2condition")
							.style("background", "#0080ff")
						d3.select("#button-timecourse")
							.style("background", "#f1f1f1")
						d3.select("#button-multicondition")
							.style("background", "#f1f1f1")
						select_groups(datatable, table);
					})
				d3.select("#button-timecourse")
					.on("click", function(d) {
						d3.select("#button-2condition")
							.style("background", "#f1f1f1")
						d3.select("#button-timecourse")
							.style("background", "#0080ff")
						d3.select("#button-multicondition")
							.style("background", "#f1f1f1")
						select_groups(datatable, table);
					})
				d3.select("#button-multicondition")
					.on("click", function(d) {
						d3.select("#button-2condition")
							.style("background", "#f1f1f1")
						d3.select("#button-timecourse")
							.style("background", "#f1f1f1")
						d3.select("#button-multicondition")
							.style("background", "#0080ff")
						select_groups(datatable, table);
					})

				



			} );    
		} 
	};
});

function read_table(file_path) {

	var data;
	try {
		data = fs.readFileSync(file_path, 'utf8')
	} catch (err) {
		console.error(err)
		alert("Unable to load input file.")
	}

	let path_array = file_path.split(".");
	var datatable;
	if (path_array[path_array.length - 1] === "csv") {
    	datatable = d3.csvParseRows(data);
    } else {
		datatable = d3.tsvParseRows(data);
    }

	return datatable

};

function parse_columns(datatable) {

	let columns = [];

	for (let x in datatable[0]) {
		columns.push({ title: datatable[0][x] });
	}

	return columns;
};

function select_groups(datatable, table) {

	var counter = 1;
	var processed_table;
	$('#define-groups').html(groups_string);
	$('input#fname').val('group' + String(counter))

	// initialize second table here with row names only 
	var processed_columns = [{ title: " " }];
	var processed_data = datatable.slice(1).map(function(x) {
        return [ x[0] ];
    });
	
	$('#process-table').html(processed_string);
	processed_table = $('#processed-table-el').DataTable( {
		data: processed_data,
		columns: processed_columns,
		scrollX: true,
		autoWidth: false,
		"aoColumnDefs": [
			{ "bSortable": true, "aTargets": [ 0 ] },
			{ "bSortable": false, "aTargets": Array.from(processed_columns.keys()).slice(1) }
		]
	} );
	$('#save-table').html(export_string);

	var control_selection = null;
	var experiment_selection = null;
	d3.select("#button-group-control")
		.on("click", function(d) {
			let output = handle_selection(datatable, table, "#text-group-control");
			control_selection = output[0];
			table = output[1];
		})
	
	d3.select("#button-group-experiment")
		.on("click", function(d) {
			let output = handle_selection(datatable, table, "#text-group-experiment");
			experiment_selection = output[0];
			table = output[1];
		})

	d3.select('#button-group-add')
		.on("click", function(d) {
			if (control_selection !== null && experiment_selection !== null) {

				let fc_array = [];
				let p_array = [];
				for (let i = 0; i < datatable.slice(1).length; i++) {
					let exp = experiment_selection.map(j => parseFloat(datatable.slice(1)[i][j]));
					let con = control_selection.map(j => parseFloat(datatable.slice(1)[i][j]));
					fc_array.push(
						Math.log2(
							exp.reduce((a, b) => a + b, 0) / con.reduce((a, b) => a + b, 0)
						)
					);
					p_array.push(
						ttest(exp, con)
					);
				}
				let bh_array = bh_corr(p_array);

				// Add new column headers for new data 
				let this_id = document.getElementById('fname').value;
				processed_columns.push({ title: this_id + "_fc" });
				processed_columns.push({ title: this_id + "_bh" });

				for (let i = 0; i < datatable.slice(1).length; i++) {
					processed_data[i].push(fc_array[i]);
					processed_data[i].push(bh_array[i]);
				}

				$('#process-table').html(processed_string);
				processed_table = $('#processed-table-el').DataTable( {
					destroy: true,
					data: processed_data,
					columns: processed_columns,
					scrollX: true,
					autoWidth: false,
					"aoColumnDefs": [
						{ "bSortable": true, "aTargets": [ 0 ] },
						{ "bSortable": false, "aTargets": Array.from(processed_columns.keys()).slice(1) }
					]
				} );

				control_selection = null;
				experiment_selection = null;
				d3.select("#text-group-control").html("")
				d3.select("#text-group-experiment").html("")
				counter = counter + 1;
				$('input#fname').val('group' + String(counter))
			}
		})

	d3.select('#button-group-export')
		.on("click", function(d) {
			write_table(processed_columns, processed_data)
		})
	//check_metabolites();
}

function check_metabolites() {
	$('#check-metabolites').html(check_string);
}

function handle_selection(datatable, table, selector) {

	//parse selection
	let selected_columns = [];
	table.columns( { selected: true } ).every(function(d) {
		selected_columns.push(d)
	} )
	selected_columns = selected_columns.filter(i => i !== 0); // prevent index selection

	//display selection
	let formatted_names = "";
	for (let x in selected_columns) {
		formatted_names = formatted_names + datatable[0][selected_columns[x]];
		if (x != selected_columns.length - 1) {
			formatted_names = formatted_names + ", "
		}
	}
	d3.select(selector)
		.html("<b>" + formatted_names + "</b>")

	//clear selections
	$( table.cells().nodes() ).removeClass( 'highlight' );
	table.columns( { selected: true } ).deselect();
	
	return [selected_columns, table];
}

function write_table(columns, data) {

	let output_data = "";
	for (let c in columns) {
		output_data = output_data + columns[c].title + "\t"
	}
	for (let d in data) {
		output_data = output_data + "\n" 
		for (let _d in data[d]) {
			output_data = output_data + data[d][_d] + "\t"
		}
	}

    filename = dialog
      	.showSaveDialog({
			defaultPath: ".." + path.sep + ".." + path.sep,
			properties: ["createDirectory"],
			filters: [
				{ name: "tab-delimited (*.tsv; *.txt)", extensions: ["txt", "tsv"] }
			]
      	})
      	.then(result => {
			let hasExtension = /\.[^\/\\]+$/.test(result.filePath);
			if (hasExtension === false) {
				result.filePath = `${ result.filePath }.${ "txt" }`;
			}
			filename = result.filePath;
			if (filename === undefined) {
				alert("File selection unsuccessful");
				return;
			}

			fs.writeFileSync(filename, output_data, function(err) {
				if (err) throw err;
				console.log("Table exported");
			});
		})
		.catch(err => {
			console.log(err);
		});
};

function ttest(arr1, arr2) {
	// function for calculating the f-statistic for two independent sample sets

	// calculate the p-value
	let p = jStat.anovaftest(arr1, arr2);

	return p;
}

function bh_corr(arr) {

	let sorted = arr.slice().sort( function(a, b) { return a - b } );
	let ranks = arr.map( function(v) { return sorted.indexOf(v)+1 } );

	let bh_array = [];
	for (let i = 0; i < arr.length; i++) {

		let _p = arr[i];
		let _r = ranks[i];
		let adj_p = (_r / arr.length) * _p;

		if (adj_p > 1) {
			bh_array.push(1.0);
		} else {
			bh_array.push(adj_p);
		}
	}

	return bh_array;
}