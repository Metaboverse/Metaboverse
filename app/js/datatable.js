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

var $ = require("jquery");
var dt = require("datatables.net")();
var fs = require("fs");
var path = require("path");
var d3 = require("d3");
const cephes = require('cephes'); // Browser


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
	Table should be tab-delimited (.txt; .tsv) or comma-delimited (.csv). <b><i>However, we recommend tab-delimited as some gene or metabolite names may contain commas which can complicate formatting.</i></b> 
	This can be created by saving the table in Microsoft Excel using the option \"Save as type\": \"Text (Tab delimited)\"; \"CSV (Comma delimited)\")
`

var table_string = `
	<table id="loaded-table-el" class="display" width="90%"></table>
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
		Add
	</div>
	<br>
`

var check_string = `Check`

var process_string = `Process`

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
				.style("height", "220px")
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
				//formatDatatable = true;
          	} catch (error) {
				console.log(error);
				alert(
					"Input is not a .csv, .txt, or .tsv file. You must upload the correct file type for the analyses to work."
				);
          	}

			// Read in datatable and show next options
			var editor, table;
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
	$('#define-groups').html(groups_string);
	$('input#fname').val('group' + String(counter))

	// initialize second table here with row names only 


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

				console.log(control_selection, "vs", experiment_selection)
				console.log(document.getElementById('fname').value)
				console.log(datatable.slice(1))
				console.log(datatable.slice(1))

				let name_array = [];
				let fc_array = [];
				let p_array = [];

				for (let i = 0; i < datatable.slice(1).length; i++) {
					let exp = experiment_selection.map(j => datatable.slice(1)[i][j]);
					let con = control_selection.map(j => datatable.slice(1)[i][j]);

					name_array.push(datatable.slice(1)[i][0]);
					fc_array.push(
						Math.log2(
							exp.reduce((a, b) => a + b, 0) / con.reduce((a, b) => a + b, 0)
						)
					);
					fc_array.push(
						ttest(exp, con).pValue()
					);


				}

				console.log(name_array)
				console.log(fc_array)
				console.log(p_array)


				control_selection = null;
				experiment_selection = null;
				d3.select("#text-group-control").html("")
				d3.select("#text-group-experiment").html("")
				counter = counter + 1;
				$('input#fname').val('group' + String(counter))
			}



			
		})


	//check_metabolites();
	process_selection();
	/*
	
	*/


	// After each control, experiment is selected, 
	// allow user to run another if timecourse or 
	// multicondition


}

function check_metabolites() {
	$('#check-metabolites').html(check_string);


}

function process_selection() {
	$('#process-table').html(process_string);


}

function process_table() {
	$('#process-table').html(process_string);


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

function write_table(file_path) {
    fs.writeFileSync(file_path, JSON.stringify(file_path), function(err) {
      if (err) throw err;
      console.log("Session data updated");
    });
    // Update session data 
    


};

function mean(arr) {
	let total = arr.reduce((a, b) => a + b, 0);
	return (total / arr.length);
}

function sd(arr){
	let avg = mean(arr);
	
	let squareDiffs = arr.map(function(value){
		let diff = value - avg;
		let sqrDiff = diff * diff;
	  	return sqrDiff;
	});
	
	let avgSquareDiff = mean(squareDiffs);
  
	let std = Math.sqrt(avgSquareDiff);
	return std;
}

function sem(std, arr) {
	// Standard error of the mean = s / √n
	return (std / Math.sqrt(arr.length))
}

function ttest(arr1, arr2, alpha=0.05) {
	// function for calculating the t-test for two independent samples

	// calculate means
	let mean1 = mean(arr1);
	let mean2 = mean(arr2);
	
	// calculate standard deviations
	let std1 = sd(arr1);
	let std2 = sd(arr2);

	// calculate standard errors
	let se1 = sem(std1, arr1);
	let se2 = sem(std2, arr2);
	
	// standard error on the difference between the samples
	sed = Math.sqrt(Math.pow(se1, 2) + Math.pow(se2, 2))
	
	// calculate the t statistic
	t_stat = (mean1 - mean2) / sed
	
	// degrees of freedom
	df = arr1.length + arr2.length - 2

	// calculate the p-value
	p = (1.0 - t.cdf(Math.abs(t_stat), df)) * 2.0
	
	return p;



}