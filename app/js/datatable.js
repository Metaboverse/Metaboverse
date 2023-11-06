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

var $ = require("jquery");
var dt = require("datatables.net")(window, $);
var fs = require("fs");
var { spawn } = require('child_process');
var path = require("path");
var d3 = require("d3");
var { ipcRenderer } = require('electron');



// -------------- //
// Stat functions //
// -------------- //
var jStat = require('jstat');

function mean(arr) {
    let sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length;
}

function variance(arr, mean) {
    return arr.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / (arr.length - 1);
}

function pValue(arr1, arr2) {
    let mean1 = mean(arr1);
    let mean2 = mean(arr2);
    let var1 = variance(arr1, mean1);
    let var2 = variance(arr2, mean2);
    let t = (mean1 - mean2) / Math.sqrt(var1 / arr1.length + var2 / arr2.length);
	let df = arr1.length + arr2.length - 2;
	// We will use the fact that the t-distribution is symmetrical around 0.
    // If t is positive, we need to calculate the right tail, so we do 1 - cdf(t).
    // If t is negative, cdf(t) gives us the left tail, which is what we want.
    if (t > 0) {
        return 2 * (1 - jStat.studentt.cdf(t, df));
    } else {
        return 2 * jStat.studentt.cdf(t, df);
    }
}

function benjaminiHochberg(pValues) {
    // Sort the p-values in ascending order and keep track of their original indices
    let sortedIndices = Array.from(pValues.keys()).sort((a, b) => pValues[a] - pValues[b]);
    let sortedPValues = sortedIndices.map(i => pValues[i]);

    // Apply the Benjamini-Hochberg procedure
    let bhValues = sortedPValues.map((p, i, arr) => Math.min(1, p * arr.length / (i + 1)));

    // Put the corrected p-values back in their original order
    let correctedPValues = Array(pValues.length);
    sortedIndices.forEach((originalIndex, i) => {
        correctedPValues[originalIndex] = bhValues[i];
    });

    return correctedPValues;
}

function standardError(arr) {
    return jStat.stdev(arr, true) / Math.sqrt(arr.length);  // use the unbiased standard deviation
}

function confidenceInterval(arr, confidence) {
    const meanValue = mean(arr);
    const margin = jStat.studentt.inv((1 + confidence) / 2, arr.length - 1) * standardError(arr);
    return [meanValue - margin, meanValue + margin];
}

function confidenceInterval_twoArray(arr1, arr2){

	let confidenceLevels = [0.9, 0.95, 0.99];
	let confidenceIntervals = confidenceLevels.map(level => {
		return [level, [confidenceInterval(arr1, level), confidenceInterval(arr2, level)]];
	});

	return confidenceIntervals;
}
// ------------------ //
// End stat functions //
// ------------------ //

// Function to show the modal
function showModalWithHMDBError() {
	// Create the modal if it doesn't exist
	if ($('#hmdbErrorModal').length === 0) {
		var $modal = $('<div id="hmdbErrorModal" class="modal"></div>');
		var $modalContent = $('<div class="modal-content"></div>');
		var $closeModal = $('<span class="closeModal">&times;</span>');
		var $message = $('<p>HMDB IDs detected. The MetaboAnalyst API has historically not worked with mapping HMDB IDs. Please use the web tool instead at: <a href="https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml" target="_blank" class="hmdbErrorLink">https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml</a>.</p>');

		$modalContent.append($closeModal);
		$modalContent.append($message);
		$modal.append($modalContent);
		$('body').append($modal);

		// Close modal event
		$('.closeModal').click(function() {
			$('#hmdbErrorModal').fadeOut();
		});

		// Click outside modal to close it
		$(window).click(function(event) {
			if ($(event.target).is('.modal')) {
				$('#hmdbErrorModal').fadeOut();
			}
		});
	}

	// Show the modal
	$('#hmdbErrorModal').fadeIn();
}

var use_stat_type;
var selected_columns = [];

var opts = { // Spinner opts from http://spin.js.org/
	lines: 10, // The number of lines to draw
	length: 19, // The length of each line
	width: 5, // The line thickness
	radius: 14, // The radius of the inner circle
	scale: 1, // Scales overall size of the spinner
	corners: 1, // Corner roundness (0..1)
	speed: 1, // Rounds per second
	rotate: 0, // The rotation offset
	animation: 'spinner-line-fade-quick', // The CSS animation name for the lines
	direction: 1, // 1: clockwise, -1: counterclockwise
	color: '#454545', // CSS color or array of colors
	fadeColor: 'transparent', // CSS color or array of colors
	top: '50%', // Top position relative to parent
	left: '50%', // Left position relative to parent
	shadow: '0 0 1px transparent', // Box-shadow for the lines
	zIndex: 2000000000, // The z-index (defaults to 2e9)
	className: 'spinner', // The CSS class to assign to the spinner
	position: 'absolute', // Element positioning
  };

var mapping_api_settings = {
	"async": true,
	"crossDomain": true,
	"url": "https://www.xialab.ca/api/mapcompounds",
	"method": "POST",
	"headers": {
	  "Content-Type": "application/json",
	  "cache-control": "no-cache"
	},
	"processData": false,
	"data": ""
  }

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
	Note: Any normalization (between groups, etc.) that needs to be performed for you dataset should have already been completed. 
	<br><br>
	The statistical procedures used below assume data are <b><i>normally distributed</i></b>. For sequencing data, which follows a negative binomial distribution, an appropriate tool, such as DESeq2 or Limma, should be used. 
	Once your data have been processed by one of these tools, or something comparible, you can export a version of the output table that is only the gene names, log<sub>2</sub>(Fold Change), and adjusted p-values.
	<br><br>
	Table should be tab-delimited (.txt; .tsv) or comma-delimited (.csv). <b><i>However, we recommend tab-delimited as some gene or metabolite names may contain commas which can complicate formatting.</i></b> 
	This can be created by saving the table in Microsoft Excel using the option \"Save as type\": \"Text (Tab delimited)\"; \"CSV (Comma delimited)\")
`

var padj_string = `
	<b><u>Statistical Values:</u></b>
	<br><br>
	<b>Adjusted P-value (normal)</b>: Sample comparisons will use a 2-group ANOVA comparison to calculate a base p-value, which assumes data are normally distributed, followed by a Benjamini-Hochberg (BH) p-value correction for multiple hypothesis 
	correction. BH correction is applied here as it is not a conservative as a Bonferroni correction procedure and is thus generally better suited for exploratory data analysis.
	<br><br>
	<b>P-value (normal)</b>: Sample comparisons will use a 2-group ANOVA comparison to calculate a base p-value, with no addition p-value correction afterwards. Assumes data are normally distributed.
	<br><br>
	<b>Confidence Intervals</b>: Confidence intervals will be calculated for each experimental/control group. This method assumes the input data array are normally distributed.
	<br><br>
	If you click this button after a sample comparison has already been performed, the change in p-value handling will only be applied to new comparisons. If you want to apply the selected procedure to other previous comparisons,
	you will need to start the data processing over. You can do this by clicking the Refresh icon on the top left of this page (<b>&#x21bb;</b>), closing this window and re-opening the data formatter, or you can refresh the page by clicking "View" -> "Reload" or CTRL + R (on Windows/Linux) or CMD + R (on MacOS).
`

var dist_info_string = `
  	This view provides a summary of the measurement distributions of the various samples in your dataset. This tool is intended to guide you in understanding the statistical tests you should consider for your dataset. If the samples appear normally distributed, you can use the Adjust P-value or P-value methods provided within this tool. If the samples do not appear normally distributed, you may want to consider an alternative tool for formatting your dataset for use within Metaboverse.
`

var group_names_string = `
	Provide a name for the current comparison in the box to the right (optional).
	<br><br>
	<b>Control</b>: Select the samples belonging to the control group (wild-type, non-treated, etc.), and click the "Control" button to add these samples to this group. These samples are used as the denominator in 
	the log<sub>2</sub>(Fold Change) calculation.
	<br><br>
	<b>Experiment</b>: Select the samples belonging to the experimental group (mutant strain, treated group, etc.), and click the "Experiment" button to add these samples to this group. These samples are used as the 
	numerator in the log<sub>2</sub>(Fold Change) calculation.
	<br><br>
	Click the <b>Add Comparison</b> button once samples have been selected for the current comparison and log<sub>2</sub>(Fold Change) values and statistical values will be calculated and formatted for Metaboverse. 
	The processed data will be shown in the table below.
	<br><br>
	For experiments where you have multiple comparisons you want to process (i.e., time-course, multi-condition data), you should continue to add additional comparison groups, which will be appended to the output table.
`

var table_string = `
	<table id="loaded-table-el" class="display" width="90%"></table>
`

var processed_string = `
	<b><u>Output table preview:</u></b>
	<br>
	<br>
	<table id="processed-table-el" class="display" width="90%"></table>
`

var datatypes_string = `
	<br>
	<b><u>2) Select data type:</u></b>
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

var dist_string = `
	<button id="dist_info_flat" class="info_enhanced">
		<i>i</i>
	</button>
	&emsp;
	<b>Data distributions:</b>
`

var groups_string = `
	<br>
	<b><u>3) Select column headers from the table above and assign to the following groups:</u></b>
	<br>
	<br>

	<button id="group_name_info_flat" class="info_enhanced">
		<i>i</i>
	</button>
	&emsp;
	Comparison Name:
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
	<div id="button-group-check" title="Check provided entity names against the MetaboAnalyst database. Processing time depends on connection to MetaboAnalyst and the number of names to check.">
		Check Names
	</div>
	&emsp;
	<div id="button-group-export">
		Export Table
	</div>
`

function test_metaboanalyst_connection() {
	var settings = mapping_api_settings;
	settings.data = "{\n\t\"queryList\": \"1,3-Diaminopropane;2-Ketobutyric acid;2-Hydroxybutyric acid;\",\n\t\"inputType\": \"name\"\n}"
	  
	$.ajax(settings).done(function (response) {
	console.log(response);
	});
}

test_metaboanalyst_connection()

function test_metaboanalyst_connection_hmdb() {
	var settings = mapping_api_settings;
	settings.data = "{\n\t\"queryList\": \"HMDB0001294;HMDB0245405;HMDB0302754;\",\n\t\"inputType\": \"hmdb\"\n}"
	  
	$.ajax(settings).done(function (response) {
	console.log(response);
	});
}

test_metaboanalyst_connection_hmdb()


window.addEventListener("load", function(event) {
    event.preventDefault();
    event.stopPropagation();

	var data_div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

	document.getElementById("menurefresh-format").onclick = function(event) {
		event.preventDefault();
		event.stopPropagation();
		clear_elements()
	  }

	d3.select("button#data_info_flat")
		.on("mouseover", function(d) {
			data_div
				.style("opacity", 0.95)
				.style("left", (d3.event.pageX + 20) + "px")
				.style("top", (d3.event.pageY - 10) + "px")
				.style("height", "340px")
				.style("width", "400px");
			data_div.html(info_string)
		})
		.on("mouseout", function(d) {
			data_div.style("opacity", 0);
			data_div.html("")
		});

	d3.select("button#padj_info_flat")
		.on("mouseover", function(d) {
			data_div
				.style("opacity", 0.95)
				.style("left", (d3.event.pageX + 20) + "px")
				.style("top", (d3.event.pageY - 10) + "px")
				.style("height", "395px")
				.style("width", "360px");
			data_div.html(padj_string)
		})
		.on("mouseout", function(d) {
			data_div.style("opacity", 0);
			data_div.html("")
		});

	// Add adj_p checkbox 
	document.getElementById("stat_type").onchange = function(event) {
		use_stat_type = document.getElementById("stat_type").value;
		clear_elements();
		console.log("Using stat type: ", use_stat_type)
	}
		
    // Format page
    document.getElementById("format-datatable-input").onclick = function(event) {
        event.preventDefault();
        event.stopPropagation();
        document.getElementById('format-input').click();
    }

    document.getElementById("format-input").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();
    
		// clear all previous elements in case new data is uploaded
		clear_elements();

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
					destroy: true,
					data: datatable.slice(1),
					columns: _columns,
					select: {
						style:    'os',
            			selector: 'td:first-child',
						items: 'column'
					},
					scrollX: true,
					autoWidth: false,
					"aoColumnDefs": [
						{ "bSortable": true, "aTargets": [ 0 ] },
						{ "bSortable": false, "aTargets": Array.from(_columns.keys()).slice(1) }
					]
				} );

				// from https://www.gyrocode.com/articles/jquery-datatables-how-to-select-columns/
				/*$('#loaded-table').on('click', 'thead th:not(:first-child)', function(e) {
					if (table.column(this, { selected: true }).length) {
						table.column(this).deselect();
						$( table.column( this ).nodes() ).removeClass( 'highlight' );
					} else {
						table.column(this).select();    
						$( table.column( this ).nodes() ).addClass( 'highlight' );  
					}
				} );
				*/

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
				
				// Add distribution plots here
				$('#loaded-distributions-title').html(dist_string);
				d3.select("button#dist_info_flat")
					.on("mouseover", function(d) {
						data_div
							.style("opacity", 0.95)
							.style("left", (d3.event.pageX + 20) + "px")
							.style("top", (d3.event.pageY - 10) + "px")
							.style("height", "110px")
							.style("width", "360px");
						data_div.html(dist_info_string)
					})
					.on("mouseout", function(d) {
						data_div.style("opacity", 0);
						data_div.html("")
					});
				makeDistributions(datatable);
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

	$('#loaded-table').on('click', 'thead th:not(:first-child)', function(e) {
		var column = table.column(this);
		var index = column.index();
	
		if ($(column.nodes()).hasClass('selected')) {
			$(column.nodes()).removeClass('selected highlight');
			
			// remove the column's index from selected_columns
			var columnIndex = selected_columns.indexOf(index);
			if (columnIndex > -1) {
				selected_columns.splice(columnIndex, 1);
			}
		} else {
			$(column.nodes()).addClass('selected highlight');
	
			// add the column's index to selected_columns
			if (selected_columns.indexOf(index) === -1) {
				selected_columns.push(index);
			}
		}
	});

	$('#loaded-table-el').on( 'click', 'tbody td', function () {
		table.cell( this ) //.edit();
	} );


	var counter = 1;
	var processed_table;
	$('#define-groups').html(groups_string);
	$('input#fname').val('group' + String(counter))

	var groups_div = d3.select("body").append("div")
		.attr("class", "tooltip")
		.style("opacity", 0);
	d3.select("button#group_name_info_flat")
		.on("mouseover", function(d) {
			groups_div
				.style("opacity", 0.95)
				.style("left", (d3.event.pageX + 20) + "px")
				.style("top", (d3.event.pageY - 10) + "px")
				.style("height", "330px")
				.style("width", "360px");
				groups_div.html(group_names_string)
		})
		.on("mouseout", function(d) {
			groups_div.style("opacity", 0);
			groups_div.html("")
		});
	
	// initialize second table here with row names only 
	var processed_columns = [{ title: " " }];
	var processed_data = datatable.slice(1).map(function(x) {
        return [ x[0] ];
    });
	//update table
	processed_table = update_table(
		'#process-table', 
		'#processed-table-el', 
		processed_string, 
		processed_data, 
		processed_columns
	);
	$('#save-table').html(export_string);

	var control_selection = null;
	var experiment_selection = null;

	d3.select("#button-group-control")
		.on("click", function(d) {
			let output = handle_selection(datatable, table, selected_columns, "#text-group-control");
			control_selection = output[0];
			table = output[1];
			selected_columns = output[2];
		})
	
	d3.select("#button-group-experiment")
		.on("click", function(d) {
			let output = handle_selection(datatable, table, selected_columns, "#text-group-experiment");
			experiment_selection = output[0];
			table = output[1];
			selected_columns = output[2];
		})

	d3.select('#button-group-add')
		.on("click", function(d) {
			console.log(selected_columns)
			if (control_selection != null && experiment_selection != null) {

				let fc_array = [];
				let p_array = [];
				let ci_array = [];
				for (let i = 0; i < datatable.slice(1).length; i++) {
					let exp = experiment_selection.map(j => parseFloat(datatable.slice(1)[i][j]));
					let con = control_selection.map(j => parseFloat(datatable.slice(1)[i][j]));
					fc_array.push(
						Math.log2(
							exp.reduce((a, b) => a + b, 0) / con.reduce((a, b) => a + b, 0)
						)
					);
					p_array.push(pValue(exp, con));

					ci_array.push(confidenceInterval_twoArray(exp, con));
				}
				let bh_array = benjaminiHochberg(p_array);

				// Add new column headers for new data 
				let stat_type;
				let this_id = document.getElementById('fname').value;
				processed_columns.push({ title: this_id + "_fc" });
	
				if (use_stat_type === "Adjusted P-value (normal)") {
					stat_type = "adj-p"
				} else if (use_stat_type === "P-value (normal)") {
					stat_type = "p"
				} else if (use_stat_type === "Confidence Intervals") {
					stat_type = "ci"
				} else {
					stat_type = "stat"
				}
				processed_columns.push({ title: this_id + "_" + stat_type });

				for (let i = 0; i < datatable.slice(1).length; i++) {
					processed_data[i].push(fc_array[i]);
					if (use_stat_type === "Adjusted P-value (normal)") {
						processed_data[i].push(bh_array[i]);
					} else if (use_stat_type === "P-value (normal)") {
						processed_data[i].push(p_array[i]);
					} else if (use_stat_type === "Confidence Intervals") {
						processed_data[i].push(ci_array[i]);
					} else {
						processed_data[i].push(p_array[i]);
					}
				}
				console.log(processed_data)
				//update table
				processed_table = update_table(
					'#process-table', 
					'#processed-table-el', 
					processed_string, 
					processed_data, 
					processed_columns
				);

				control_selection = null;
				experiment_selection = null;
				d3.select("#text-group-control").html("")
				d3.select("#text-group-experiment").html("")
				counter = counter + 1;
				$('input#fname').val('group' + String(counter))
			}
		})

	//check_metabolites();
	d3.select('#button-group-check')
		.on("click", function(d) {

			// Print error message that this feature is disabled for the time being 
			//alert("This feature is currently disabled as the MetaboAnalyst API is no longer available. We are actively working on a solution. Please check back in a later version.")
			// Skip the rest of the function 
			//return

			//get current names
			var entity_string = "";
			var entity_names = datatable.slice(1).map(function(x) {
				entity_string = entity_string + x[0] + ";"
			});

			// inject in data string
			var settings = mapping_api_settings;
			var dataObject = {
				queryList: entity_string,
			  };

			if (entity_string && entity_string.toLowerCase().includes('hmdb')) {
				showModalWithHMDBError();
				dataObject.inputType = "hmdb";  // Add the inputType property
			} else {
				dataObject.inputType = "name";
				settings.data = JSON.stringify(dataObject, null, '\t');
				console.log(settings);

				// run ajax query and update table with names ("Match")
				d3.select("#button-group-check")
					.html("Checking...");
				var target = document.getElementById('processed-table-el')
				var spinner = new Spinner(opts);
				if (d3.event.button === 0) { spinner.spin(target) };
				$.ajax(settings).done(function (response) {
					// for each Match, if not NA, use; else ignore and use Query
					// update row names
					for (let i in processed_data) {
						if (response.Match[i] !== "NA") {
							processed_data[i][0] = response.Match[i];
						}
					}
					console.log('NameMapperResponse: ', response)

					//update table
					processed_table = update_table(
						'#process-table', 
						'#processed-table-el', 
						processed_string, 
						processed_data, 
						processed_columns
					);
					spinner.stop();
					d3.select("#button-group-check")
						.html("Check Names");
				});
			}
		})

	d3.select('#button-group-export')
		.on("click", function(d) {
			saveTableData(processed_columns, processed_data)
		})
}

function update_table(_identifier, _element, _display, _data, _cols) {

	console.log(_cols)

	$(_identifier).html(_display);
	let _processed_table = $(_element).DataTable( {
		destroy: true,
		data: _data,
		columns: _cols,
		scrollX: true,
		autoWidth: false,
		"aoColumnDefs": [
			{ "bSortable": true, "aTargets": [ 0 ] },
			{ "bSortable": false, "aTargets": Array.from(_cols.keys()).slice(1) }
		]
	} );
	return _processed_table;
}


function handle_selection(datatable, table, selected_columns, selector) {

	console.log(selector)
	console.log(selected_columns)
	this_selection = selected_columns.filter(i => i !== 0); // prevent index selection

	//display selection
	let formatted_names = "";
	for (let x in this_selection) {
		formatted_names = formatted_names + datatable[0][this_selection[x]];
		if (x != this_selection.length - 1) {
			formatted_names = formatted_names + ", "
		}
	}
	d3.select(selector)
		.html("<b>" + formatted_names + "</b>")

	//clear selections
	$(table.cells().nodes()).removeClass('selected highlight');
	selected_columns = [];

	return [this_selection, table, selected_columns];
}

async function saveTableData(processed_columns, processed_data) {

	let output_data = "";
	for (let c in processed_columns) {
		output_data = output_data + processed_columns[c].title + "\t"
	}
	for (let d in processed_data) {
		output_data = output_data + "\n" 
		for (let _d in processed_data[d]) {
			let write_data = processed_data[d][_d];
			if (Array.isArray(write_data)) {
				output_data = output_data + JSON.stringify(write_data) + "\t";
			} else {
				output_data = output_data + write_data + "\t";
			}
		}
	}

	try {
	  const filename = await ipcRenderer.invoke('save-file-dialog-table', output_data);
	  if (filename === undefined) {
		alert("File selection unsuccessful");
	  }
	} catch (err) {
	  console.log(err);
	}
  }

function clear_elements() {
	$('#define-groups').html("");
	$('#process-table').html("");
	$('#save-table').html("");
	d3.select("#button-2condition")
		.style("background", "#f1f1f1")
	d3.select("#button-timecourse")
		.style("background", "#f1f1f1")
	d3.select("#button-multicondition")
		.style("background", "#f1f1f1")
	
}

function makeDistributions(datatable) {
	// See https://d3-graph-gallery.com/graph/density_basic.html
	let _d_table = datatable.slice(1);
	let sample_n = _d_table[0].length;
	let _d_arrays = [];
	
	for (let i = 1; i < sample_n; i++) {
		_this_sample = [];
		for (let j = 0; j < _d_table.length; j++) {
			_this_sample.push(parseFloat(_d_table[j][i]));
		}
		_d_arrays.push(_this_sample);
	}

	var margin = {top: 30, right: 200, bottom: 30, left: 50},
		width = window.innerWidth - margin.left - margin.right,
		height = 250;

	d3.select("svg").remove();
	var svg = d3.select("#loaded-distributions")
		.append("svg")
			.attr("width", width + margin.left + margin.right)
			.attr("height", height + margin.top + margin.bottom)
		.append("g")
			.attr("transform",
				"translate(" + margin.left + "," + margin.top + ")");
	// add the x Axis
	// See https://stackoverflow.com/a/39343864
	function getMin(a){
		return Math.min(...a.map(e => Array.isArray(e) ? getMin(e) : e));
	}
	function getMax(a){
		return Math.max(...a.map(e => Array.isArray(e) ? getMax(e) : e));
	}

	var x = d3.scaleLinear()
		.domain([getMin(_d_arrays), getMax(_d_arrays)])
		.range([0, width]);
	svg.append("g")
		.attr("transform", "translate(0," + height + ")")
		.call(d3.axisBottom(x));

	// add the y Axis
	var y = d3.scaleLinear()
		.range([height, 0])
		.domain([0, 0.12]);
	svg.append("g")
		.call(d3.axisLeft(y));

	// Compute kernel density estimation
	var kde = kernelDensityEstimator(kernelEpanechnikov(7), x.ticks(40))

	for (let z in _d_arrays) {
		_this_array = _d_arrays[z];

		// Plot the area
		svg.append("path")
			.attr("class", "mypath")
			.datum(kde(_this_array))
			.attr("fill", "none")
			.attr("opacity", ".8")
			.attr("stroke", "#000")
			.attr("stroke-width", 1)
			.attr("stroke-linejoin", "round")
			.attr("d",  d3.line()
			.curve(d3.curveBasis)
			.x(function(d) { return x(d[0]); })
			.y(function(d) { return y(d[1]); })
		);
	}
}

// Function to compute density
function kernelDensityEstimator(kernel, X) {
	// See https://d3-graph-gallery.com/graph/density_basic.html
	return function(V) {
		return X.map(function(x) {
		return [x, d3.mean(V, function(v) { return kernel(x - v); })];
		});
	};
}
function kernelEpanechnikov(k) {
	// See https://d3-graph-gallery.com/graph/density_basic.html
	return function(v) {
		return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
	};
}

// Unit tests 
function test_pValue() {
	let exp = [1, 2, 3, 4, 5];
	let con = [5, 6, 7, 8, 9];
	let p = pValue(exp, con);
	return p.toFixed(5);
}

function test_benjaminiHochberg() {
	let p = [0.5, 0.1, 0.2, 0.3, 0.4, 0.001, 0.04];
	let bh = benjaminiHochberg(p);
	// round bh values in array to 2 decimals 
	for (let i in bh) {
		bh[i] = bh[i].toFixed(2);
	}
	return bh;
}

function test_confidenceInterval_twoArray() {
	let exp = [1, 2, 3, 4, 5];
	let con = [5, 6, 7, 8, 9];
	let ci = confidenceInterval_twoArray(exp, con);
	
	// Convert all values in arrays to 2 decimals
	for (let i in ci) {
		ci[i][1] = ci[i][1].map(arr => arr.map(x => parseFloat(x.toFixed(2))));
	}  
	return ci;
}

function arraysEqual(a, b) {
    return a.length === b.length && a.every((val, i) => val === b[i]);
}

function multiDimensionalArrayEqual(a, b) {
    return JSON.stringify(a) === JSON.stringify(b);
}

function test() {
	var assert = require('assert');
	var { JSDOM } = require('jsdom');
	var jsdom = new JSDOM('<!doctype html><html><body></body></html>');
	var { window } = jsdom;
	$ = global.jQuery = require('jquery')(window);

	let valid_p = 0.00395;
	let valid_bh = ["0.50","0.23","0.35","0.42","0.47","0.01","0.14"];
	let valid_ci = [
		[0.9, [[1.49, 4.51], [5.49, 8.51]]],
		[0.95, [[1.04, 4.96], [5.04, 8.96]]],
		[0.99, [[-0.26, 6.26], [3.74, 10.26]]]
	];
	
	assert(test_pValue() == valid_p)
	assert(arraysEqual(test_benjaminiHochberg(), valid_bh))
	assert(multiDimensionalArrayEqual(test_confidenceInterval_twoArray(), valid_ci))
	
}
module.exports = test
