/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) Jordan A. Berg

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
var fs = require("fs");
var $ = require("jquery");
var d3 = require("d3");

window.addEventListener("load", function(event) {

  var data_div = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  d3.select("button#data_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "240px")
        .style("width", "250px");
      data_div.html("<b><u>Example Data Format:</u></b><br><br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;log2(fc)&emsp;&emsp;stat<br>metabolite name 1&emsp;&emsp;2.43&emsp;&emsp;&emsp;0.003737<br>metabolite name 2&emsp;&emsp;1.72&emsp;&emsp;&emsp;0.009739<br>metabolite name 3&emsp;&emsp;0.49&emsp;&emsp;&emsp;0.080173<br>metabolite name 4&emsp;&emsp;-2.43&emsp;&ensp;&emsp;0.000173<br>...&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;...&emsp;&emsp;&emsp;&emsp;...<br><br>Table should be tab-delimited. This can be created by saving the table in Microsoft Excel using the option \"Save as type\": \"Text (Tab delimited)\"<br><br>For most consistent behavior, you should only symbolize a decimal with a period (.), NOT a comma (,).")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

  d3.select("button#format_info")
  .on("mouseover", function(d) {
    data_div
      .style("opacity", 0.95)
      .style("left", (d3.event.pageX + 20) + "px")
      .style("top", (d3.event.pageY - 10) + "px")
      .style("height", "40px")
      .style("width", "250px");
    data_div.html("Open a non-Metaboverse-formatted datatable and have it prepared for usage in Metaboverse.")
  })
  .on("mouseout", function(d) {
    data_div.style("opacity", 0);
    data_div.html("")
  });

  d3.select("button#transcript_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "115px")
        .style("width", "250px");
      data_div.html("This datatype will accept any tab-delimited dataset with mappings to gene identifiers. Currently, there should only be an index column (no title) and a column named fold change with the calculated fold changes between experimental conditions. Please see documentation for more information on correct formatting of this table.")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

  d3.select("button#protein_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "115px")
        .style("width", "250px");
      data_div.html("This datatype will accept any tab-delimited dataset with mappings to protein identifiers. Currently, there should only be an index column (no title) and a column named fold change with the calculated fold changes between experimental conditions. Please see documentation for more information on correct formatting of this table.")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

  d3.select("button#metabolite_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "115px")
        .style("width", "250px");
      data_div.html("This datatype will accept any tab-delimited dataset with mappings to metabolite identifiers. Currently, there should only be an index column (no title) and a column named fold change with the calculated fold changes between experimental conditions. Please see documentation for more information on correct formatting of this table.")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

  d3.select("button#experiment_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "60px")
        .style("width", "250px");
      data_div.html("Enter a name for your experiment. For best stability, we recommend using underscores (\"_\") instead of spaces in the experiment name.")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

  d3.select("button#experiment_type_info")
    .on("mouseover", function(d) {
      data_div
        .style("opacity", 0.95)
        .style("left", (d3.event.pageX + 20) + "px")
        .style("top", (d3.event.pageY - 10) + "px")
        .style("height", "45px")
        .style("width", "250px");
      data_div.html("Select the experiment type you would like to analyze. Please see the documentation for more information on selecting this parameter.")
    })
    .on("mouseout", function(d) {
      data_div.style("opacity", 0);
      data_div.html("")
    });

    d3.select("button#broadcast_gene_info")
      .on("mouseover", function(d) {
        data_div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", "30px")
          .style("width", "250px");
        data_div.html("Check to broadcast gene expression values to proteins when protein values are not available.")
      })
      .on("mouseout", function(d) {
        data_div.style("opacity", 0);
        data_div.html("")
      });

    d3.select("button#broadcast_metabolite_info")
      .on("mouseover", function(d) {
        data_div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", "30px")
          .style("width", "250px");
        data_div.html("Check to broadcast metabolite values to proteins complexes.")
      })
      .on("mouseout", function(d) {
        data_div.style("opacity", 0);
        data_div.html("")
      });

    d3.select("button#modifier_collapse_info")
      .on("mouseover", function(d) {
        data_div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", "60px")
          .style("width", "250px");
        data_div.html("Check to include modifiers in reaction collapsing. Catalysts are included as outputs, inhibitors are included as inputs. Please refer to documentation for more information.")
      })
      .on("mouseout", function(d) {
        data_div.style("opacity", 0);
        data_div.html("")
      });

    d3.select("button#blocklist_info")
      .on("mouseover", function(d) {
        data_div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", "100px")
          .style("width", "250px");
        data_div.html("Enter the names for metabolites, etc., to ignore in the analysis and visualization. Names should be separated by a comma. If a spelling of the entity does not seem to be working, try searching for the correct formatting in the total entity list in the Visualization page or in the Reactome Knowledgebase.")
      })
      .on("mouseout", function(d) {
        data_div.style("opacity", 0);
        data_div.html("")
      });

    d3.select("button#collapse_threshold_info")
      .on("mouseover", function(d) {
        data_div
          .style("opacity", 0.95)
          .style("left", (d3.event.pageX + 20) + "px")
          .style("top", (d3.event.pageY - 10) + "px")
          .style("height", "55px")
          .style("width", "250px");
        data_div.html("Change the percentage of matching entites between two reactions to be considered related enough to collapse into a single collapsed reaction.")
      })
      .on("mouseout", function(d) {
        data_div.style("opacity", 0);
        data_div.html("")
      });
})

// Hyperlinks listener
const tableBrowserSettings = 'top=500,left=200,frame=true,nodeIntegration=yes,enableRemoteModule=yes,worldSafeExecuteJavaScript=yes,contextIsolation=no';

window.addEventListener("load", function(event) {
    event.preventDefault();
    event.stopPropagation();

    var user_path = window.location.pathname;
    var page = user_path.split('/').pop();
  
    document.getElementById("format-dropDatabase").onclick = function(event) {
      event.preventDefault();
      event.stopPropagation();
  
      window.open(
        'datatable.html',
        '_blank',
        tableBrowserSettings
      )
    }
})

window.addEventListener("load", function(event) {
  event.preventDefault();
  event.stopPropagation();

  // If user goes back to this page, force re-curation
  update_session_info("processed", false);

  // Format page
  var formatURL = "";
  var defaultFormat = "No file selected";
  if (formatURL !== null) {
    defaultFormat = formatURL;
  }
  $('#selectedFormat').append('<font size="2">' + defaultFormat + '</font>');

  // Transcriptomics
  var transcriptomicsURL = get_session_info("transcriptomics");
  var defaultTranscriptomics = "No file selected";
  if (transcriptomicsURL !== null) {
    defaultTranscriptomics = transcriptomicsURL;
  }
  $('#selectedTranscriptomics').append('<font size="2">' + defaultTranscriptomics + '</font>');

  // Proteomics
  var proteomicsURL = get_session_info("proteomics");
  var defaultProteomics = "No file selected";
  if (proteomicsURL !== null) {
    defaultProteomics = proteomicsURL;
  }
  $('#selectedProteomics').append('<font size="2">' + defaultProteomics + '</font>');

  // Metabolomics
  var metabolomicsURL = get_session_info("metabolomics");
  var defaultMetabolomics = "No file selected";
  if (metabolomicsURL !== null) {
    defaultMetabolomics = metabolomicsURL;
  }
  $('#selectedMetabolomics').append('<font size="2">' + defaultMetabolomics + '</font>');

  document.getElementById("menurefresh").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();

    refresh_session()
  }

  console.log("Collapse with modifiers:", get_session_info("collapseWithModifiers"))
  if (get_session_info("collapseWithModifiers") === true) {
    document.getElementById("use_modifiers_in_collapse").checked = true;
  } else {
    document.getElementById("use_modifiers_in_collapse").checked = false;
  }
  document.getElementById("use_modifiers_in_collapse").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("collapseWithModifiers") === false) {
      update_session_info("collapseWithModifiers", true);
    } else {
      update_session_info("collapseWithModifiers", false);
    }
    console.log("Reaction collapse evaluation with modifiers: ", get_session_info("collapseWithModifiers"))
  }

  console.log("Broadcast gene expression:", get_session_info("broadcastGeneExpression"))
  if (get_session_info("broadcastGeneExpression") === true) {
    document.getElementById("broadcast_gene_expression").checked = true;
  } else {
    document.getElementById("broadcast_gene_expression").checked = false;
  }
  document.getElementById("broadcast_gene_expression").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("broadcastGeneExpression") === false) {
      update_session_info("broadcastGeneExpression", true);
    } else {
      update_session_info("broadcastGeneExpression", false);
    }
    console.log("Broadcast gene expression values: ", get_session_info("broadcastGeneExpression"))
  }

  console.log("Broadcast metabolites:", get_session_info("broadcastMetabolites"))
  if (get_session_info("broadcastMetabolites") === true) {
    document.getElementById("broadcast_metabolite_expression").checked = true;
  } else {
    document.getElementById("broadcast_metabolite_expression").checked = false;
  }
  document.getElementById("broadcast_metabolite_expression").onclick = function(event) {
    event.stopPropagation();
    if (get_session_info("broadcastMetabolites") === false) {
      update_session_info("broadcastMetabolites", true);
    } else {
      update_session_info("broadcastMetabolites", false);
    }
    console.log("Broadcast metabolite values: ", get_session_info("broadcastMetabolites"))
  }

  document.getElementById("transcriptomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('transcriptomics-input').click();
  }

  document.getElementById("transcriptomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document
      .getElementById("transcriptomics-input")
      .value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedTranscriptomics').html('<font size="2">' + f.path + '</font>');
        update_session_info("transcriptomics", f.path);

        transcriptomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  document.getElementById("proteomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('proteomics-input').click();
  }

  document.getElementById("proteomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("proteomics-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedProteomics').html('<font size="2">' + f.path + '</font>');
        update_session_info("proteomics", f.path);

        proteomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };

  document.getElementById("metabolomics-dropDatabase").onclick = function(event) {
    event.preventDefault();
    event.stopPropagation();
    document.getElementById('metabolomics-input').click();
  }

  document.getElementById("metabolomics-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("metabolomics-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedMetabolomics').html('<font size="2">' + f.path + '</font>');

        update_session_info("metabolomics", f.path);

        metabolomics = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  };
  /*
  document.getElementById("reactions-input").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("reactions-input").value.split(".");

    if (
      (inputVal[inputVal.length - 1] !== "tsv") &
      (inputVal[inputVal.length - 1] !== "txt")
    ) {
      alert(
        "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
      );
    } else {
      try {
        f = event.srcElement.files[0];
        console.log("The file you dragged: ", f);
        $('#selectedReactions').html('<font size="2">' + f.path + '</font>');

        update_session_info("additional_reactions", f.path);

        reactions = true;
      } catch (error) {
        console.log(error);
        alert(
          "Input is not a .txt or .tsv file. You must upload the correct file type for the analyses to work."
        );
      }
    }
  }
  */

  if (get_session_info("experiment_name") !== null) {
    $('#updateExperimentName').val(get_session_info("experiment_name"));
  }

  if (get_session_info("experiment_type") !== null) {
    $('#updateExperiment').val(get_session_info("experiment_type"));

    if ((get_session_info("experiment_type") === "timecourse") | (get_session_info("experiment_type") === "multiple_conditions")) {
      $("#nameField").html(
        "<form>" +
        "Sample labels: " +
        "<button class='info' title='Enter the names for each condition or timepoint for you dataset in the order that they appear in the data table. Labels should be separated by a comma.'><i>i</i></button>" +
        "<br />" +
        "<br />" +
        "<input type='text' class='experimentName' id='updateExperimentLabels'></input>" +
        "</form>" +
        "<br />"
      );
      $('#updateExperimentLabels').val(get_session_info("labels"));

      document.getElementById("updateExperimentLabels").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();

        var inputVal = document.getElementById("updateExperimentLabels").value;

        try {
          console.log("Your provided labels: ", inputVal);

          update_session_info("labels", inputVal);
        } catch (error) {
          console.log(error);
          alert(
            "Labels are not valid."
          );
        }
      }

    } else {
      $("#nameField").html('');
      update_session_info("labels", "0");
    }
  }

  document.getElementById("updateExperiment").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var experiment_type = document.getElementById("updateExperiment").value;
    if (experiment_type === "null" || experiment_type === null) {
      experiment_type = null;
    } else {}

    try {
      update_session_info("experiment_type", experiment_type);
    } catch (error) {
      console.log(error);
      alert(error);
    }

    // If timecourse or multiple conditions, have user input labels in correct order as listed  in dataframe
    if ((experiment_type === "timecourse") | (experiment_type === "multiple_conditions")) {
      $("#nameField").html(
        "<form>" +
        "Sample labels: " +
        "<button class='info' title='Enter the names for each condition or timepoint for you dataset in the order that they appear in the data table. Labels should be separated by a comma.'><i>i</i></button>" +
        "<br />" +
        "<br />" +
        "<input type='text' class='experimentName' id='updateExperimentLabels'></input>" +
        "</form>" +
        "<br />"
      );

      document.getElementById("updateExperimentLabels").onchange = function(event) {
        event.preventDefault();
        event.stopPropagation();

        var inputVal = document.getElementById("updateExperimentLabels").value;

        try {
          console.log("Your provided labels: ", inputVal);

          update_session_info("labels", inputVal);
        } catch (error) {
          console.log(error);
          alert(
            "Labels are not valid."
          );
        }
      }

    } else {
      $("#nameField").html('');
      update_session_info("labels", "0");
    }
  };

  document.getElementById("updateExperimentName").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateExperimentName").value
    //.split(".");

    try {
      var replacer = String(String.fromCharCode(92) + " ");
      var inputVal = inputVal.replace(/ /g, replacer)
      console.log("Your provided experiment name: ", inputVal);
      update_session_info("experiment_name", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "Experiment name is not valid."
      );
    }
  }

  $('#updateBlocklist').val(get_session_info("blocklist"));

  document.getElementById("updateBlocklist").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var inputVal = document.getElementById("updateBlocklist").value;

    try {
      console.log("Your provided blocklisted entities: ", inputVal);

      update_session_info("blocklist", inputVal);
    } catch (error) {
      console.log(error);
      alert(
        "IDs are not valid."
      );
    }
  }

  $('#collapse_button').val(get_session_info("collapse_threshold"));

  document.getElementById("collapse_button").onchange = function(event) {
    event.preventDefault();
    event.stopPropagation();

    var collapseVal = document.getElementById("collapse_button").value;

    try {
      console.log("Collapse threshold: ", collapseVal);

      update_session_info("collapse_threshold", collapseVal);
    } catch (error) {
      console.log(error);
      alert(
        "Collapse threshold is not valid."
      );
    }
  }
});
