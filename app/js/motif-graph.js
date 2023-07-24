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
var fs = require("fs");
var savePNG = require("save-svg-as-png");

var _width = (0.45 * window.innerWidth) + 50;
var _height = 675;

var sample = 0;
var exclude_idx = -1;
var last_click = 0;
var cov_threshold = 0.1;
var significance_weight = 3;
var nonsignificance_weight = 1;
var cofactor = "";
var selected_pattern = "";

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
  position: 'relative', // Element positioning
};

get_session_info("database_url", (err, value) => {
  if (err) {
    console.log(err)
  } else {
    var database_url = value;
  }
});

class MetaGraph {
  constructor(data) {

    console.log(data.metadata)
    var stat_type = data.metadata.stat_type;
    set_stat_button(stat_type);
    
    // Get the data
    var update_output = update_nodes_links(
      data.nodes,
      data.links
    );
    data.nodes = update_output[0];
    data.links = update_output[1];

    this.data = data;
    this.nodes = data.nodes;
    this.links = data.links;
    this.stat_type = stat_type;

    // Generate expression and stats dictionaries
    var dict_output = create_dictionaries(this.nodes);
    this.expression_dict = dict_output[0];
    this.stats_dict = dict_output[1];
    this.inferred_dict = dict_output[2];

    this.link_neighbors = create_link_neighbors(
      this.nodes,
      this.links
    );

    if (this.data.metadata.blocklist === null) {
      this.data.metadata.blocklist = "";
    }
    this.blocklist = this.data.species_blocklist;
    this.blocklist = complete_blocklist(
      this.blocklist,
      this.data.metadata.blocklist,
      this.nodes
    )
    
    this.reaction_dict = data.reaction_dictionary;
    console.log(this.reaction_dict)
    this.mod_pathway_dictionary = {};
    for (let p in data.pathway_dictionary) {
      this.mod_pathway_dictionary[data.pathway_dictionary[p].id] = data.pathway_dictionary[p];
    }
    this.pathway_dictionary = make_pathway_dictionary(
      data,
      "pathway_dictionary"
    );

    this.collapsed_reaction_dict = data.collapsed_reaction_dictionary;
    console.log(this.collapsed_reaction_dict)
    this.mod_collapsed_pathways = data.mod_collapsed_pathways;
    this.collapsed_pathway_dict = make_pathway_dictionary(
      data,
      "collapsed_pathway_dictionary"
    );
    this.path_mapper = data.motif_reaction_dictionary;
    this.degree_dict = data.degree_dictionary;
    this.categories = data.categories;
    this.labels = data.labels;
    let neighbors_output = make_neighbors_dictionary(
      data,
      this.degree_dict
    );
    this.neighbors_dictionary = neighbors_output[0];
    this.collapsed_neighbors_dictionary = neighbors_output[1];

    this.metabolite_species_dictionary = make_metabolite_species_dictionary(
      data
    );
    this.entity_species_reverse_dictionary = make_entity_species_r_dictionary(
      data
    );

    if (this.labels === null) {
      this.names = ['0']
    } else {
      this.names = this.labels.split(',');
    }

    timecourse = checkCategories(this.categories, this.names);
    if (timecourse === true) {
      populateExclusions(this.categories, this.names);
    } else {
      var select = document.getElementById("exclude_type");
      var option = document.createElement('option');
      option.text = option.value = "No exclusion";
      select.add(option, 0);
    }

    let superPaths = make_superPathway_dictionary(data);
    let superPathList = [];
    for (let k in superPaths) {
      superPathList.push(superPaths[k].pathway_id);
    }
    this.superPathwayDict = superPathList;

    // Generate stamp view output
    try {
      this.stamp_svg = d3.select("#stamp-view-svg");
      this.stamp_svg_width = parseFloat(this.stamp_svg.style("width"));
      this.stamp_svg_height = parseFloat(this.stamp_svg.style("height"));
      this.stamp_svg_margin = {
        "top": 5,
        "left": 5,
        "horizontal": 10,
        "vertical": 10
      };
      this.stamp_svg_frame_group = this.stamp_svg.append("g")
        .attr("id", "stamp-frame-group");
      this.stamp_svg_link_group = this.stamp_svg.append("g")
        .attr("id", "stamp-link-group");
      this.stamp_svg_circle_group = this.stamp_svg.append("g")
        .attr("id", "stamp-circle-group");
      this.stamp_svg_selection_group = this.stamp_svg.append("g")
        .attr("id", "stamp-selection-group");

      // Generate motif-pathway box in stamp view
      this.mp_svg = d3.select("#motif-pathway-svg");
      this.mp_svg_width = parseFloat(this.mp_svg.style("width"));
      this.mp_motif_group = this.mp_svg.append("g")
        .attr("id", "mp-motif-group");
      this.mp_motif_link_group = this.mp_svg.append("g")
        .attr("id", "mp-motif-link-group");
      this.mp_motif_circle_group = this.mp_svg.append("g")
        .attr("id", "mp-motif-circle-group");
      this.mp_pathway_group = this.mp_svg.append("g")
        .attr("id", "mp-pathway-group");
      this.mp_selection_group = this.mp_svg.append("g")
        .attr("id", "mp-selection-group");

    } catch (e) {}
    this.initLines();
    this.motifSearch();
  }

  // Add line plot space if time-course/multi-condition
  initLines() {
    if (timecourse === true) {
      $("#frame-line-plots").append(
        `
          <div class='frame frame-motif-3'>
            <div class='line-sub-frame'>
              <h6><b>Selected Reaction Motif</b></h6>
              <div id='line-plot-container'></div>
            </div>
            <i>
              <font size='2' class='line-description'>
                View how <b>measured</b> components of a given reaction change across the time-course or across conditions.
                <br>
                Dashed lines indicate an entity with only one representative measurement (i.e., a single timepoint was provided for that datapoint).
              </font>
            </i>
            <div title="Click to save the current graph view">
              <button id="saveLineGraph" class="option_button">Export PNG</button>
            </div>
            <div title="Click to save the current graph view as an SVG. Edges may only show up in Inkscape, not Adobe Illustrator.">
              <button id="saveLineSVG" class="option_button">Export SVG</button>
            </div>
          </div>
          <br>
          <br>
          <br>
        `  
      );
      d3.select("#saveLineGraph").on("click", function() {
        savePNG.saveSvgAsPng(
          d3.select("#line-graph")._groups[0][0],
          "plot.png", {
            encoderOptions: 1,
            scale: 10,
            encoderType: "image/png"
          }
        );
      });
      d3.select("#saveLineSVG").on("click", function() {

        var xmlns = "http://www.w3.org/2000/xmlns/";
        var xlinkns = "http://www.w3.org/1999/xlink";
        var svgns = "http://www.w3.org/2000/svg";
 
        var _this_svg = d3.select("#line-graph")._groups[0][0].cloneNode(true);
        var _Serializer = new XMLSerializer();
        var svg_string = _Serializer.serializeToString(_this_svg);
  
        var filename = dialog
          .showSaveDialog({
            title: "plot",
            defaultPath: ".." + path.sep + ".." + path.sep,
            properties: ["createDirectory"],
            filters: [{
              name: "svg",
              extensions: ["svg"]
            }]
          })
          .then(result => {
            let hasExtension = /\.[^\/\\]+$/.test(result.filePath);
            if (hasExtension === false) {
              result.filePath = `${ result.filePath }.${ "svg" }`;
            }
            console.log(result)

            filename = result.filePath;
            if (filename === undefined) {
              alert("File selection unsuccessful");
              return;
            }
            fs.writeFileSync(filename, svg_string, 'utf-8');
            console.log(filename);
          })
          .catch(err => {
            console.log(err);
          });
      });    
    }
  }

  watchSlider() {
    if (this.motif !== undefined) {
      if (timecourse === true) {
        d3.select("svg#slide")
          .on("click", () => {
            console.log(this)
            let sample_idx = d3.select("circle#dot").attr("x");
            this.exclude_type_dropdown = document.getElementById("exclude_type");
            exclude_idx = this.names.indexOf(this.exclude_type_dropdown.value);
            if (sample_idx !== last_click) {
              reset_objects();
              this.drawMotifSearchResult(
                this.motif, sample_idx, exclude_idx);
              last_click = sample_idx;
            }
          })
      }
    }
  }

  watchMenu() {
    if (this.motif !== undefined) {
      let filtered_motifs = [];
      d3.select("#pathwayMenu-motif")
        .on("change", () => {
          let exclude_type_dropdown = document.getElementById("exclude_type");
          let exclude_idx = this.names.indexOf(exclude_type_dropdown);
          let sample_idx = 0;

          // get filtering cofactor
          var filter_cofactor = document.getElementById("pathwayMenu-motif").value;
          if (filter_cofactor === "No metabolite co-factor selection...") {
            filtered_motifs = this.motif;
          } else {
            // get species ID for cofactor
            let cofactor_id = this.metabolite_species_dictionary[filter_cofactor];

            // If "No metabolite co-factor selection...", no selection
            filtered_motifs = [];
            for (let condition in this.motif) {
              filtered_motifs.push([]);
              for (let motif in this.motif[condition]) {
                let entities = new Set();
                for (let r in this.motif[condition][motif].reactants) {
                  entities.add(this.motif[condition][motif].reactants[r])
                }
                for (let p in this.motif[condition][motif].products) {
                  entities.add(this.motif[condition][motif].products[p])
                }
                for (let m in this.motif[condition][motif].modifiers) {
                  entities.add(this.motif[condition][motif].modifiers[m][0])
                }
                //console.log(entities)
                for (let c in cofactor_id) {
                  if (entities.has(cofactor_id[c])) {
                    filtered_motifs[condition].push(this.motif[condition][motif])
                    break;
                  }
                }
              }
            }
          }

          reset_objects();
          this.drawMotifSearchResult(
            filtered_motifs, sample_idx, exclude_idx);
          this.motif = filtered_motifs;
        }
      );
    }
  }

  watchType() {
    d3.select("#sort_type")
      .on("change", () => {
        reset_objects();
        this.sort_type_dropdown = document.getElementById("sort_type");
        let sample_idx = 0;
        try {
          let sample_idx = d3.select("circle#dot").attr("x");
        } catch(err) {}
        this.exclude_type_dropdown = document.getElementById("exclude_type");
        exclude_idx = this.names.indexOf(this.exclude_type_dropdown.value);
        this.drawMotifSearchResult(
          this.motif, sample_idx, exclude_idx);
      })
  }

  watchExclude() {
    d3.select("#exclude_type")
      .on("change", () => {
        reset_objects();
        this.exclude_type_dropdown = document.getElementById("exclude_type");
        exclude_idx = this.names.indexOf(this.exclude_type_dropdown.value);
        let sample_idx = d3.select("circle#dot").attr("x");
        this.drawMotifSearchResult(
          this.motif, sample_idx, exclude_idx);
      })
  }

  watchExport() {
    d3.select("#saveTable")
      .on("click", () => {
        // this.motif, sample_idx, exclude_idx

        // for each item (condition) in this.motif, output with column describing the condition
        var i = 0;
        var rows = [];
        let column_names = [
          "index", 
          "condition", 
          "id", 
          "name", 
          "collapsed", 
          "reactants", 
          "products", 
          "modifiers", 
          "additional_components",
          "magnitude_change",
          "source_p_value",
          "target_p_value",
          "averaged_p_value"
        ];
        rows.push(column_names);

        for (let condition in this.motif) {
          for (let motif in this.motif[condition]) {
            let this_entry = this.motif[condition][motif];

            let these_reactants = [];
            for (let r in this_entry["reactants"]) {
              if (this_entry["reactants"][r] in this.entity_species_reverse_dictionary) {
                these_reactants.push(this.entity_species_reverse_dictionary[this_entry["reactants"][r]][0]);
              } else {
                these_reactants.push(this_entry["reactants"][r]);
              }
            }
            let these_products = [];
            for (let p in this_entry["products"]) {
              if (this_entry["products"][p] in this.entity_species_reverse_dictionary) {
                these_products.push(this.entity_species_reverse_dictionary[this_entry["products"][p]][0]);
              } else {
                these_products.push(this_entry["products"][p]);
              }
            }
            let these_modifiers = [];
            for (let m in this_entry["modifiers"]) {
              if (this_entry["modifiers"][m][0] in this.entity_species_reverse_dictionary) {
                these_modifiers.push([
                  this.entity_species_reverse_dictionary[this_entry["modifiers"][m][0]][0],
                  this_entry["modifiers"][m][1]
                ]);
              } else {
                these_modifiers.push(this_entry["modifiers"][m]);
              }
            }

            let entry = [
              String(i), 
              String(condition), 
              String(this_entry["id"]), 
              String(this_entry["name"]), 
              String(this_entry["collapsed"]), 
              String(these_reactants),
              String(these_products),
              String(these_modifiers), 
              String(this_entry["additional_components"]), 
              String(this_entry["magnitude_change"]), 
              String(this_entry["p_values"]["source"]),
              String(this_entry["p_values"]["target"]),
              String(this_entry["p_values"]["agg"])
            ];
            rows.push(entry);
            i += 1;
          }
        }

        // Source: https://stackoverflow.com/a/14966131/9571488
        let csvContent = "data:text/tab-separated-values;charset=utf-8,";
        rows.forEach(function(rowArray) {
          let row = rowArray.join("\t");
          csvContent += row + "\r\n";
        });

        var encodedUri = encodeURI(csvContent);
        // End code snippet

        // Source: https://stackoverflow.com/a/14966131/9571488
        let timestamp = formatDate(new Date());
        let filename = this.data.metadata.experiment_name + "_" + selected_pattern + "_patterns_" + timestamp + ".tsv";
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", filename);
        document.body.appendChild(link); // Required for FF
        link.click(); // This will download the data file named "my_data.csv".
        // End code snippet
      })
  }

  processIdenticals(exclude=true) {

    for (let i in this.motif) {
      for (let m = 0; m < this.motif[i].length; m++) {
        let r = this.motif[i][m].reactants;
        let p = this.motif[i][m].products;

        let r_names = [];
        let p_names = [];

        for (let rr in r) {
          r_names.push(this.nodes[r[rr]].map_id);

        }
        for (let pp in p) {
          p_names.push(this.nodes[p[pp]].map_id);
        }

        if (
            (JSON.stringify(r_names.sort()) === JSON.stringify(p_names.sort()))
            &&
            (exclude === true)) {
          this.motif[i].splice(m, 1);
          m--;
        }
      }
    }
  }

  motifSearch() {

    var target = document.getElementById('stamp-view-container')
    var spinner = new Spinner(opts);

    console.log(this)

    d3.select("#motif1")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "average";
        highlight_selection("#avg_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#avg_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_Avg(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif2")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "sustained";
        highlight_selection("#sustained_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#sustained_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_Sustained(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif3")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "modreg";
        highlight_selection("#modreg_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#modreg_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = modifierReg(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif4")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "transreg";
        highlight_selection("#transreg_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#transreg_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = modifierTransport(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals(false);
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif5")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "enzyme";
        highlight_selection("#enzyme_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#enzyme_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let eval_neighbors_dictionary;
        if (eval_collapsed === false) {
          eval_neighbors_dictionary = this.neighbors_dictionary;
        } else {
          eval_neighbors_dictionary = this.collapsed_neighbors_dictionary;
        }
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = enzymeMotif(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          eval_neighbors_dictionary,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories,
          this.nodes,
          this.link_neighbors)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawTwoReactionSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif6")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "activity";
        highlight_selection("#activity_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#activity_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let eval_neighbors_dictionary;
        if (eval_collapsed === false) {
          eval_neighbors_dictionary = this.neighbors_dictionary;
        } else {
          eval_neighbors_dictionary = this.collapsed_neighbors_dictionary;
        }
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = activityMotif(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          eval_neighbors_dictionary,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories,
          this.nodes,
          this.link_neighbors)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawTwoReactionSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif7_1")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "maxmax";
        highlight_selection("#maxmax_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#maxmax_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_MaxMax(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif7_2")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "minmin";
        highlight_selection("#minmin_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#minmin_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_MinMin(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif8_1")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "maxmin";
        highlight_selection("#maxmin_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#maxmin_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_MaxMin(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
      })

    d3.select("#motif8_2")
      .on("mousedown", function() {
        if (d3.event.button === 0) { spinner.spin(target) };
      })
      .on("click", () => {
        selected_pattern = "minmax";
        highlight_selection("#minmax_num");
        reset_dot();
        reset_filter();
        reset_objects();
        let threshold = d3.select("#minmax_num").node().value;
        this.sort_type_dropdown = document.getElementById("sort_type");
        let this_reaction_dict = get_reaction_dict(this);
        this.motif = motifSearch_MinMax(
          threshold,
          this_reaction_dict,
          this.expression_dict,
          this.stats_dict,
          this.stat_type,
          stat_value,
          this.inferred_dict,
          this.link_neighbors,
          this.path_mapper,
          this.degree_dict,
          this.blocklist,
          this.categories)
        this.processIdenticals();
        this.watchSlider();
        this.watchMenu();
        this.watchType();
        this.watchExclude();
        this.drawMotifSearchResult(this.motif, 0, exclude_idx);
        this.watchExport();
        spinner.stop();
        console.log(this.motif)
      })
  }

  drawMotifSearchResult(motifs, indexer, exclusion_indexer) {

    let current_motifs = motifs[indexer];

    // Exclude any motifs in multiple time-points/conditions
    let motif_list = get_included_motifs(
      current_motifs,
      motifs,
      exclusion_indexer);

    if (motif_list.length > 1000) {
      alert(String(motif_list.length) + " reaction patterns identified. This may take a while, or you could try increasing the threshold to a more stringent level.")
    }

    // Sort motifs based on user input
    let this_stat_type = this.stat_type;
    let sort_type = "";
    if (this_stat_type === "array") {
      sort_type = "Sort Magnitude Change";
    } else {
      sort_type = this.sort_type_dropdown.value;
    }
    console.log(sort_type)
    let motif_significance = {
      'both': [],
      'one': [],
      'none': []
    };

    let sort_output = sort_motifs(
      motif_list,
      motif_significance,
      sort_type
    );
    motif_list = sort_output[0];
    motif_significance = sort_output[1];

    // Remove any duplicated items
    motif_list = remove_duplicate_motifs(motif_list);
    motif_significance = remove_duplicate_motifs(motif_significance);
    
    let stamp_height = 50;
    this.stamp_svg_height = Math.ceil(
      motif_list.length / 3) *
      (stamp_height + this.stamp_svg_margin.vertical);
    this.stamp_svg.attr("height", this.stamp_svg_height);

    let stamp_width = this.stamp_svg_width / 3 -
      this.stamp_svg_margin.horizontal;

    let clicked_stamp_reaction;
    let sg = this.stamp_svg_selection_group.selectAll("rect")
      .data(motif_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x", (d, i) => this.stamp_svg_margin.left + i % 3 * (stamp_width + this.stamp_svg_margin.horizontal))
      .attr("y", (d, i) => this.stamp_svg_margin.top + Math.floor(i / 3) * (stamp_height + this.stamp_svg_margin.vertical))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("fill", d => {
        if (sort_type === "Sort Reaction FDR") {
          if (d.p_values.agg <= stat_value) {
            return "green";
          } else {
            return "orange";
          }
        } else {
          if (d.p_values.source <= stat_value && d.p_values.target <= stat_value) {
            return "green";
          } else {
            return "orange";
          }
        }
      })
      .attr("id", (d) => "stamp-cover-" + d.id)
      .style("opacity", 0)
      .on("click", (d) => {
        if (timecourse === true) {
          this.generateLines(d);
        }
        document.getElementById("pathway_name").innerHTML = "<h6><b></b></h6>";
        d3.select("#pathway-view-svg").style("visibility", "hidden");
        d3.select(".network-panel").style("visibility", "hidden");

        for (let rxn in motif_list) {
          if (motif_list[rxn].id !== d.id) {
            d3.select("#stamp-cover-" + motif_list[rxn].id).style("opacity", 0);
          }
        }
        clicked_stamp_reaction = d.id;
        this.drawMotifPathway(d, indexer, motif_list);
        d3.select("#motif-pathway-svg").style("visibility", "visible");

        document.getElementById("pathway_name").innerHTML = "<h6><b>" + d.name + "</b></h6>";
        this.previewReaction(d, indexer, "#pathway-view-svg");
        d3.select("#pathway-view-svg").style("visibility", "visible");
        d3.select(".network-panel").style("visibility", "visible");
      })
      .on("mouseover", (d) => {
        d3.select("#stamp-cover-" + d.id).style("opacity", 0.5);
      })
      .on("mouseout", (d) => {
        if (d.id !== clicked_stamp_reaction) {
          d3.select("#stamp-cover-" + d.id).style("opacity", 0);
        }
      });

    let fg = this.stamp_svg_frame_group.selectAll("rect")
      .data(motif_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x", (d, i) => this.stamp_svg_margin.left + i % 3 *
        (stamp_width + this.stamp_svg_margin.horizontal))
      .attr("y", (d, i) => this.stamp_svg_margin.top + Math.floor(i / 3) *
        (stamp_height + this.stamp_svg_margin.vertical))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("stroke", "lightgray")
      .attr("fill", "white")

    let cg = this.stamp_svg_circle_group.selectAll("g")
      .data(motif_list);
    cg.exit().remove();
    cg = cg.enter().append("g").merge(cg)
      .attr("id", (d, i) => "stamp-circle-" + i);

    let lg = this.stamp_svg_link_group.selectAll("g")
      .data(motif_list);
    lg.exit().remove();
    lg = lg.enter().append("g").merge(lg)
      .attr("id", (d, i) => "stamp-link-" + i);

    for (let i = 0; i < motif_list.length; i++) {
      let mnodes = [];
      let mnodes_id = [];
      let mlinks = [];
      let r_idx = 0;
      let p_idx = 0;

      // Add reaction node
      this.nodes[motif_list[i].id].current_type = "reaction";
      mnodes.push(this.nodes[motif_list[i].id]);

      // Add reactant nodes
      motif_list[i].reactants.forEach(l => {
        if (mnodes_id.indexOf(l) === -1) {
          this.nodes[l].current_type = "reactant";
          this.nodes[l].r_idx = r_idx;
          mnodes.push(this.nodes[l]);
          mnodes_id.push(l);
          mlinks.push({
            'source': l,
            'target': motif_list[i].id
          });
          r_idx += 1;
        }
      })

      // Add product nodes
      motif_list[i].products.forEach(l => {
        if (mnodes_id.indexOf(l) === -1) {
          this.nodes[l].current_type = "product";
          this.nodes[l].p_idx = p_idx;
          mnodes.push(this.nodes[l]);
          mnodes_id.push(l);
          mlinks.push({
            'source': motif_list[i].id,
            'target': l
          });
          p_idx += 1;
        }
      })

      let start_x = this.stamp_svg_margin.left + i % 3 * (stamp_width + this.stamp_svg_margin.horizontal);
      let start_y = this.stamp_svg_margin.top + Math.floor(i / 3) * (stamp_height + this.stamp_svg_margin.vertical);
      let x_interval = stamp_width / 3;
      let y_interval_r = (stamp_height - 20) / r_idx;
      let y_interval_p = (stamp_height - 20) / p_idx;

      let mg = d3.select("#stamp-circle-" + i).selectAll("circle")
        .data(mnodes);
      mg.exit().remove();
      mg = mg.enter().append("circle").merge(mg)
        .attr("cx", (d) => {
          if (d.current_type === "reactant") {
            return start_x + 15;
          } else if (d.current_type === "reaction") {
            return start_x + x_interval + 15;
          } else if (d.current_type === "product") {
            return start_x + x_interval * 2 + 15;
          }
        })
        .attr("cy", (d) => {
          if (d.current_type === "reactant") {
            return start_y + y_interval_r * (d.r_idx) + 10;
          } else if (d.current_type === "reaction") {
            return start_y + 20;
          } else if (d.current_type === "product") {
            return start_y + y_interval_p * (d.p_idx) + 10;
          }
        })
        .attr("class", (d) => d.current_type)
        .attr("fill", (d) => {
          if (d.values_js[indexer] === undefined) {
            return "rgba(191, 191, 191, 1)";
          } else {
            return "rgba(" + d.values_js[indexer].toString() + ")";
          }
        })
        .attr("r", 4)
        .attr("stroke", "black")
        .attr("id", d => "stamp-" + i + "-" + d.id)

      let mlg = d3.select("#stamp-link-" + i).selectAll("line")
        .data(mlinks);
      mlg.exit().remove();
      mlg = mlg.enter().append("line").merge(mlg)
        .attr("x1", (d) => d3.select("#stamp-" + i + "-" + d.source).attr("cx"))
        .attr("y1", (d) => d3.select("#stamp-" + i + "-" + d.source).attr("cy"))
        .attr("x2", (d) => d3.select("#stamp-" + i + "-" + d.target).attr("cx"))
        .attr("y2", (d) => d3.select("#stamp-" + i + "-" + d.target).attr("cy"))
        .attr("stroke", "gray")

      let tg = d3.select("#stamp-circle-" + i).selectAll("text")
        .data([motif_list[i]]);
      tg.exit().remove();
      tg = tg.enter().append("text").merge(tg)
        .attr("x", start_x + 10)
        .attr("y", start_y + 45)
        .text(function(d) {
          if (sort_type === "Sort Number of Pathways") {
            return d.pathways.length + " pathways";
          } else if (sort_type === "Sort Magnitude Change" || this_stat_type === "array") {
            return "Change: " + parseFloat(d.magnitude_change).toFixed(4);
          } else if (sort_type === "Sort Statistical Significance") {
            return "Stats: " + parseFloat(d.p_values.source).toExponential(1) + " / " + parseFloat(d.p_values.target).toExponential(1);
          } else if (sort_type === "Sort Reaction FDR") {
            return "FDR: " + parseFloat(d.p_values.agg).toExponential(1);
          }
        })
        .style("font-size", "11px")
        .style("font-weight", "normal")
    }
  }

  drawTwoReactionSearchResult(motifs, indexer, exclusion_indexer) {

    let current_motifs = motifs[indexer];

    // Exclude any motifs in multiple time-points/conditions
    let motif_list = get_included_motifs(
      current_motifs,
      motifs,
      exclusion_indexer);

    if (motif_list.length > 1000) {
      alert(String(motif_list.length) + " reaction patterns identified. This may take a while, or you could try increasing the threshold to a more stringent level.")
    }

    // Sort motifs based on user input
    let sort_type = this.sort_type_dropdown.value;
    let this_stat_type = this.stat_type;
    let motif_significance = {
      'both': [],
      'one': [],
      'none': []
    };

    let sort_output = sort_motifs(
      motif_list,
      motif_significance,
      sort_type
    );
    motif_list = sort_output[0];
    motif_significance = sort_output[1];

    // Remove any duplicated items
    motif_list = remove_duplicate_motifs_twoReactions(motif_list);
    motif_significance = remove_duplicate_motifs_twoReactions(motif_significance);

    // Remove if no shared components
    motif_list = remove_noshare_twoReactions(motif_list);
    motif_significance = remove_noshare_twoReactions(motif_significance);

    let stamp_height = 50;
    this.stamp_svg_height = Math.ceil(
      motif_list.length / 3) *
      (stamp_height + this.stamp_svg_margin.vertical);
    this.stamp_svg.attr("height", this.stamp_svg_height);

    let stamp_width = this.stamp_svg_width / 3 -
      this.stamp_svg_margin.horizontal;

    let clicked_stamp_reaction;
    let sg = this.stamp_svg_selection_group.selectAll("rect")
      .data(motif_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x", (d, i) => this.stamp_svg_margin.left + i % 3 * (stamp_width + this.stamp_svg_margin.horizontal))
      .attr("y", (d, i) => this.stamp_svg_margin.top + Math.floor(i / 3) * (stamp_height + this.stamp_svg_margin.vertical))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("fill", d => {
        if (d.p_values.source <= stat_value && d.p_values.target <= stat_value) {
          return "green";
        } else {
          return "orange";
        }
      })
      .attr("id", (d) => "stamp-cover-" + d.id)
      .style("opacity", 0)
      .on("click", (d) => {
        if (timecourse === true) {
          this.generateLines(d);
        }
        document.getElementById("pathway_name").innerHTML = "<h6><b></b></h6>";
        d3.select("#pathway-view-svg").style("visibility", "hidden");
        d3.select(".network-panel").style("visibility", "hidden");

        for (let rxn in motif_list) {
          if (motif_list[rxn].id !== d.id) {
            d3.select("#stamp-cover-" + motif_list[rxn].id).style("opacity", 0);
          }
        }
        clicked_stamp_reaction = d.id;
        this.drawMotifPathwayTwoReaction(d, indexer, motif_list);
        d3.select("#motif-pathway-svg").style("visibility", "visible");

        document.getElementById("pathway_name").innerHTML = "<h6><b>" + d['rxn1'].name + " // " + d['rxn2'].name + "</b></h6>";

        this.previewTwoReaction(d, indexer, "#pathway-view-svg");
        d3.select("#pathway-view-svg").style("visibility", "visible");
        d3.select(".network-panel").style("visibility", "visible");
      })
      .on("mouseover", (d) => {
        d3.select("#stamp-cover-" + d.id).style("opacity", 0.5);
      })
      .on("mouseout", (d) => {
        if (d.id !== clicked_stamp_reaction) {
          d3.select("#stamp-cover-" + d.id).style("opacity", 0);
        }
      });

    let fg = this.stamp_svg_frame_group.selectAll("rect")
      .data(motif_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x", (d, i) => this.stamp_svg_margin.left + i % 3 *
        (stamp_width + this.stamp_svg_margin.horizontal))
      .attr("y", (d, i) => this.stamp_svg_margin.top + Math.floor(i / 3) *
        (stamp_height + this.stamp_svg_margin.vertical))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("stroke", "lightgray")
      .attr("fill", "white")

    let cg = this.stamp_svg_circle_group.selectAll("g")
      .data(motif_list);
    cg.exit().remove();
    cg = cg.enter().append("g").merge(cg)
      .attr("id", (d, i) => "stamp-circle-" + i);

    let lg = this.stamp_svg_link_group.selectAll("g")
      .data(motif_list);
    lg.exit().remove();
    lg = lg.enter().append("g").merge(lg)
      .attr("id", (d, i) => "stamp-link-" + i);

    for (let i = 0; i < motif_list.length; i++) {
      let mnodes = [];
      let mnodes_id = [];
      let mlinks = [];
      let r_idx = 0;
      let p_idx = 0;

      // Add reaction node
      this.nodes[motif_list[i]['rxn1'].id].current_type = "reaction";
      mnodes.push(this.nodes[motif_list[i]['rxn1'].id]);

      // Add reactant nodes
      motif_list[i]['rxn1'].reactants.forEach(l => {
        if (mnodes_id.indexOf(l) === -1) {
          this.nodes[l].current_type = "reactant";
          this.nodes[l].r_idx = r_idx;
          mnodes.push(this.nodes[l]);
          mnodes_id.push(l);
          mlinks.push({
            'source': l,
            'target': motif_list[i]['rxn1'].id
          });
          r_idx += 1;
        }
      })

      // Add product nodes
      motif_list[i]['rxn1'].products.forEach(l => {
        if (mnodes_id.indexOf(l) === -1) {
          this.nodes[l].current_type = "product";
          this.nodes[l].p_idx = p_idx;
          mnodes.push(this.nodes[l]);
          mnodes_id.push(l);
          mlinks.push({
            'source': motif_list[i]['rxn1'].id,
            'target': l
          });
          p_idx += 1;
        }
      })

      let start_x = this.stamp_svg_margin.left + i % 3 * (stamp_width + this.stamp_svg_margin.horizontal);
      let start_y = this.stamp_svg_margin.top + Math.floor(i / 3) * (stamp_height + this.stamp_svg_margin.vertical);
      let x_interval = stamp_width / 3;
      let y_interval_r = (stamp_height - 20) / r_idx;
      let y_interval_p = (stamp_height - 20) / p_idx;

      let mg = d3.select("#stamp-circle-" + i).selectAll("circle")
        .data(mnodes);
      mg.exit().remove();
      mg = mg.enter().append("circle").merge(mg)
        .attr("cx", (d) => {
          if (d.current_type === "reactant") {
            return start_x + 15;
          } else if (d.current_type === "reaction") {
            return start_x + x_interval + 15;
          } else if (d.current_type === "product") {
            return start_x + x_interval * 2 + 15;
          }
        })
        .attr("cy", (d) => {
          if (d.current_type === "reactant") {
            return start_y + y_interval_r * (d.r_idx) + 10;
          } else if (d.current_type === "reaction") {
            return start_y + 20;
          } else if (d.current_type === "product") {
            return start_y + y_interval_p * (d.p_idx) + 10;
          }
        })
        .attr("class", (d) => d.current_type)
        .attr("fill", (d) => {
          if (d.values_js[indexer] === undefined) {
            return "rgba(191, 191, 191, 1)";
          } else {
            return "rgba(" + d.values_js[indexer].toString() + ")";
          }
        })
        .attr("r", 4)
        .attr("stroke", "black")
        .attr("id", d => "stamp-" + i + "-" + d.id)

      let mlg = d3.select("#stamp-link-" + i).selectAll("line")
        .data(mlinks);
      mlg.exit().remove();
      mlg = mlg.enter().append("line").merge(mlg)
        .attr("x1", (d) => d3.select("#stamp-" + i + "-" + d.source).attr("cx"))
        .attr("y1", (d) => d3.select("#stamp-" + i + "-" + d.source).attr("cy"))
        .attr("x2", (d) => d3.select("#stamp-" + i + "-" + d.target).attr("cx"))
        .attr("y2", (d) => d3.select("#stamp-" + i + "-" + d.target).attr("cy"))
        .attr("stroke", "gray")

      let tg = d3.select("#stamp-circle-" + i).selectAll("text")
        .data([motif_list[i]]);
      tg.exit().remove();
      tg = tg.enter().append("text").merge(tg)
        .attr("x", start_x + 10)
        .attr("y", start_y + 45)
        .text(function(d) {
          console.log(this_stat_type)
          if (sort_type === "Sort Number of Pathways") {
            return d.pathways.length + " pathways";
          } else if (sort_type === "Sort Magnitude Change" || this_stat_type === "array") {
            return "Change: " + parseFloat(d.magnitude_change).toFixed(4);
          } else if (sort_type === "Sort Statistical Significance") {
            return "Stats: " + parseFloat(d.p_values.source).toExponential(1) + " / " + parseFloat(d.p_values.target).toExponential(1);
          } else if (sort_type === "Sort Reaction FDR") {
            return "FDR: " + parseFloat(d.p_values.agg).toExponential(1);
          }
        })
        .style("font-size", "11px")
        .style("font-weight", "normal")
    }
  }

  drawMotifSearchResultPathway(motif_list) {
    // modifier for path coverage motif
    var coverage_true = false;
    var motif_dict = {};
    if (motif_list[0].constructor === Array) {
      coverage_true = true;
    }

    if (coverage_true === true) {
      var original_list = motif_list;
      motif_list = [];
      for (let x in original_list) {
        motif_list.push(original_list[x][0])
        motif_dict[original_list[x][0]["id"]] = String(original_list[x][2]).concat("/", String(original_list[x][3]));
      }
    }

    let stamp_height = 30;
    this.stamp_svg_height = Math.ceil(
        motif_list.length) *
      (stamp_height + this.stamp_svg_margin.vertical);
    this.stamp_svg.attr("height", this.stamp_svg_height);

    let stamp_width = this.stamp_svg_width -
      this.stamp_svg_margin.horizontal - 5;

    let clicked_stamp_reaction;
    let sg = this.stamp_svg_selection_group.selectAll("rect")
      .data(motif_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x", 5)
      .attr("y", (d, i) => this.stamp_svg_margin.top + (i * (stamp_height + this.stamp_svg_margin.vertical)))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("fill", "lightgray")
      .attr("id", (d) => "stamp-cover-" + d.id)
      .style("opacity", 0)
      .on("click", (d) => {
        for (let rxn in motif_list) {
          if (motif_list[rxn].id !== d.id) {
            d3.select("#stamp-cover-" + motif_list[rxn].id).style("opacity", 0);
          }
        }
        clicked_stamp_reaction = d.id;
        document.getElementById("pathway_name").innerHTML = "<h6><b>" + d.name + "</b></h6>";
        if (typeof(d.id) === "string" && d.id !== "Collapsed") {
          this.drawPathwayView(d.id, "#pathway-view-svg", motif_list);
          d3.select("#pathway-view-svg").style("visibility", "visible");
          d3.select(".network-panel").style("visibility", "visible");
        }
      })
      .on("mouseover", (d) => {
        d3.select("#stamp-cover-" + d.id).style("opacity", 0.4);
      })
      .on("mouseout", (d) => {
        if (d.id !== clicked_stamp_reaction) {
          d3.select("#stamp-cover-" + d.id).style("opacity", 0);
        }
      })

    let fg = this.stamp_svg_frame_group.selectAll("rect")
      .data(motif_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x", 5)
      .attr("y", (d, i) => this.stamp_svg_margin.top + Math.floor(i) *
        (stamp_height + this.stamp_svg_margin.vertical))
      .attr("width", stamp_width)
      .attr("height", stamp_height)
      .attr("stroke", "lightgray")
      .attr("fill", "white")

    // remove previous elements, ignoring timecourse slider
    d3.selectAll("circle:not(#dot)").remove();
    d3.selectAll("line:not(#track)").remove();
    d3.selectAll("text:not(#tick)").remove();

    let cg = this.stamp_svg_circle_group.selectAll("g")
      .data(motif_list);
    cg.exit().remove();
    cg = cg.enter().append("g").merge(cg)
      .attr("id", (d, i) => "stamp-circle-" + i);
    console.log(cg)
    for (let i = 0; i < motif_list.length; i++) {

      let start_x = 5;
      let start_y = this.stamp_svg_margin.top +
        Math.floor(i) * (stamp_height + this.stamp_svg_margin.vertical);

      let tg = d3.select("#stamp-circle-" + i).selectAll("text")
        .data([motif_list[i]]);
      tg.exit().remove();
      tg = tg.enter().append("text").merge(tg)
        .attr("x", start_x + 10)
        .attr("y", start_y + 21)
        .text(d => {
          if (coverage_true === true) {
            if (d.name.length < 40) {
              return d.name + " (" + motif_dict[d.id] + ")";
            } else {
              return d.name.substring(0, 40) + "... (" + motif_dict[d.id] + ")";
            }
          } else {
            if (d.name.length < 45) {
              return d.name;
            } else {
              return d.name.substring(0, 45) + " ...";
            }
          }
        })
        .style("font-size", "14px")
        .style("font-weight", "bold")
    }
  }

  drawMotifPathway(motif, indexer, motif_list) {

    // recover graph
    let mnodes = [motif];
    let mlinks = [];
    let mnodes_id = [];
    let r_idx = 0;
    let p_idx = 0;

    motif.current_type = "reaction";
    motif.reactants.forEach(l => {
      if (mnodes_id.indexOf(l) === -1) {
        this.nodes[l].current_type = "reactant";
        this.nodes[l].r_idx = r_idx;
        mnodes.push(this.nodes[l]);
        mnodes_id.push(l);
        mlinks.push({
          'source': l,
          'target': motif.id
        });
        r_idx += 1;
      }
    })
    motif.products.forEach(l => {
      if (mnodes_id.indexOf(l) === -1) {
        this.nodes[l].current_type = "product";
        this.nodes[l].p_idx = p_idx;
        mnodes.push(this.nodes[l]);
        mnodes_id.push(l);
        mlinks.push({
          'source': motif.id,
          'target': l
        });
        p_idx += 1;
      }
    })

    let motif_height = 120;

    // ***** draw motif glyph *****
    let x_interval = this.mp_svg_width / 3;
    let y_interval_r = (motif_height - 20) / r_idx;
    let y_interval_p = (motif_height - 20) / p_idx;
    let ng = this.mp_motif_circle_group.selectAll("circle")
      .data(mnodes)
    ng.exit().remove();
    ng = ng.enter().append("circle").merge(ng)
      .attr("cx", (d) => {
        if (d.current_type === "reactant") {
          return 75;
        } else if (d.current_type === "reaction") {
          return x_interval + 75;
        } else if (d.current_type === "product") {
          return x_interval * 2 + 75;
        }
      })
      .attr("cy", (d) => {
        if (d.current_type === "reactant") {
          return y_interval_r * (d.r_idx) + 20;
        } else if (d.current_type === "reaction") {
          return motif_height / 2;
        } else if (d.current_type === "product") {
          return y_interval_p * (d.p_idx) + 20;
        }
      })
      .attr("class", (d) => d.current_type)
      .attr("fill", (d) => {
        if (d.current_type === "reaction") {
          return "rgba(191, 191, 191, 1)";
        } else {
          return "rgba(" + d.values_js[indexer].toString() + ")";
        }
      })
      .attr("r", 8)
      .attr("stroke", "black")
      .attr("id", d => "mp-circle-" + d.id)

    let lg = this.mp_motif_link_group.selectAll("line")
      .data(mlinks);
    lg.exit().remove();
    lg = lg.enter().append("line").merge(lg)
      .attr("x1", (d) => d3.select("#mp-circle-" + d.source).attr("cx"))
      .attr("y1", (d) => d3.select("#mp-circle-" + d.source).attr("cy"))
      .attr("x2", (d) => d3.select("#mp-circle-" + d.target).attr("cx"))
      .attr("y2", (d) => d3.select("#mp-circle-" + d.target).attr("cy"))
      .attr("stroke", "gray")
      .attr("marker-end", "url(#end)");

    let tg = this.mp_motif_circle_group.selectAll("text")
      .data([motif]);
    tg.exit().remove();
    tg = tg.enter().append("text").merge(tg)
      .attr("x", this.mp_svg_width / 12 - 10)
      .attr("y", (motif_height - 5))
      .text(function(d) {
        if (d.collapsed === "true") {
          var string_length = 40;
          var string_suffix = " (collapsed)";
        } else {
          var string_length = 52;
          var string_suffix = "";
        }
        if (d.name.length < string_length) {
          return d.name + string_suffix;
        } else {
          return d.name.substring(0, string_length) + " ..." + string_suffix;
        }
      })
      .style("font-size", "12px")
      .style("font-weight", "bold")

    // **** draw pathway glyph ****
    let pathway_list_pre = [...new Set(motif.pathways)];
    let pathway_list = [];
    for (let path in pathway_list_pre) {
      if (!this.superPathwayDict.includes(pathway_list_pre[path])) {
        pathway_list.push(pathway_list_pre[path]);
      }
    }

    if (pathway_list.length === 0) {
      pathway_list.push("Collapsed")
    }

    let pathway_height = 25;
    let margin = {
      "horizontal": 15,
      "vertical": 10,
      "top": 5,
      "left": 5
    };
    this.mp_svg_height = Math.ceil(pathway_list.length) * (pathway_height + margin.vertical) + motif_height + margin.top;
    this.mp_svg.attr("height", this.mp_svg_height);

    let pathway_width = (this.mp_svg_width) - margin.horizontal;

    let clicked_stamp_pathway;
    let sg = this.mp_selection_group.selectAll("rect")
      .data(pathway_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x", 5)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical)))
      .attr("width", pathway_width)
      .attr("height", pathway_height)
      .attr("fill", "lightblue")
      .attr("id", (d) => "mp-cover-" + d)
      .style("opacity", 0)
      .on("click", (d) => {
        for (let path in pathway_list) {
          if (pathway_list[path] !== d) {
            d3.select("#mp-cover-" + pathway_list[path]).style("opacity", 0);
          }
        }
        clicked_stamp_pathway = d;
        if (typeof(d) === "string" && d !== "Collapsed") {
          if (d.length !== 0) {
            this.drawPathwayView(d, "#pathway-view-svg", motif_list);
            this.findAllMotif(d, this.motif[indexer]);
          } else {
            document.getElementById("pathway_name").innerHTML = "<h6><b>Collapsed Reaction</h6></b>";
            this.drawPathwayView(motif, "#pathway-view-svg", motif_list);
            this.findAllMotif(motif, this.motif[indexer]);
          }
          d3.select("#pathway-view-svg").style("visibility", "visible");
          d3.select(".network-panel").style("visibility", "visible");
        }
      })
      .on("mouseover", (d) => {
        d3.select("#mp-cover-" + d).style("opacity", 0.4);
      })
      .on("mouseout", (d) => {
        if (d !== clicked_stamp_pathway) {
          d3.select("#mp-cover-" + d).style("opacity", 0);
        }
      })

    let fg = this.mp_pathway_group.selectAll("rect")
      .data(pathway_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x", 5)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical)))
      .attr("width", pathway_width)
      .attr("height", pathway_height)
      .attr("fill", "lightgrey")
      .style("opacity", 0.5);

    // get pathway names
    let ptg = this.mp_pathway_group.selectAll("text")
      .data(pathway_list);
    ptg.exit().remove();
    ptg = ptg.enter().append("text").merge(ptg)
      .attr("x", 20)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical) + 17))
      .text(d => {
        if (d === "Collapsed") {
          return "Cross-pathway pattern";
        } else if (this.mod_collapsed_pathways[d] !== undefined) {
          return this.mod_collapsed_pathways[d].name.substring(0, 62);
        } else {
          return "";
        }
      })
      .style("font-size", 12)
  }

  drawMotifPathwayTwoReaction(motif, indexer, motif_list) {

    // recover graph
    let mnodes = [motif['rxn1']];
    let mlinks = [];
    let mnodes_id = [];
    let r_idx = 0;
    let p_idx = 0;

    motif['rxn1'].current_type = "reaction";
    motif['rxn1'].reactants.forEach(l => {
      if (mnodes_id.indexOf(l) === -1) {
        this.nodes[l].current_type = "reactant";
        this.nodes[l].r_idx = r_idx;
        mnodes.push(this.nodes[l]);
        mnodes_id.push(l);
        mlinks.push({
          'source': l,
          'target': motif['rxn1'].id
        });
        r_idx += 1;
      }
    })
    motif['rxn1'].products.forEach(l => {
      if (mnodes_id.indexOf(l) === -1) {
        this.nodes[l].current_type = "product";
        this.nodes[l].p_idx = p_idx;
        mnodes.push(this.nodes[l]);
        mnodes_id.push(l);
        mlinks.push({
          'source': motif['rxn1'].id,
          'target': l
        });
        p_idx += 1;
      }
    })

    let motif_height = 120;

    // ***** draw motif glyph *****
    let x_interval = this.mp_svg_width / 3;
    let y_interval_r = (motif_height - 20) / r_idx;
    let y_interval_p = (motif_height - 20) / p_idx;
    let ng = this.mp_motif_circle_group.selectAll("circle")
      .data(mnodes)
    ng.exit().remove();
    ng = ng.enter().append("circle").merge(ng)
      .attr("cx", (d) => {
        if (d.current_type === "reactant") {
          return 75;
        } else if (d.current_type === "reaction") {
          return x_interval + 75;
        } else if (d.current_type === "product") {
          return x_interval * 2 + 75;
        }
      })
      .attr("cy", (d) => {
        if (d.current_type === "reactant") {
          return y_interval_r * (d.r_idx) + 20;
        } else if (d.current_type === "reaction") {
          return motif_height / 2;
        } else if (d.current_type === "product") {
          return y_interval_p * (d.p_idx) + 20;
        }
      })
      .attr("class", (d) => d.current_type)
      .attr("fill", (d) => {
        if (d.current_type === "reaction") {
          return "rgba(191, 191, 191, 1)";
        } else {
          return "rgba(" + d.values_js[indexer].toString() + ")";
        }
      })
      .attr("r", 8)
      .attr("stroke", "black")
      .attr("id", d => "mp-circle-" + d.id)

    let lg = this.mp_motif_link_group.selectAll("line")
      .data(mlinks);
    lg.exit().remove();
    lg = lg.enter().append("line").merge(lg)
      .attr("x1", (d) => d3.select("#mp-circle-" + d.source).attr("cx"))
      .attr("y1", (d) => d3.select("#mp-circle-" + d.source).attr("cy"))
      .attr("x2", (d) => d3.select("#mp-circle-" + d.target).attr("cx"))
      .attr("y2", (d) => d3.select("#mp-circle-" + d.target).attr("cy"))
      .attr("stroke", "gray")
      .attr("marker-end", "url(#end)");

    let tg = this.mp_motif_circle_group.selectAll("text")
      .data([motif]);
    tg.exit().remove();
    tg = tg.enter().append("text").merge(tg)
      .attr("x", this.mp_svg_width / 12 - 10)
      .attr("y", (motif_height - 5))
      .text(function(d) {
        if (d.collapsed === "true") {
          var string_length = 40;
          var string_suffix = " (collapsed)";
        } else {
          var string_length = 52;
          var string_suffix = "";
        }
        let this_name = d['rxn1'].name + " // " + d['rxn2'].name;
        if (this_name.length < string_length) {
          return this_name + string_suffix;
        } else {
          return this_name.substring(0, string_length) + " ..." + string_suffix;
        }
      })
      .style("font-size", "12px")
      .style("font-weight", "bold")

    // **** draw pathway glyph ****
    let pathway_list_pre = [...new Set(motif.pathways)];
    let pathway_list = [];
    for (let path in pathway_list_pre) {
      if (!this.superPathwayDict.includes(pathway_list_pre[path])) {
        pathway_list.push(pathway_list_pre[path]);
      }
    }

    if (pathway_list.length === 0) {
      pathway_list.push("Collapsed")
    }

    let pathway_height = 25;
    let margin = {
      "horizontal": 15,
      "vertical": 10,
      "top": 5,
      "left": 5
    };
    this.mp_svg_height = Math.ceil(pathway_list.length) * (pathway_height + margin.vertical) + motif_height + margin.top;
    this.mp_svg.attr("height", this.mp_svg_height);

    let pathway_width = (this.mp_svg_width) - margin.horizontal;

    let clicked_stamp_pathway;
    let sg = this.mp_selection_group.selectAll("rect")
      .data(pathway_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x", 5)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical)))
      .attr("width", pathway_width)
      .attr("height", pathway_height)
      .attr("fill", "lightblue")
      .attr("id", (d) => "mp-cover-" + d)
      .style("opacity", 0)
      .on("click", (d) => {
        for (let path in pathway_list) {
          if (pathway_list[path] !== d) {
            d3.select("#mp-cover-" + pathway_list[path]).style("opacity", 0);
          }
        }
        clicked_stamp_pathway = d;
        if (typeof(d) === "string" && d !== "Collapsed") {
          if (d.length !== 0) {
            this.drawPathwayView(d, "#pathway-view-svg", motif_list);
            this.findAllMotif(d, this.motif[indexer]);
          } else {
            document.getElementById("pathway_name").innerHTML = "<h6><b>Collapsed Reaction</h6></b>";
            this.drawPathwayView(motif, "#pathway-view-svg", motif_list);
            this.findAllMotif(motif, this.motif[indexer]);
          }
          d3.select("#pathway-view-svg").style("visibility", "visible");
          d3.select(".network-panel").style("visibility", "visible");
        }
      })
      .on("mouseover", (d) => {
        d3.select("#mp-cover-" + d).style("opacity", 0.4);
      })
      .on("mouseout", (d) => {
        if (d !== clicked_stamp_pathway) {
          d3.select("#mp-cover-" + d).style("opacity", 0);
        }
      })

    let fg = this.mp_pathway_group.selectAll("rect")
      .data(pathway_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x", 5)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical)))
      .attr("width", pathway_width)
      .attr("height", pathway_height)
      .attr("fill", "lightgrey")
      .style("opacity", 0.5);

    // get pathway names
    let ptg = this.mp_pathway_group.selectAll("text")
      .data(pathway_list);
    ptg.exit().remove();
    ptg = ptg.enter().append("text").merge(ptg)
      .attr("x", 20)
      .attr("y", (d, i) => motif_height + margin.top + (i * (pathway_height + margin.vertical) + 17))
      .text(d => {
        if (d === "Collapsed") {
          return "Cross-pathway pattern";
        } else if (this.mod_collapsed_pathways[d] !== undefined) {
          return this.mod_collapsed_pathways[d].name.substring(0, 62);
        } else {
          return "";
        }
      })
      .style("font-size", 12)
  }

  previewReaction(d, indexer, selector) {

    let components = [];
    var rxn = 0;
    components.push(d.id);
    for (let x in d["reactants"]) {
      components.push(d["reactants"][x]);
    }
    for (let x in d["products"]) {
      components.push(d["products"][x]);
    }
    for (let x in d["modifiers"]) {
      components.push(d["modifiers"][x][0]);
    }
    for (let x in d["additional_components"]) {
      components.push(d["additional_components"][x]);
    }

    var elements = get_nodes_links(this.data, components);
    var new_nodes = elements[0];
    var new_links = elements[1];

    // Initialize variables
    var node_dict = {};
    var type_dict = {};

    var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
    var node_dict = node_elements[0];
    var type_dict = node_elements[1];
    var display_analytes_dict = node_elements[2];
    var display_reactions_dict = node_elements[3];
    var entity_id_dict = node_elements[4];

    make_graph(
      this.data,
      new_nodes,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      display_analytes_dict,
      display_reactions_dict,
      selector,
      this.stat_type,
      _width,
      _height,
      [d]
    );
  }

  previewTwoReaction(d, indexer, selector) {

    let components = [];
    var rxn = 0;
    components.push(d['rxn1'].id);
    for (let x in d['rxn1']["reactants"]) {
      components.push(d['rxn1']["reactants"][x]);
    }
    for (let x in d['rxn1']["products"]) {
      components.push(d['rxn1']["products"][x]);
    }
    for (let x in d['rxn1']["modifiers"]) {
      components.push(d['rxn1']["modifiers"][x][0]);
    }
    for (let x in d['rxn1']["additional_components"]) {
      components.push(d['rxn1']["additional_components"][x]);
    }
    components.push(d['rxn2'].id);
    for (let x in d['rxn2']["reactants"]) {
      components.push(d['rxn2']["reactants"][x]);
    }
    for (let x in d['rxn2']["products"]) {
      components.push(d['rxn2']["products"][x]);
    }
    for (let x in d['rxn2']["modifiers"]) {
      components.push(d['rxn2']["modifiers"][x][0]);
    }
    for (let x in d['rxn2']["additional_components"]) {
      components.push(d['rxn2']["additional_components"][x]);
    }

    var elements = get_nodes_links(this.data, components);
    var new_nodes = elements[0];
    var new_links = elements[1];

    // Initialize variables
    var node_dict = {};
    var type_dict = {};

    var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
    var node_dict = node_elements[0];
    var type_dict = node_elements[1];
    var display_analytes_dict = node_elements[2];
    var display_reactions_dict = node_elements[3];
    var entity_id_dict = node_elements[4];

    make_graph(
      this.data,
      new_nodes,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      display_analytes_dict,
      display_reactions_dict,
      selector,
      this.stat_type,
      _width,
      _height,
      [d]
    );
  }

  drawPathwayView(p, selector, motif_list) {

    graph_genes = true;
    collapse_reactions = true;
    if (typeof(p) === "string") {
      var motif_reactions = this.mod_collapsed_pathways[p]["reactions"];
      update_session_info("current_pathway", p);
    } else {
      var motif_reactions = [p.id];
      update_session_info("current_pathway", null);
    }

    // Parse through each reaction listed and get the component parts
    let components = [];
    var rxn = 0;
    for (rxn in motif_reactions) {

      var target_rxns = this.collapsed_reaction_dict[motif_reactions[rxn]];
      components.push(motif_reactions[rxn]);
      for (let x in target_rxns["reactants"]) {
        components.push(target_rxns["reactants"][x]);
      }
      for (let x in target_rxns["products"]) {
        components.push(target_rxns["products"][x]);
      }
      for (let x in target_rxns["modifiers"]) {
        components.push(target_rxns["modifiers"][x][0]);
      }
      for (let x in target_rxns["additional_components"]) {
        components.push(target_rxns["additional_components"][x]);
      }
    }

    var elements = get_nodes_links(this.data, components);
    var new_nodes = elements[0];
    var new_links = elements[1];

    // Initialize variables
    var node_dict = {};
    var type_dict = {};

    var node_elements = initialize_nodes(new_nodes, node_dict, type_dict);
    var node_dict = node_elements[0];
    var type_dict = node_elements[1];
    var display_analytes_dict = node_elements[2];
    var display_reactions_dict = node_elements[3];
    var entity_id_dict = node_elements[4];

    make_graph(
      this.data,
      new_nodes,
      new_links,
      type_dict,
      node_dict,
      entity_id_dict,
      display_analytes_dict,
      display_reactions_dict,
      selector,
      this.stat_type,
      _width,
      _height,
      motif_list
    );
  }

  findAllMotif(pathway, motif_list) {

    d3.select(".network-panel").style("visibility", "visible");

    let current_motif = [];
    if (typeof(pathway) === "string") {
      let current_pathway = this.mod_collapsed_pathways[pathway].name;
      document.getElementById("pathway_name").innerHTML = "<h6><b>" + current_pathway + "</b></h6>";

      let motif_id = [];
      motif_list.forEach(m => {
        motif_id.push(m.id);
      })
      let pathway_list = this.mod_collapsed_pathways[pathway].reactions;

      pathway_list.forEach(react_id => {
        if (motif_id.indexOf(react_id) != -1) {
          current_motif.push(react_id);
        }
      })
      current_motif.forEach(react_id => {
        d3.selectAll("circle#" + react_id)
          .style("r", "16px")
          .style("stroke", "purple")
          .style("stroke-width", "5px")

      })
    }
  }

  createId(id) {
    return id.replace(/[^a-zA-Z0-9]/g, "")
  }

  generateLines(d, stat_threshold) {

    d3.select("svg#line-graph").remove();
    let margin = {
      top: 50,
      right: 500,
      bottom: 50,
      left: 75
    };
    let width = 1110;
    let height = 405;
    let x_offset = 30;
    let y_offset = ((height / 2) - 40);
    let names = this.labels.split(',');
    let n = this.categories.length - 1;

    let xScale = d3.scaleLinear()
      .domain([0, n])
      .range([0, width - 170]);

    function extract_values(nodes, items, all_values, modifier=false) {
      // Extract measurements values from reaction
      for (let i in items) {
        if (modifier === false) {
          for (let x in nodes[items[i]].values) {
            let _x = nodes[items[i]].values[x];
            if (_x != null) {
              all_values.push(_x);
            }
          }
        } else {
          for (let x in nodes[items[i][0]].values) {
            let _x = nodes[items[i][0]].values[x];
            if (_x != null) {
              all_values.push(_x);
            }
          }
        }
      }
      return all_values;
    } 

    let all_values = [];
    all_values = extract_values(this.nodes, d.reactants, all_values);
    all_values = extract_values(this.nodes, d.products, all_values);
    all_values = extract_values(this.nodes, d.modifiers, all_values, true);
    all_values = extract_values(this.nodes, d.additional_components, all_values);
    
    let lower_bound = Math.min(...all_values) - 1;
    let upper_bound = Math.max(...all_values) + 1;

    let yScale = d3.scaleLinear()
      .domain([lower_bound, upper_bound])
      .range([height - 80, 0]);

    let line_svg = d3.select("#line-plot-container")
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("id", "line-graph")
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    line_svg
      .append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(" + x_offset + "," + y_offset * 2 + ")")
      .call(d3.axisBottom(xScale)
        .tickValues(this.categories)
        .tickFormat(function(d, i) {
          return names[i]
        })
      );
    line_svg
      .append("g")
      .attr("class", "y axis")
      .call(d3.axisLeft(yScale));

    for (let a in d.additional_components) {
      let _i = this.nodes[d.additional_components[a]];
      let _n = _i.name.replace(/,/g, '_').replace(/:/g, '_');
      let _n_type = _i.name;
      let _t = _i.type;
      let _st = _i.sub_type;
      let _v = _i.values;
      let _s = _i.stats;
      let _inf = _i.inferred;
      let _v_ = [_v, _s];
      if (_v[0] != null && _t !== "complex_component" && _inf !== "true") {
        let dash_instruction;
        let _set = new Set(_v)
        if (_set.size === 1) {
          dash_instruction = "15 5"
        } else {
          dash_instruction = "1 0"
        }
        let current_line = d3.line()
          .x(function(d, i) {
            return xScale(i) + x_offset;
          })
          .y(function(d, i) {
            return yScale(d);
          })
          .curve(d3.curveMonotoneX);
        line_svg
          .append("path")
          .datum(_v)
          .attr("id", function() {
            return "line_" + _n + "_" + _t
          })
          .attr("class", "line")
          .style("fill", "none")
          .style("stroke-width", "3px")
          .style("stroke", function() {
            if (_t === "gene_component" || _st === "gene_component") {
              return "rgba(165, 55, 253, 1)";
            } else if (_t === "protein_component" || _st === "protein_component") {
              return "orange";
            } else if (_t === "metabolite_component" || _st === "metabolite_component") {
              return "blue";
            } else if (_t === "complex_component" || _st === "complex_component") {
              return "black";
            } else {
              return "#808080";
            }
          })
          .style("stroke-dasharray", dash_instruction)
          .attr("d", current_line)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        line_svg
          .append("text")
          .attr("id", _n + "_" + _t + "_text")
          .attr("x", xScale(_v.length - 1) + x_offset + 15)
          .attr("y", yScale(_v[_v.length - 1]) + 5)
          .style("fill", function() {
            if (_t === "gene_component" || _st === "gene_component") {
              return "rgba(165, 55, 253, 1)";
            } else if (_t === "protein_component" || _st === "protein_component") {
              return "orange";
            } else if (_t === "metabolite_component" || _st === "metabolite_component") {
              return "blue";
            } else if (_t === "complex_component" || _st === "complex_component") {
              return "black";
            } else {
              return "#808080";
            }
          })
          .text(_n_type)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        // add nodes and info on hover
        for (let _i_ in this.categories) {
          let _v_ = _v[_i_];
          let _s_ = _s[_i_];
          line_svg
            .append("path")
            .attr("id", function() {
              return _n + "_" + _i_
            })
            .style("fill", "white")
            .style("stroke", "black")
            .style("stroke-width", function() {
              if ((_s_ === undefined) || (_s_ === null)) {
                return nonsignificance_weight;
              } else if (_s_ < stat_value) {
                return significance_weight;
              } else {
                return nonsignificance_weight;
              }
            })
            .attr("d", d3.symbol()
              .size(function() {
                return 100;
              })
              .type(function() {
                if (_t === "gene_component" || _st === "gene_component") {
                  return d3.symbolTriangle;
                } else if (_t === "protein_component" || _st === "protein_component") {
                  return d3.symbolDiamond;
                } else if (_t === "metabolite_component" || _st === "metabolite_component") {
                  return d3.symbolCircle;
                } else if (_t === "complex_component" || _st === "complex_component") {
                  return d3.symbolSquare;
                } else {
                  return d3.symbolCross;
                }
              }))
            .attr("transform", "translate(" + (xScale(_i_) + x_offset) + "," + (yScale(_v_)) + ")")
        }
      }
    }
    for (let m in d.modifiers) {
      let _m = d.modifiers[m][1];
      let _i = this.nodes[d.modifiers[m][0]];
      let _n = _i.name.replace(/,/g, '_').replace(/:/g, '_');
      let _n_type = _i.name;
      let _t = _i.type;
      let _st = _i.sub_type;
      let _v = _i.values;
      let _s = _i.stats;
      let _inf = _i.inferred;
      let _v_ = [_v, _s];
      if (_v[0] != null && _st !== "complex_component" && _inf !== "true") {
        let dash_instruction;
        let _set = new Set(_v)
        if (_set.size === 1) {
          dash_instruction = "15 5"
        } else {
          dash_instruction = "1 0"
        }
        let current_line = d3.line()
          .x(function(d, i) {
            return xScale(i) + x_offset;
          })
          .y(function(d, i) {
            return yScale(d);
          })
          .curve(d3.curveMonotoneX);
        if (_m === "catalyst") {
          line_svg
            .append("path")
            .datum(_v)
            .attr("id", function() {
              return "line_" + _n + "_" + _t
            })
            .attr("class", "line")
            .style("fill", "none")
            .style("stroke-width", "3px")
            .style("stroke", "#1b9e77")
            .style("stroke-dasharray", dash_instruction)
            .attr("d", current_line)
            .on("mouseover", function() {
              let _this;
              _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
              _this.parentNode.appendChild(_this);
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "5px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "18px")
            })
            .on("mouseout", function() {
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "3px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "12px")
            });
          line_svg
            .append("text")
            .attr("id", _n + "_" + _t + "_text")
            .attr("x", xScale(_v.length - 1) + x_offset + 15)
            .attr("y", yScale(_v[_v.length - 1]) + 5)
            .style("fill", "#1b9e77")
            .text(_n_type)
            .on("mouseover", function() {
              let _this;
              _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
              _this.parentNode.appendChild(_this);
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "5px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "18px")
            })
            .on("mouseout", function() {
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "3px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "12px")
            });
        } else {
          line_svg
            .append("path")
            .datum(_v)
            .attr("id", function() {
              return "line_" + _n + "_" + _t
            })
            .attr("class", "line")
            .style("fill", "none")
            .style("stroke-width", "3px")
            .style("stroke", "red")
            .style("stroke-dasharray", dash_instruction)
            .attr("d", current_line)
            .on("mouseover", function() {
              let _this;
              _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
              _this.parentNode.appendChild(_this);
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "5px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "18px")
            })
            .on("mouseout", function() {
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "3px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "12px")
            });
          line_svg
            .append("text")
            .attr("id", _n + "_" + _t + "_text")
            .attr("x", xScale(_v.length - 1) + x_offset + 15)
            .attr("y", yScale(_v[_v.length - 1]) + 5)
            .style("fill", "red")
            .text(_n_type)
            .on("mouseover", function() {
              let _this;
              _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
              _this.parentNode.appendChild(_this);
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "5px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "18px")
            })
            .on("mouseout", function() {
              d3.select("path#" + "line_" + _n + "_" + _t)
                .style("stroke-width", "3px")
              d3.select("text#" + _n + "_" + _t + "_text")
                .style("font-size", "12px")
            });
        }
        // add nodes and info on hover
        for (let _i_ in this.categories) {
          let _v_ = _v[_i_];
          let _s_ = _s[_i_];
          line_svg
            .append("path")
            .attr("id", function() {
              return _n + "_" + _i_
            })
            .style("fill", "white")
            .style("stroke", "black")
            .style("stroke-width", function() {
              if ((_s_ === undefined) || (_s_ === null)) {
                return nonsignificance_weight;
              } else if (_s_ < stat_value) {
                return significance_weight;
              } else {
                return nonsignificance_weight;
              }
            })
            .attr("d", d3.symbol()
              .size(function() {
                return 100;
              })
              .type(function() {
                if (_t === "gene_component" || _st === "gene_component") {
                  return d3.symbolTriangle;
                } else if (_t === "protein_component" || _st === "protein_component") {
                  return d3.symbolDiamond;
                } else if (_t === "metabolite_component" || _st === "metabolite_component") {
                  return d3.symbolCircle;
                } else if (_t === "complex_component" || _st === "complex_component") {
                  return d3.symbolSquare;
                } else {
                  return d3.symbolCross;
                }
              }))
            .attr("transform", "translate(" + (xScale(_i_) + x_offset) + "," + (yScale(_v_)) + ")")
        }
      }
    }
    for (let r in d.reactants) {
      let _i = this.nodes[d.reactants[r]];
      let _n = _i.name.replace(/,/g, '_').replace(/:/g, '_');
      let _n_type = _i.name + " (React.)";
      let _t = _i.type;
      let _st = _i.sub_type;
      let _v = _i.values;
      let _s = _i.stats;
      let _inf = _i.inferred;
      if (_v[0] != null && _t !== "complex_component" && _inf !== "true") {
        let dash_instruction;
        let _set = new Set(_v)
        if (_set.size === 1) {
          dash_instruction = "15 5"
        } else {
          dash_instruction = "1 0"
        }
        let current_line = d3.line()
          .x(function(d, i) {
            return xScale(i) + x_offset;
          })
          .y(function(d, i) {
            return yScale(d);
          })
          .curve(d3.curveMonotoneX);
        line_svg
          .append("path")
          .datum(_v)
          .attr("id", function() {
            return "line_" + _n + "_" + _t
          })
          .attr("class", "line")
          .style("fill", "none")
          .style("stroke-width", "3px")
          .style("stroke", "#808080")
          .style("stroke-dasharray", dash_instruction)
          .attr("d", current_line)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
              .style("z-index", "1")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        line_svg
          .append("text")
          .attr("id", _n + "_" + _t + "_text")
          .attr("x", xScale(0) + x_offset + 15)
          .attr("y", yScale(_v[0]) + 5)
          .style("fill", "#808080")
          .text(_n_type)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        // add nodes and info on hover
        for (let _i_ in this.categories) {
          let _v_ = _v[_i_];
          let _s_ = _s[_i_];
          line_svg
            .append("path")
            .attr("id", function() {
              return _n + "_" + _i_
            })
            .style("fill", "white")
            .style("stroke", "black")
            .style("stroke-width", function() {
              if ((_s_ === undefined) || (_s_ === null)) {
                return nonsignificance_weight;
              } else if (_s_ < stat_value) {
                return significance_weight;
              } else {
                return nonsignificance_weight;
              }
            })
            .attr("d", d3.symbol()
              .size(function() {
                return 100;
              })
              .type(function() {
                if (_t === "gene_component" || _st === "gene_component") {
                  return d3.symbolTriangle;
                } else if (_t === "protein_component" || _st === "protein_component") {
                  return d3.symbolDiamond;
                } else if (_t === "metabolite_component" || _st === "metabolite_component") {
                  return d3.symbolCircle;
                } else if (_t === "complex_component" || _st === "complex_component") {
                  return d3.symbolSquare;
                } else {
                  return d3.symbolCross;
                }
              }))
            .attr("transform", "translate(" + (xScale(_i_) + x_offset) + "," + (yScale(_v_)) + ")")
        }
      }
    }
    for (let p in d.products) {
      let _i = this.nodes[d.products[p]];
      let _n = _i.name.replace(/,/g, '_').replace(/:/g, '_');
      let _n_type = _i.name + " (Prod.)";
      let _t = _i.type;
      let _st = _i.sub_type;
      let _v = _i.values;
      let _s = _i.stats;
      let _inf = _i.inferred;
      let _v_ = [_v, _s];
      if (_v[0] != null && _t !== "complex_component" && _inf !== "true") {
        let dash_instruction;
        let _set = new Set(_v)
        if (_set.size === 1) {
          dash_instruction = "15 5"
        } else {
          dash_instruction = "1 0"
        }
        let current_line = d3.line()
          .x(function(d, i) {
            return xScale(i) + x_offset;
          })
          .y(function(d, i) {
            return yScale(d);
          })
          .curve(d3.curveMonotoneX);
        line_svg
          .append("path")
          .datum(_v)
          .attr("id", function() {
            return "line_" + _n + "_" + _t
          })
          .attr("class", "line")
          .style("fill", "none")
          .style("stroke-width", "3px")
          .style("stroke", "#808080")
          .style("stroke-dasharray", dash_instruction)
          .attr("d", current_line)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        line_svg
          .append("text")
          .attr("id", _n + "_" + _t + "_text")
          .attr("x", xScale(0) + x_offset + 15)
          .attr("y", yScale(_v[0]) + 5)
          .style("fill", "#808080")
          .text(_n_type)
          .on("mouseover", function() {
            let _this;
            _this = d3.select("text#" + _n + "_" + _t + "_text")['_groups'][0][0];
            _this.parentNode.appendChild(_this);
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "5px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "18px")
          })
          .on("mouseout", function() {
            d3.select("path#" + "line_" + _n + "_" + _t)
              .style("stroke-width", "3px")
            d3.select("text#" + _n + "_" + _t + "_text")
              .style("font-size", "12px")
          });
        // add nodes and info on hover
        for (let _i_ in this.categories) {
          let _v_ = _v[_i_];
          let _s_ = _s[_i_];
          line_svg
            .append("path")
            .attr("id", function() {
              return _n + "_" + _i_
            })
            .style("fill", "white")
            .style("stroke", "black")
            .style("stroke-width", function() {
              if ((_s_ === undefined) || (_s_ === null)) {
                return nonsignificance_weight;
              } else if (_s_ < stat_value) {
                return significance_weight;
              } else {
                return nonsignificance_weight;
              }
            })
            .attr("d", d3.symbol()
              .size(function() {
                return 100;
              })
              .type(function() {
                if (_t === "gene_component" || _st === "gene_component") {
                  return d3.symbolTriangle;
                } else if (_t === "protein_component" || _st === "protein_component") {
                  return d3.symbolDiamond;
                } else if (_t === "metabolite_component" || _st === "metabolite_component") {
                  return d3.symbolCircle;
                } else if (_t === "complex_component" || _st === "complex_component") {
                  return d3.symbolSquare;
                } else {
                  return d3.symbolCross;
                }
              }))
            .attr("transform", "translate(" + (xScale(_i_) + x_offset) + "," + (yScale(_v_)) + ")")
        }
      }
    }
  }
}

function reset_dot() {
  if (timecourse === true) {
    d3.select("circle#dot")
      .attr("cx", 90)
  }
}

function reset_filter() {
  document.getElementById("pathwayMenu-motif").value = "No metabolite co-factor selection...";
}

function reset_objects() {
  d3.select("#motif-pathway-svg")
    .style("visibility", "hidden");
  d3.select("#pathway-view-svg")
    .style("visibility", "hidden");
  d3.select(".network-panel")
    .style("visibility", "hidden");
  if (timecourse === true) {
    d3.select("svg#line-graph").remove();
  }
}

function reset_all() {
  d3.selectAll("circle:not(#dot)").remove();
  d3.selectAll("line:not(#track)").remove();
  d3.selectAll("text:not(#tick)").remove();
}

function highlight_selection(_selector) {

  if (timecourse === true) {
    d3.select("svg#line-graph").remove();
  }
  document.getElementById("pathway_name").innerHTML = "<h6><b> </b></h6>"

  let _selectors = [
    "#avg_num", //avg_num
    "#maxmax_num", //maxmax_num
    "#maxmin_num", //maxmin_num
    "#minmax_num", //minmax_num
    "#minmin_num", //minmin_num
    "#sustained_num", //sustained_num
    "#modreg_num",
    "#transreg_num",
    "#enzyme_num",
    "#activity_num"
  ]

  for (s in _selectors) {
    d3.select(_selectors[s])
      .style("background-color", "white");
  }

  d3.select(_selector)
    .style("background-color", "#FF7F7F");


}

function get_eval_motifs(eval_collapsed, motifs, indexer) {
  let eval_motifs;
  if (eval_collapsed === false) {
    eval_motifs = [];
    for (let m in motifs[indexer]) {
      if (motifs[indexer][m].collapsed === false || motifs[indexer][m].collapsed === "false") {
        eval_motifs.push(motifs[indexer][m]);
      }
    }
  } else {
    eval_motifs = motifs[indexer];
  }
  return eval_motifs;
}

function get_included_motifs(current_motifs, motifs, exclusion_indexer) {
  let motif_list = [];
  let subtracting_motifs = motifs[exclusion_indexer];
  let subtracting_ids = [];
  if (subtracting_motifs !== undefined) {
    for (let s in subtracting_motifs) {
      subtracting_ids.push(subtracting_motifs[s].id);
    }
    for (let m in current_motifs) {
      if (!subtracting_ids.includes(current_motifs[m].id)) {
        motif_list.push(current_motifs[m]);
      }
    }
  } else {
    motif_list = current_motifs;
  }
  return motif_list;
}

function sort_motifs(motif_list, motif_significance, sort_type) {
  if (sort_type === "Sort Number of Pathways") {
    motif_list.sort(function(a, b) {
      return d3.descending(a.pathways.length, b.pathways.length);
    })
  } else if (sort_type === "Sort Magnitude Change") {
    motif_list.sort(function(a, b) {
      return d3.descending(a.magnitude_change, b.magnitude_change);
    })
  } else if (sort_type === "Sort Statistical Significance") {
    motif_list.forEach(m => {
      if (m.p_values === undefined) {
        //pass
      } else if (m.p_values.source <= stat_value && m.p_values.target <= stat_value) {
        m.significance_type = 'Both significant';
        motif_significance.both.push(m);
      } else if (m.p_values.source <= stat_value) {
        m.significance_type = 'Source significant';
        motif_significance.one.push(m);
      } else if (m.p_values.target <= stat_value) {
        m.significance_type = 'Target significant';
        motif_significance.one.push(m);
      } else { // both > stat_value
        m.significance_type = 'Both not significant';
        motif_significance.none.push(m);
      }
    })
    motif_significance.both.sort(function(a, b) {
      return d3.ascending(a.p_values.source * a.p_values.target, b.p_values.source * b.p_values.target);
    })
    motif_significance.one.sort(function(a, b) {
      return d3.ascending(a.p_values.source * a.p_values.target, b.p_values.source * b.p_values.target);
    })
    motif_significance.none.sort(function(a, b) {
      return d3.ascending(a.p_values.source * a.p_values.target, b.p_values.source * b.p_values.target);
    })
    motif_list = motif_significance.both.concat(
      motif_significance.one,
      motif_significance.none);
  } else if (sort_type === "Sort Reaction FDR") {
    motif_list.sort(function(a, b) {
      return d3.ascending(a.p_values.agg, b.p_values.agg);
    })
  } return [motif_list, motif_significance];
}

function remove_duplicate_motifs(motif_list) {
  let checked_rxns = {};
  let remove_idx = [];
  for (let i = 0; i < motif_list.length - 1; i++) {
    checked_rxns[motif_list[i].id] = motif_list[i].compartment;
    let i2 = parseInt(i) + 1;
    if (checked_rxns[motif_list[i2].id] === motif_list[i2].compartment) {
      remove_idx.push(i2);
    }
  }
  for (var i = remove_idx.length -1; i >= 0; i--) {
    motif_list.splice(remove_idx[i], 1);
  }
  return motif_list;
}

function remove_duplicate_motifs_twoReactions(motif_list) {
  let checked_rxns = {};
  let remove_idx = [];
  for (let i = 0; i < motif_list.length - 1; i++) {
    checked_rxns[motif_list[i].rxn1.id] = motif_list[i].rxn2.id;
    let i2 = parseInt(i) + 1;
    if (checked_rxns[motif_list[i2].rxn1.id] === motif_list[i2].rxn2.id
    || checked_rxns[motif_list[i2].rxn2.id] === motif_list[i2].rxn1.id
    || motif_list[i2].rxn1.id === motif_list[i2].rxn2.id) {
      remove_idx.push(i2);
    }
  }
  for (var i = remove_idx.length -1; i >= 0; i--) {
    motif_list.splice(remove_idx[i], 1);
  }
  return motif_list;
}

function remove_noshare_twoReactions(motif_list) {
  let remove_idx = [];
  for (let i = 0; i < motif_list.length - 1; i++) {
    let components_rxn1 = motif_list[i].rxn1.reactants.concat(
      motif_list[i].rxn1.products,
      motif_list[i].rxn1.modifiers.map(x => x[0])
    );
    let components_rxn2 = motif_list[i].rxn2.reactants.concat(
      motif_list[i].rxn2.products,
      motif_list[i].rxn2.modifiers.map(x => x[0])
    );
    let shared_components = components_rxn1.filter(
      value => components_rxn2.includes(value));
    if (shared_components.length === 0) {
      remove_idx.push(i);
    }
  }
  for (var i = remove_idx.length -1; i >= 0; i--) {
    motif_list.splice(remove_idx[i], 1);
  }
  return motif_list;
}

function get_reaction_dict(data) {
  if (eval_collapsed === false) {
    return data.reaction_dict;
  } else {
    return data.collapsed_reaction_dict;
  }
}