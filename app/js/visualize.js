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

let current_path;
let database_url;
var selection = null;
var superSelection = null;
var selector = "#graph";
var _width = window.innerWidth;
var _height = window.innerHeight - 75;

// MAIN
function set_visualization(database_url) {
  try {
    var data = JSON.parse(fs.readFileSync(database_url).toString());
  } catch (e) {
    alert('Failed to open: \n' + database_url)
  }
  console.log(data)
  var stat_type = data.metadata.stat_type;
  set_stat_button(stat_type);

  var pathway_dict = make_pathway_dictionary(
    data,
    'pathway_dictionary');
  var collapsed_pathway_dict = make_pathway_dictionary(
    data,
    'collapsed_pathway_dictionary');
  var superPathwayDict = make_superPathway_dictionary(data);
  var entity_dictionary = parseEntities(data.nodes);

  var motif_outputs = gatherMotifs(data, data.categories);
  var collapsed_global_motifs = motif_outputs[0];
  var global_motifs = motif_outputs[1];

  let names;
  if (data.labels === null) {
    names = [];
  } else {
    names = data.labels.split(',');
  }


  timecourse = checkCategories(data.categories, names); //, data.names);

  make_menu(
    superPathwayDict,
    "superPathwayMenu",
    "Select a super pathway...",
    (provide_all = true)
  );

  var update_nodes = {};
  for (n in data.nodes) {
    update_nodes[data.nodes[n].id] = data.nodes[n];
  }
  data.nodes = update_nodes;

  var update_links = {};
  for (l in data.links) {
    let link_id = data.links[l].source + "," + data.links[l].target;
    update_links[link_id] = data.links[l];
  }
  data.links = update_links;

  data.blocklist = data.species_blocklist;
  data.blocklist = complete_blocklist(
    data.blocklist,
    data.metadata.blocklist,
    data.nodes
  )

  d3.select("#superPathwayMenu").on("change", function() {
    changeSuper(data, pathway_dict, superPathwayDict);
  });
  d3.select("#pathwayMenu").on("change", function() {
    change(data, collapsed_pathway_dict, pathway_dict);
  });
  d3.select("#kNN_button").on("change", function() {
    kNN_input(data, collapsed_pathway_dict, pathway_dict)
  });
  d3.select("#hub_button").on("change", function() {
    hub_input(data, collapsed_pathway_dict, pathway_dict)
  });
  d3.select("#stat_button").on("change", function() {
    stat_input(data, collapsed_pathway_dict, pathway_dict)
  });
}


async function main() {
  try {
    current_path = await get_session_info_async("current_pathway")
    if ((current_path != null) & (current_path != "null")) {
      // set back after opening
      update_session_info("current_pathway", null);
    } else {}

    database_url = await get_session_info_async("database_url");
    console.log("Database path: " + database_url);
    set_visualization(database_url);
  } catch (err) {
    console.log(err);
  }
}


main();
