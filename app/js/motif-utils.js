/*
Metaboverse
Visualizing and Analyzing Metabolic Networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

MIT License

Copyright (c) 2022 Metaboverse

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

function update_nodes_links(_nodes, _links) {
  var update_nodes = {};
  let n;
  for (n in _nodes) {
    update_nodes[_nodes[n].id] = _nodes[n];
  }

  var update_links = {};
  let l;
  for (l in _links) {
    let link_id = _links[l].source + "," + _links[l].target;
    update_links[link_id] = _links[l];
  }
  
  return [update_nodes, update_links];
}

function complete_blocklist(blocklist, blocklist_names, _nodes) {

  let split_blocklist = blocklist_names.split(',');
  let add_blocklist = [];

  for (let n in _nodes) {
    if (split_blocklist.includes(_nodes[n].name)) {
      add_blocklist.push(_nodes[n].id);
    }
  }
  
  return [... new Set(blocklist.concat(add_blocklist))];
}

function create_dictionaries(_nodes) {
  
  let expression_dict = {};
  let stats_dict = {};
  let inferred_dict = {};
  for (let x in _nodes) {
    let id = _nodes[x]['id'];
    let expression = _nodes[x]['values'];
    let stats = _nodes[x]['stats'];
    expression_dict[id] = expression;
    stats_dict[id] = stats;
    inferred_dict[id] = _nodes[x]['inferred']
  }

  return [expression_dict, stats_dict, inferred_dict];
}

function create_link_neighbors(_nodes, _links) {

  let link_neighbors = {};
  for (let l in _links) {
    let _source = _links[l].source;
    let _target = _links[l].target;
    
    if (_nodes[_source].type === "reaction" 
    || _nodes[_target].type === "reaction"
    || _nodes[_source].type === "collapsed"
    || _nodes[_target].type === "collapsed") {
    } else {
      if (!(_source in link_neighbors)) {
        link_neighbors[_source] = [];
      }
      link_neighbors[_source].push(_target);
  
      if (!(_target in link_neighbors)) {
        link_neighbors[_target] = [];
      }
      link_neighbors[_target].push(_source);
    }
  }

  return link_neighbors;
}

function make_metabolite_species_dictionary(data) {

  let species_dictionary = {};

  for (let n in data.nodes) {
    if (data.nodes[n].type === "metabolite_component" || data.nodes[n].sub_type === "metabolite_component") {
      let this_name;
      if (data.nodes[n].user_label !== undefined) {
        this_name = data.nodes[n].user_label;
      } else {
        this_name = data.nodes[n].name;
      }
      if (this_name in species_dictionary) {
        species_dictionary[this_name].push(data.nodes[n].id);
      } else {
        species_dictionary[this_name] = [data.nodes[n].id];
      }
    }
  }

  return species_dictionary;
}

function make_entity_species_r_dictionary(data) {

  let species_dictionary = {};

  for (let n in data.nodes) {
    let this_name;
    if (data.nodes[n].user_label !== undefined) {
      this_name = data.nodes[n].user_label;
    } else {
      this_name = data.nodes[n].name;
    }
    if (data.nodes[n].id in species_dictionary) {
      species_dictionary[data.nodes[n].id].push(this_name);
    } else {
      species_dictionary[data.nodes[n].id] = [this_name];
    }
  }

  return species_dictionary;
}