/*
Metaboverse
Metaboverse is designed for analysis of metabolic networks
https://github.com/Metaboverse/Metaboverse/
alias: metaboverse

Copyright (C) 2019  Youjia Zhou, Jordan A. Berg
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

var sample = 0;

database_url = get_session_info("database_url");
console.log("Database path: " + database_url);

class MetaGraph{
  constructor(graphdata){

    // Get the data
    console.log(graphdata)
    this.data = graphdata;
    this.nodes = graphdata.nodes;
    this.collapsed_reaction_dict = graphdata.collapsed_reaction_dictionary;
    this.mod_collapsed_pathways = graphdata.mod_collapsed_pathways;
    this.collapsed_pathway_dict = make_pathway_dictionary(
      graphdata,
      "collapsed_pathway_dictionary"
    )
    this.path_mapper = graphdata.motif_reaction_dictionary

    this.node_dict = {};
    for (let node in this.nodes) {
      this.node_dict[this.nodes[node]['id']] = this.nodes[node];
    }

    // Generate expression and stats dictionaries
    let expression_dict = {};
    let stats_dict = {};
    for (let x in this.nodes) {
      let id = this.nodes[x]['id'];
      let expression = this.nodes[x]['values'];
      let stats = this.nodes[x]['stats'];

      expression_dict[id] = expression;
      stats_dict[id] = stats;
    }
    this.expression_dict = expression_dict;
    this.stats_dict = stats_dict;

    // Generate stamp view output
    try {
    this.stamp_svg = d3.select("#stamp-view-svg");
    this.stamp_svg_width = parseFloat(this.stamp_svg.style("width"));
    this.stamp_svg_height = parseFloat(this.stamp_svg.style("height"));
    this.stamp_svg_margin = {
      "top": 5,
      "left": 5,
      "horizontal": 10,
      "vertical": 10};
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
      .attr("id","mp-pathway-group");
    this.mp_selection_group = this.mp_svg.append("g")
      .attr("id","mp-selection-group");

    // Generate pathway viewer
    // Just get pathway ID and let the viz script do the rest
    this.pathway_svg = d3.select("#pathway-view-svg");
    this.pathway_svg_width = parseFloat(this.pathway_svg.style("width", "45vw"));
    this.pathway_svg_height = parseFloat(this.pathway_svg.style("height", "490px"));
  } catch(e) {}

    this.motifSearch();
    }

    motifSearch(){
      d3.select("#motif1")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg")
            .style("visibility", "hidden");
          d3.select("#pathway-view-svg")
            .style("visibility", "hidden");
          d3.select(".network-panel")
            .style("visibility", "hidden");
          let threshold = d3.select("#avg_num").node().value;
          this.motif = motifSearch_Avg(
            threshold,
            this.collapsed_reaction_dict,
            this.expression_dict,
            this.path_mapper)
          if (this.motif !== undefined) {
            this.drawMotifSearchResult(this.motif);
          }
        })

      d3.select("#motif2")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg")
            .style("visibility","hidden");
          d3.select("#pathway-view-svg")
            .style("visibility","hidden");
          d3.select(".network-panel")
            .style("visibility","hidden");
          let threshold = d3.select("#maxmax_num").node().value;
          this.motif = motifSearch_MaxMax(
            threshold,
            this.collapsed_reaction_dict,
            this.expression_dict,
            this.path_mapper)
          if (this.motif !== undefined) {
            this.drawMotifSearchResult(this.motif);
          }
        })

      d3.select("#motif3")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg")
            .style("visibility","hidden");
          d3.select("#pathway-view-svg")
            .style("visibility","hidden");
          d3.select(".network-panel")
            .style("visibility","hidden");
          let threshold = d3.select("#maxmin_num").node().value;
          this.motif = motifSearch_MaxMin(
            threshold,
            this.collapsed_reaction_dict,
            this.expression_dict,
            this.path_mapper)
          if (this.motif !== undefined) {
            this.drawMotifSearchResult(this.motif);
          }
        }
      )

      d3.select("#motif4")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg")
            .style("visibility","hidden");
          d3.select("#pathway-view-svg")
            .style("visibility","hidden");
          d3.select(".network-panel")
            .style("visibility","hidden");
          let threshold = d3.select("#sustained_num").node().value;
          this.motif = motifSearch_Sustained(
            threshold,
            this.collapsed_reaction_dict,
            this.expression_dict,
            this.path_mapper)
          if (this.motif !== undefined) {
            this.drawMotifSearchResult(this.motif);
          }
        }
      )

      d3.select("#motif5")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg")
            .style("visibility","hidden");
          d3.select("#pathway-view-svg")
            .style("visibility","hidden");
          d3.select(".network-panel")
            .style("visibility","hidden");
          let threshold = d3.select("#pathmax_num").node().value;
          this.motif = motifSearch_PathMax(
            threshold,
            this.mod_collapsed_pathways,
            this.collapsed_reaction_dict,
            this.expression_dict,
            this.path_mapper)

          if (this.motif !== undefined) {
            this.drawMotifSearchResultPathway(this.motif);
          }
        }
      )
    }

    drawMotifSearchResult(motif_list){
      //Change this to sort by p-values
      //motif_list.sort(function(a,b){
      //  console.log(a.pathways)
      //  return d3.descending(a.pathways.length, b.pathways.length)
      //})
      let stamp_height = 50;
      this.stamp_svg_height = Math.ceil(
        motif_list.length / 3)
        * (stamp_height + this.stamp_svg_margin.vertical);
      this.stamp_svg.attr("height",this.stamp_svg_height);

      let stamp_width = this.stamp_svg_width / 3
        - this.stamp_svg_margin.horizontal;

      let sg = this.stamp_svg_selection_group.selectAll("rect")
        .data(motif_list);
      sg.exit().remove();
      sg = sg.enter().append("rect").merge(sg)
        .attr("x",(d,i)=>this.stamp_svg_margin.left + i%3*(stamp_width+this.stamp_svg_margin.horizontal))
        .attr("y",(d,i)=>this.stamp_svg_margin.top + Math.floor(i/3)*(stamp_height+this.stamp_svg_margin.vertical))
        .attr("width",stamp_width)
        .attr("height",stamp_height)
        .attr("fill","lightgray")
        .attr("id",(d)=>"stamp-cover-"+d.id)
        .style("opacity",0)
        .on("click",(d)=>{
          this.drawMotifPathway(d);
          d3.select("#motif-pathway-svg").style("visibility","visible");
        })
        .on("mouseover",(d)=>{
          d3.select("#stamp-cover-"+d.id).style("opacity",0.5);
        })
        .on("mouseout",(d)=>{
          d3.select("#stamp-cover-"+d.id).style("opacity",0);
        });

      let fg = this.stamp_svg_frame_group.selectAll("rect")
        .data(motif_list);
      fg.exit().remove();
      fg = fg.enter().append("rect").merge(fg)
        .attr("x",(d,i)=>this.stamp_svg_margin.left + i % 3
          * (stamp_width + this.stamp_svg_margin.horizontal))
        .attr("y",(d,i)=>this.stamp_svg_margin.top + Math.floor(i / 3 )
          * (stamp_height + this.stamp_svg_margin.vertical))
        .attr("width",stamp_width)
        .attr("height",stamp_height)
        .attr("stroke","lightgray")
        .attr("fill","white")

      let cg = this.stamp_svg_circle_group.selectAll("g")
        .data(motif_list);
      cg.exit().remove();
      cg = cg.enter().append("g").merge(cg)
        .attr("id",(d,i)=>"stamp-circle-"+i);

      let lg = this.stamp_svg_link_group.selectAll("g")
        .data(motif_list);
      lg.exit().remove();
      lg = lg.enter().append("g").merge(lg)
        .attr("id",(d,i)=>"stamp-link-"+i);

      for(let i=0; i<motif_list.length; i++){
        let mnodes = [];
        let mnodes_id = [];
        let mlinks = [];
        let r_idx = 0;
        let p_idx = 0;
        // Add reaction node
        this.node_dict[motif_list[i].id].current_type = "reaction";
        mnodes.push(this.node_dict[motif_list[i].id]);

        // Add reactant nodes
        motif_list[i].reactants.forEach(l=>{
          if(mnodes_id.indexOf(l)===-1){
            this.node_dict[l].current_type = "reactant";
            this.node_dict[l].r_idx = r_idx;
            mnodes.push(this.node_dict[l]);
            mnodes_id.push(l);
            mlinks.push({'source': l, 'target': motif_list[i].id});
            r_idx += 1;
          }
        })

        // Add product nodes
        motif_list[i].products.forEach(l=>{
          if(mnodes_id.indexOf(l)===-1){
            this.node_dict[l].current_type = "product";
            this.node_dict[l].p_idx = p_idx;
            mnodes.push(this.node_dict[l]);
            mnodes_id.push(l);
            mlinks.push({'source': motif_list[i].id, 'target': l});
            p_idx += 1;
          }
        })

        let start_x = this.stamp_svg_margin.left + i%3*(stamp_width+this.stamp_svg_margin.horizontal);
        let start_y = this.stamp_svg_margin.top + Math.floor(i/3)*(stamp_height+this.stamp_svg_margin.vertical);
        let x_interval = stamp_width/3;
        let y_interval_r = (stamp_height-20)/r_idx;
        let y_interval_p = (stamp_height-20)/p_idx;

        let mg = d3.select("#stamp-circle-"+i).selectAll("circle")
          .data(mnodes);
        mg.exit().remove();
        mg = mg.enter().append("circle").merge(mg)
          .attr("cx",(d)=>{
            if(d.current_type==="reactant"){
              return start_x + 15;
            } else if(d.current_type==="reaction"){
              return start_x + x_interval+15;
            } else if (d.current_type ==="product"){
              return start_x + x_interval*2+15;
            }
          })
          .attr("cy",(d)=>{
            if(d.current_type==="reactant"){
              return start_y + y_interval_r*(d.r_idx)+10;
            } else if(d.current_type==="reaction"){
              return start_y + 20;
            } else if (d.current_type ==="product"){
              return start_y + y_interval_p*(d.p_idx)+10;
            }
          })
          .attr("class",(d)=>d.current_type)
          .attr("fill", (d)=>{
            if (d.values_js === undefined) {
              return "rgba(191, 191, 191, 1)";
            } else {
              return "rgba(" + d.values_js.toString() + ")";
            }
          })
          .attr("r",4)
          .attr("stroke","black")
          .attr("id",d=>"stamp-"+i+"-"+d.id)

        let mlg = d3.select("#stamp-link-"+i).selectAll("line")
          .data(mlinks);
        mlg.exit().remove();
        mlg = mlg.enter().append("line").merge(mlg)
          .attr("x1",(d)=>d3.select("#stamp-"+i+"-"+d.source).attr("cx"))
          .attr("y1",(d)=>d3.select("#stamp-"+i+"-"+d.source).attr("cy"))
          .attr("x2",(d)=>d3.select("#stamp-"+i+"-"+d.target).attr("cx"))
          .attr("y2",(d)=>d3.select("#stamp-"+i+"-"+d.target).attr("cy"))
          .attr("stroke","gray")

        let tg = d3.select("#stamp-circle-"+i).selectAll("text")
          .data([motif_list[i]]);
        tg.exit().remove();
        tg = tg.enter().append("text").merge(tg)
          .attr("x",start_x + 10)
          .attr("y",start_y + 45)
          .text(d=>d.pathways.length + " pathways")
          .style("font-size","11px")
      }
    }

    drawMotifSearchResultPathway(motif_list){
      //Change this to sort by p-values
      //motif_list.sort(function(a,b){
      //  console.log(a.pathways)
      //  return d3.descending(a.pathways.length, b.pathways.length)
      //})

      let stamp_height = 30;
      this.stamp_svg_height = Math.ceil(
        motif_list.length)
        * (stamp_height + this.stamp_svg_margin.vertical);
      this.stamp_svg.attr("height",this.stamp_svg_height);

      let stamp_width = this.stamp_svg_width
        - this.stamp_svg_margin.horizontal - 5;

      let sg = this.stamp_svg_selection_group.selectAll("rect")
        .data(motif_list);
        sg.exit().remove();
        sg = sg.enter().append("rect").merge(sg)
          .attr("x",5)
          .attr("y",(d,i)=>this.stamp_svg_margin.top + Math.floor(i)*(stamp_height+this.stamp_svg_margin.vertical))
          .attr("width",stamp_width)
          .attr("height",stamp_height)
          .attr("fill","lightgray")
          .attr("id",(d)=>"stamp-cover-"+d.id)
          .style("opacity",0)
          .on("click",(d)=>{
            document.getElementById("pathway_name").innerHTML = "<h6><b>" + d.name + "</b></h6>" ;
            this.drawPathwayView(d.id, "#pathway-view-svg");
            d3.select("#pathway-view-svg").style("visibility","visible");
            d3.select(".network-panel").style("visibility","visible");
          })
          .on("mouseover",(d)=>{
            d3.select("#stamp-cover-" + d.id).style("opacity",0.4);
          })
          .on("mouseout",(d)=>{
            d3.select("#stamp-cover-" + d.id).style("opacity",0);
          })

      let fg = this.stamp_svg_frame_group.selectAll("rect")
        .data(motif_list);
        fg.exit().remove();
        fg = fg.enter().append("rect").merge(fg)
          .attr("x",5)
          .attr("y",(d,i)=>this.stamp_svg_margin.top + Math.floor(i)
            * (stamp_height + this.stamp_svg_margin.vertical))
          .attr("width",stamp_width)
          .attr("height",stamp_height)
          .attr("stroke","lightgray")
          .attr("fill","white")

      d3.selectAll("circle").remove();
      d3.selectAll("line").remove();
      d3.selectAll("text").remove();

      let cg = this.stamp_svg_circle_group.selectAll("g")
        .data(motif_list);
      cg.exit().remove();
      cg = cg.enter().append("g").merge(cg)
        .attr("id",(d,i)=>"stamp-circle-"+i);
      console.log(cg)
      for(let i=0; i<motif_list.length; i++){

        let start_x = 5;
        let start_y = this.stamp_svg_margin.top + Math.floor(i)*(stamp_height+this.stamp_svg_margin.vertical);

        let tg = d3.select("#stamp-circle-"+i).selectAll("text")
          .data([motif_list[i]]);
        tg.exit().remove();
        tg = tg.enter().append("text").merge(tg)
          .attr("x",start_x + 10)
          .attr("y",start_y + 21)
          .text(d=> {
            if (d.name.length < 48) {
              return d.name;
            } else {
              return d.name.substring(0,45) + " ...";
            }
          })
          .style("font-size","14px")
          .style("font-weight","bold")
      }
    }

    drawMotifPathway(motif){
      // recover graph
      let mnodes = [motif];
      let mlinks = [];
      let mnodes_id = [];
      let r_idx = 0;
      let p_idx = 0;

      motif.current_type = "reaction";
      motif.reactants.forEach(l=>{
        if(mnodes_id.indexOf(l)===-1){
          this.node_dict[l].current_type = "reactant";
          this.node_dict[l].r_idx = r_idx;
          mnodes.push(this.node_dict[l]);
          mnodes_id.push(l);
          mlinks.push({'source': l, 'target': motif.id});
          r_idx += 1;
        }
      })
      motif.products.forEach(l=>{
        if(mnodes_id.indexOf(l)===-1){
          this.node_dict[l].current_type = "product";
          this.node_dict[l].p_idx = p_idx;
          mnodes.push(this.node_dict[l]);
          mnodes_id.push(l);
          mlinks.push({'source': motif.id, 'target': l});
          p_idx += 1;
        }
      })

      let motif_height = 120;

      // ***** draw motif glyph *****
      let x_interval = this.mp_svg_width/3;
      let y_interval_r = (motif_height-20)/r_idx;
      let y_interval_p = (motif_height-20)/p_idx;
      let ng = this.mp_motif_circle_group.selectAll("circle")
        .data(mnodes)
      ng.exit().remove();
      ng = ng.enter()
        .append("circle").merge(ng)
          .attr("cx",(d)=>{
            if(d.current_type==="reactant"){
              return 75;
            } else if(d.current_type==="reaction"){
              return x_interval+75;
            } else if (d.current_type ==="product"){
              return x_interval*2+75;
            }
          })
          .attr("cy",(d)=>{
            if(d.current_type==="reactant"){
              return y_interval_r*(d.r_idx)+20;
            } else if(d.current_type==="reaction"){
              return motif_height/2;
            } else if (d.current_type ==="product"){
              return y_interval_p*(d.p_idx)+20;
            }
          })
          .attr("class",(d)=>d.current_type)
          .attr("fill", (d)=>{
            if (d.values_js === undefined) {
              return "rgba(191, 191, 191, 1)";
            } else {
              return "rgba(" + d.values_js.toString() + ")";
            }
          })
          .attr("r",8)
          .attr("stroke","black")
          .attr("id",d=>"mp-circle-"+d.id)

      let lg = this.mp_motif_link_group.selectAll("line")
        .data(mlinks);
      lg.exit().remove();
      lg = lg.enter().append("line").merge(lg)
        .attr("x1",(d)=>d3.select("#mp-circle-"+d.source).attr("cx"))
        .attr("y1",(d)=>d3.select("#mp-circle-"+d.source).attr("cy"))
        .attr("x2",(d)=>d3.select("#mp-circle-"+d.target).attr("cx"))
        .attr("y2",(d)=>d3.select("#mp-circle-"+d.target).attr("cy"))
        .attr("stroke","gray")
        .attr("marker-end", "url(#end)");

      let tg = this.mp_motif_circle_group.selectAll("text")
        .data([motif]);
      tg.exit().remove();
      tg = tg.enter().append("text").merge(tg)
        .attr("x",this.mp_svg_width/12 - 10)
        .attr("y",(motif_height-5))
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
            return d.name.substring(0,string_length) + " ..." + string_suffix;
          }
        })
        .style("font-size","12px")
        .style("font-weight","bold")

    // **** draw pathway glyph ****
    let pathway_list = motif.pathways;
    let pathway_height = 20;
    let margin = {"horizontal":9, "vertical":10, "top":10, "left":5};
    this.mp_svg_height = Math.ceil(pathway_list.length/3) * (pathway_height+margin.vertical) + motif_height + margin.top;
    this.mp_svg.attr("height",this.mp_svg_height);

    let pathway_width = (this.mp_svg_width / 3) - margin.horizontal;

    let sg = this.mp_selection_group.selectAll("rect")
      .data(pathway_list);
    sg.exit().remove();
    sg = sg.enter().append("rect").merge(sg)
      .attr("x",(d,i)=>margin.left + i%3*(pathway_width+margin.horizontal))
      .attr("y",(d,i)=>motif_height+margin.top + Math.floor(i/3)*(pathway_height+margin.vertical))
      .attr("width",pathway_width)
      .attr("height",pathway_height)
      .attr("fill","blue")
      .attr("id",(d)=>"mp-cover-"+d)
      .style("opacity",0)
      .on("click",(d)=>{
        this.drawPathwayView(d, "#pathway-view-svg");
        this.findAllMotif(d, this.motif);
        d3.select("#pathway-view-svg").style("visibility","visible");
        d3.select(".network-panel").style("visibility","visible");

      })
      .on("mouseover",(d)=>{
        d3.select("#mp-cover-" + d).style("opacity",0.4);
      })
      .on("mouseout",(d)=>{
        d3.select("#mp-cover-" + d).style("opacity",0);
      })

    let fg = this.mp_pathway_group.selectAll("rect")
      .data(pathway_list);
    fg.exit().remove();
    fg = fg.enter().append("rect").merge(fg)
      .attr("x",(d,i)=>margin.left + i%3*(pathway_width+margin.horizontal))
      .attr("y",(d,i)=>motif_height+margin.top + Math.floor(i/3)*(pathway_height+margin.vertical))
      .attr("width",pathway_width)
      .attr("height",pathway_height)
      .attr("fill","lightblue")
      .style("opacity",0.5)

    // get pathway names
    let ptg = this.mp_pathway_group.selectAll("text")
      .data(pathway_list);
    ptg.exit().remove();
    ptg = ptg.enter().append("text").merge(ptg)
      .attr("x",(d,i)=>margin.left + i%3*(pathway_width+margin.horizontal)+10)
      .attr("y",(d,i)=>motif_height+margin.top + Math.floor(i/3)*(pathway_height+margin.vertical)+12)
      .text(d=>this.mod_collapsed_pathways[d].name.substring(0,24))
      .style("font-size",9)

  }

  drawPathwayView(p, selector) {

    graph_genes = true;
    collapse_reactions = true;
    var motif_reactions = this.mod_collapsed_pathways[p]["reactions"];

    update_session_info("current_pathway", p);

    // Parse through each reaction listed and get the component parts
    let components = [];
    var rxn = 0;
    for (rxn in motif_reactions) {

      var target_rxns = this.collapsed_reaction_dict[motif_reactions[rxn]];
      components.push(motif_reactions[rxn]);
      for (x in target_rxns["reactants"]) {
        components.push(target_rxns["reactants"][x]);
      }
      for (x in target_rxns["products"]) {
        components.push(target_rxns["products"][x]);
      }
      for (x in target_rxns["modifiers"]) {
        components.push(target_rxns["modifiers"][x][0]);
      }
      for (x in target_rxns["additional_components"]) {
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

    var _width = 0.45 * (window.innerWidth);
    var _height = 570;

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
      _width,
      _height
    );
  }

  findAllMotif(pathway, motif_list){
    d3.select(".network-panel").style("visibility","visible");
    let current_pathway = this.mod_collapsed_pathways[pathway].name;
    document.getElementById("pathway_name").innerHTML = "<h6><b>" + current_pathway + "</b></h6>" ;

    let motif_id = [];
    motif_list.forEach(m=>{
      motif_id.push(m.id);
    })
    let pathway_list = this.mod_collapsed_pathways[pathway].reactions;
    let current_motif = [];
    pathway_list.forEach(react_id=>{
      if(motif_id.indexOf(react_id)!=-1){
        current_motif.push(react_id);
      }
    })
    current_motif.forEach(react_id=>{
      d3.selectAll("circle#" + react_id)
        .style("r", "16px")
        .style("stroke", "purple")
        .style("stroke-width", "5px")

    })
    this.showMotifNames(current_motif);
  }

  showMotifNames(current_motif) {
    let tg = d3.select("#all-motif-list")
      .select("ul")
      .selectAll("li")
      .data(current_motif);
    tg.exit().remove();
    tg = tg.enter().append("li").merge(tg)
      .html(d => "- " + this.collapsed_reaction_dict[d].name)
  }

  createId(id){
    return id.replace(/[^a-zA-Z0-9]/g, "")
  }
}

function runPathway() {
  console.log(current_pathway)
}
