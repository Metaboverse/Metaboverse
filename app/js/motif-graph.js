class MetaGraph{
  constructor(graphdata){

    // Get the data
    console.log(graphdata)
    this.react_nodes = graphdata.react_nodes;
    this.other_nodes = graphdata.other_nodes;
    this.pathways = graphdata.pathways_motif;

    // Generate stamp view output
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
    this.pathway_link_svg = this.pathway_svg.append("g")
      .attr("id", "pathway-link-group");
    this.pathway_circle_svg = this.pathway_svg.append("g")
      .attr("id", "pathway-circle-group");
    this.pathway_svg_width = parseFloat(this.pathway_svg.style("width"));
    this.pathway_svg_height = parseFloat(this.pathway_svg.style("height"));

    this.mp_svg.append("svg:defs").selectAll("marker")
        .data(["end"])      // Different link/path types can be defined here
      .enter().append("svg:marker")    // This section adds in the arrows
        .attr("id", String)
        .attr("refX", 15)
        .attr("refY", 0)
        .attr("markerWidth", 6)
        .attr("markerHeight", 6)
        .attr("orient", "auto")
      .append("svg:path")
        .attr("d", "M0,-5L10,0L0,5");

    this.pathway_svg.append("svg:defs").selectAll("marker")
        .data(["end"])      // Different link/path types can be defined here
      .enter().append("svg:marker")    // This section adds in the arrows
        .attr("id", String)
        // .attr("viewBox", "0 -5 10 10")
        .attr("refX", 15)
        .attr("refY", 0)
        .attr("markerWidth", 6)
        .attr("markerHeight", 6)
        .attr("orient", "auto")
      .append("svg:path")
        .attr("d", "M0,-5L10,0L0,5");

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
          this.motifSearch_Avg();
          this.drawMotifSearchResult(this.motif);
        })

      d3.select("#motif2")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg").style("visibility","hidden");
          d3.select("#pathway-view-svg").style("visibility","hidden");
          d3.select(".network-panel").style("visibility","hidden");
          this.motifSearch_MaxMax();
          this.drawMotifSearchResult(this.motif);
        })

      d3.select("#motif3")
        .on("click", ()=>{
          d3.select("#motif-pathway-svg").style("visibility","hidden");
          d3.select("#pathway-view-svg").style("visibility","hidden");
          d3.select(".network-panel").style("visibility","hidden");
          this.motifSearch_MaxMin();
          this.drawMotifSearchResult(this.motif);
        }
      )
    }

    motifSearch_Avg(){ // search for motif 1
      console.log("motif search 1")
      this.motif = [];
      let threshold = d3.select("#avg_num").node().value;
      for(let node_id in this.react_nodes){
        let react = this.react_nodes[node_id];
        let source_links = react.links.reactant;
        let target_links = react.links.product;
        let source_expression = [];
        let target_expression = [];
        source_links.forEach(l=>{
          let source_node = this.other_nodes[l.source];
          if(source_node.expression != "None"){
            source_expression.push(parseFloat(source_node.expression));
          }
        })
        target_links.forEach(l=>{
          let target_node = this.other_nodes[l.target];
          if(target_node.expression != "None"){
            target_expression.push(parseFloat(target_node.expression));
          }
        })
        if(source_expression.length>0 && target_expression.length>0){
          let source_avg = this.computeAvg(source_expression);
          let target_avg = this.computeAvg(target_expression);
          if(Math.abs(source_avg - target_avg)>=threshold){
            this.motif.push(react);
          }
        }
      }
      console.log(this.motif);
    }

    motifSearch_MaxMax(){ // MaxMax
      console.log("motif search 2")
      this.motif = [];
      let threshold = d3.select("#maxmax_num").node().value;
      for(let node_id in this.react_nodes){
        let react = this.react_nodes[node_id];
        let source_links = react.links.reactant;
        let target_links = react.links.product;
        let source_expression = [];
        let target_expression = [];
        source_links.forEach(l=>{
          let source_node = this.other_nodes[l.source];
          if(source_node.expression != "None"){
            source_expression.push(parseFloat(source_node.expression));
          }
        })
        target_links.forEach(l=>{
          let target_node = this.other_nodes[l.target];
          if(target_node.expression != "None"){
            target_expression.push(parseFloat(target_node.expression));
          }
        })
        if(source_expression.length>0 && target_expression.length>0){
          let source_max = Math.max(...source_expression);
          let target_max = Math.max(...target_expression);
          if(Math.abs(source_max - target_max)>=threshold){
            this.motif.push(react);
          }
        }
      }
    }

    motifSearch_MaxMin(){ //MaxMin
      console.log("motif search 3")
      this.motif = [];
      let threshold = d3.select("#maxmin_num").node().value;
      for(let node_id in this.react_nodes){
        let react = this.react_nodes[node_id];
        let source_links = react.links.reactant;
        let target_links = react.links.product;
        let source_expression = [];
        let target_expression = [];
        source_links.forEach(l=>{
          let source_node = this.other_nodes[l.source];
          if(source_node.expression != "None"){
            source_expression.push(parseFloat(source_node.expression));
          }
        })
        target_links.forEach(l=>{
          let target_node = this.other_nodes[l.target];
          if(target_node.expression != "None"){
            target_expression.push(parseFloat(target_node.expression));
          }
        })
        if(source_expression.length>0 && target_expression.length>0){
          let source_max = Math.max(...source_expression);
          let target_min = Math.min(...target_expression);
          if(Math.abs(source_max - target_min)>=threshold){
            this.motif.push(react);
          }
        }
      }
    }

    drawMotifSearchResult(motif_list){
      motif_list.sort(function(a,b){
        return d3.descending(a.pathways.length, b.pathways.length)
      })
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
          let mlinks = [];
          let mnodes_id = [];
          mnodes.push(motif_list[i]);
          let r_idx = 0;
          let p_idx = 0;
          motif_list[i].current_type = "reaction";
          motif_list[i].links.reactant.forEach(l=>{
            if(mnodes_id.indexOf(l.source)===-1){
              this.other_nodes[l.source].current_type = "reactant";
              this.other_nodes[l.source].r_idx = r_idx;
              mnodes.push(this.other_nodes[l.source]);
              mnodes_id.push(l.source);
              mlinks.push(l);
              r_idx += 1;
            }
          })
          motif_list[i].links.product.forEach(l=>{
            if(mnodes_id.indexOf(l.target)===-1){
              this.other_nodes[l.target].current_type = "product";
              this.other_nodes[l.target].p_idx = p_idx;
              mnodes.push(this.other_nodes[l.target]);
              mnodes_id.push(l.target);
              mlinks.push(l);
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
            .attr("r",4)
            .attr("stroke","gray")
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
            .text(d=>"#P: "+d.pathways.length)
            .style("font-size","11px")
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
      motif.links.reactant.forEach(l=>{
        if(mnodes_id.indexOf(l.source)===-1){
          this.other_nodes[l.source].current_type = "reactant";
          this.other_nodes[l.source].r_idx = r_idx;
          mnodes.push(this.other_nodes[l.source]);
          mnodes_id.push(l.source);
          mlinks.push(l);
          r_idx += 1;
        }
      })
      motif.links.product.forEach(l=>{
        if(mnodes_id.indexOf(l.target)===-1){
          this.other_nodes[l.target].current_type = "product";
          this.other_nodes[l.target].p_idx = p_idx;
          mnodes.push(this.other_nodes[l.target]);
          mnodes_id.push(l.target);
          mlinks.push(l);
          p_idx += 1;
        }
      })

      let motif_height = 80;

      // ***** draw motif glyph *****
      let x_interval = this.mp_svg_width/3;
      let y_interval_r = (motif_height-20)/r_idx;
      let y_interval_p = (motif_height-20)/p_idx;
      let ng = this.mp_motif_circle_group.selectAll("circle")
        .data(mnodes)
      ng.exit().remove();
      ng = ng.enter().append("circle").merge(ng)
        .attr("cx",(d)=>{
          if(d.current_type==="reactant"){
            return 50;
          } else if(d.current_type==="reaction"){
            return x_interval+50;
          } else if (d.current_type ==="product"){
            return x_interval*2+50;
          }
        })
        .attr("cy",(d)=>{
          if(d.current_type==="reactant"){
            return y_interval_r*(d.r_idx)+10;
          } else if(d.current_type==="reaction"){
            return motif_height/2;
          } else if (d.current_type ==="product"){
            return y_interval_p*(d.p_idx)+10;
          }
        })
        .attr("class",(d)=>d.current_type)
        .attr("r",8)
        .attr("stroke","gray")
        .attr("id",d=>"mp-circle-"+d.id);

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
        .attr("x",this.mp_svg_width/4)
        .attr("y",(motif_height-5))
        .text(d=>"Motif Id: "+d.id)
        .style("font-size","11px")
        .style("font-weight","bold")

    // **** draw pathway glyph ****
    let pathway_list = motif.pathways;
    let pathway_height = 20;
    let margin = {"horizontal":5, "vertical":10, "top":10, "left":5};
    this.mp_svg_height = Math.ceil(pathway_list.length/3) * (pathway_height+margin.vertical) + motif_height + margin.top;
    this.mp_svg.attr("height",this.mp_svg_height);

    let pathway_width = this.mp_svg_width/3-margin.horizontal;

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
        this.drawPathwayView(d);
        this.findAllMotif(d, this.motif);
        d3.select("#pathway-view-svg").style("visibility","visible");
        d3.select(".network-panel").style("visibility","visible");
      })
      .on("mouseover",(d)=>{
        d3.select("#mp-cover-"+d).style("opacity",0.4);
      })
      .on("mouseout",(d)=>{
        d3.select("#mp-cover-"+d).style("opacity",0);
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

      let ptg = this.mp_pathway_group.selectAll("text")
        .data(pathway_list);
      ptg.exit().remove();
      ptg = ptg.enter().append("text").merge(ptg)
        .attr("x",(d,i)=>margin.left + i%3*(pathway_width+margin.horizontal)+15)
        .attr("y",(d,i)=>motif_height+margin.top + Math.floor(i/3)*(pathway_height+margin.vertical)+12)
        .text(d=>d)
        .style("font-size",9)

    }

    drawPathwayView(p){
      // recover pathway
      let react_list = this.pathways[p].reactions;
      let pnodes = [];
      let plinks = [];
      let pnodes_id = [];

      react_list.forEach(rid=>{
        let react = this.react_nodes[rid];
        if(react){
          pnodes.push(react);
          react.links.reactant.forEach(l=>{
            if(pnodes_id.indexOf(l.source)===-1){
              pnodes.push(this.other_nodes[l.source]);
              pnodes_id.push(l.source);
            }
            plinks.push(l);
          })
          react.links.product.forEach(l=>{
            if(pnodes_id.indexOf(l.target)===-1){
              pnodes.push(this.other_nodes[l.target]);
              pnodes_id.push(l.target);
            }
            plinks.push(l);
          })
        }

      })
      console.log(pnodes, plinks)

      let simulation = d3.forceSimulation(pnodes)
        .force("link", d3.forceLink(plinks).id(function(d) { return d.id; }))
        .force("charge", d3.forceManyBody().strength(-1))
        .force("center", d3.forceCenter(this.pathway_svg_width / 2, this.pathway_svg_height / 2));

      let links = this.pathway_link_svg.selectAll("line").data(plinks);
      links.exit().remove();
      links = links.enter().append("line").merge(links);
      links
        .attr("stroke", (d)=>d.color)
        .attr("stroke-width",1)
        .attr("marker-end", "url(#end)");

      let nodes = this.pathway_circle_svg.selectAll("circle").data(pnodes);
      nodes.exit().remove();
      nodes = nodes.enter().append("circle").merge(nodes);
      nodes
        .attr("r",8)
        .attr("id",d=>"circle_"+this.createId(d.id))
        .attr("class",(d)=>{
          if(['reactant','product','reaction'].indexOf(d.type)!=-1){
            return d.type
          } else {
            return "other";
          }
        })
        .attr("stroke","gray")
        .attr("stroke-width",1)
        .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended))
        .on("mouseover",(d)=>{
          d3.select("#pathway-label-"+d.id)
            .style("visibility","visible");
        })
        .on("mouseout",(d)=>{
          d3.select("#pathway-label-"+d.id)
            .style("visibility","hidden");
        })
        .attr("id",(d)=>"pathway-circle-"+d.id);
      //     .on("mouseover",(d)=>{
      //         if(d.type === "reaction"){
      //             d3.select("#label_"+this.createId(d.id)).text(d=>d.id)
      //         }
      //         d3.select("#circle_"+this.createId(d.id)).attr("r",15)
      //         d3.select("#circle_frame_"+this.createId(d.id)).attr("r",15)
      //     })
      //     .on("mouseout",(d)=>{
      //         if(d.type === "reaction"){
      //             d3.select("#label_"+this.createId(d.id)).text("")
      //         }
      //         d3.select("#circle_"+this.createId(d.id)).attr("r",10)
      //         d3.select("#circle_frame_"+this.createId(d.id)).attr("r",10)
      //     });


      let labels = this.pathway_circle_svg.selectAll("text").data(pnodes);
      labels.exit().remove();
      let newLabels = labels.enter().append("text");
      labels = newLabels.merge(labels);
      labels
        .text(d=> {
          return d.id;
        })
        .attr("id",d=>"pathway-label-"+d.id)
        .style("font-size",10)
        .style("font-weight","bold")
        .style("visibility","hidden")

      simulation
        .nodes(pnodes)
          .on("tick", ticked);

      simulation.force("link")
        .links(plinks);

      let that = this;
      function ticked() {
        links
          .attr("x1", d => d.source.x)
          .attr("y1", d => d.source.y)
          .attr("x2", d => d.target.x)
          .attr("y2", d => d.target.y);

        let radius = 8;
        nodes
          .attr("cx", function(d) {
            return (d.x = Math.max(radius, Math.min(that.pathway_svg_width - radius, d.x)));
          })
          .attr("cy", function(d) {
            return (d.y = Math.max(radius, Math.min(that.pathway_svg_height/8*7 - radius, d.y)));
          })

        labels
          .attr("x",d=>d.x+5)
          .attr("y",d=>d.y-10);
      }

      function dragstarted(d) {
        if (!d3.event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
      }

      function dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
      }

      function dragended(d) {
        if (!d3.event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
      }
    }

    findAllMotif(pathway, motif_list){
      d3.select(".network-panel").style("visibility","visible");
      console.log(pathway, motif_list)
      let motif_id = [];
      motif_list.forEach(m=>{
        motif_id.push(m.id);
      })
      let pathway_list = this.pathways[pathway].reactions;
      let current_motif = [];
      pathway_list.forEach(react_id=>{
        if(motif_id.indexOf(react_id)!=-1){
          current_motif.push(react_id);
        }
      })
      current_motif.forEach(react_id=>{
        d3.select("#pathway-circle-"+react_id)
          .attr("stroke-width", 5)
          .attr("stroke", "gold")
      })

      let tg = d3.select("#all-motif-list").select("ul").selectAll("li").data(current_motif);
      tg.exit().remove();
      tg = tg.enter().append("li").merge(tg)
        .html(d=>"Motif Id: "+d)

    }


    computeAvg(arr){
      let arr_sum = arr[0];
      for(let i=1; i<arr.length; i++){
        arr_sum += arr[i];
      }
      let arr_avg = Math.round(arr_sum/arr.length*1000)/1000;
      return arr_avg;
    }



    computeMax(node_array){
      let node_max;
      node_array.forEach(nd=>{
        let nd_value = parseFloat(nd.expression)
        if(!Number.isNaN(nd_value)){
          if(node_max===undefined){
            node_max = nd_value;
          } else if(nd_value > node_max){
            node_max = nd_value;
          }
        }
      })
      return node_max;
    }

  computeMin(node_array){
    // if the length of the node array is 1, do not return min value (because max_val = min_val in this situation)
    let node_min;
    let node_length = 0;
    node_array.forEach(nd=>{
      let nd_value = parseFloat(nd.expression)
      if(!Number.isNaN(nd_value)){
        if(node_min===undefined){
          node_min = nd_value;
        } else if(nd_value > node_min){
          node_min = nd_value;
        }
        node_length += 1;
      }
    })
    if(node_length > 1){
      return node_min;
    } else {
      return;
    }
  }

  createId(id){
    return id.replace(/[^a-zA-Z0-9]/g, "")
  }
}
