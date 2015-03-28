// useful function to see if an array contains an object
Array.prototype.contains = function(obj) {
	var i = this.length;
	while (i--) {
		if (this[i] == obj) {
			return true;
		}
	}
	return false;
};

// set up the margins and such
var W = parseInt(d3.select('#tree').style('width'))
var margin = {top: 1, right: 1, bottom: 1, left: 1},
    width = W - margin.left - margin.right,
	height = 1.2*W - margin.top - margin.bottom;

var zoo = "2";
var image_offset;
var tree;
function set_zoo() {
    d3.json("./config/zoo"+zoo+"_offset.json", function(d){
	    image_offset = d;
	    draw_tree();
    });
}

function draw_tree() {
    d3.select("#pdf_link")
        .attr("href", "./images/gz"+zoo+"_tree_crop.pdf")
    d3.json("./config/gz"+zoo+"_tree.json", function(d){
	    tree = d;
	    updateData(tree);
    });
}

set_zoo();
d3.selectAll("#zoo_buttons > label").on("click", function() {
    val=d3.select(this).select("input").property("value");
    if (zoo!=val) {
        zoo = val;
        set_zoo();
    }
});

var ky = height/16.0;
var kx = width/18.0;
//set_zoo();

// function that takes in a galaxy id and makes the node tree
function updateData(answers){
    // clear the page
    d3.selectAll("svg").remove();

    function update_charge(new_val){
	    d3.select("#slider_charge_value").text(new_val);
	    d3.select("#slider_charge").property("value", new_val);
	    force.charge(function(n) {return -1 * new_val * 1700 * n.value});
	    force.stop();
	    force.start();
    }

    function update_strength(new_val){
	    d3.select("#slider_strength_value").text(new_val);
	    d3.select("#slider_strength").property("value", new_val);
	    force.linkStrength(new_val);
	    force.stop();
	    force.start();
    }

    function update_friction(new_val){
	    d3.select("#slider_friction_value").text(new_val);
	    d3.select("#slider_friction").property("value", new_val);
	    // use 1-new_val to make 0 frictionless instead of 1!
	    force.friction(1-new_val);
	    force.stop();
	    force.start();
    }

    // add the draw window
    var svg = d3.select("#tree").append("svg")
	    .attr("width", width + margin.left + margin.right)
	    .attr("height", height + margin.top + margin.bottom)
        .append("g")
	    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .attr("id", "canvus");

    // create the node tree object
    var force = d3.layout.force()
	    .size([width, height]);
    
    // update the sliders to default values
    update_charge(1);
    update_strength(1);
    update_friction(0.35);
    json_callback(answers);

    // now that the basics are set up read in the json file
    function json_callback(answers) { 
	    // make sure to minpulate data *before* the update loop
	    // add a list of source and target Links to each node
        root = answers;
	    root.nodes.forEach(function(node) {
	        node.sourceLinks = [];
	        node.targetLinks = [];
	    });
	    root.links.forEach(function(link, i) {
	        // give each link a unique id
	        link.link_id = i;
	        var source = link.source,
		        target = link.target;
	        if (typeof source === "number") source = link.source = root.nodes[link.source];
	        if (typeof target === "number") target = link.target = root.nodes[link.target];
	        source.sourceLinks.push(link);
	        target.targetLinks.push(link);
	    });
	    // Get the number of votes for each node
	    root.nodes.forEach(function(node) {
            if (!node.name) {
                if (node.question=="gather" | node.question=="cp") {
                    node.value = 0;
                } else {
                    node.value = .05;
                }
            } else {
                node.value = .5;
            };
	    });
	    
	    // Normalize votes by total number
	    Total_value=root.nodes[0].value
	    root.nodes.forEach(function(node, i) {
	        // set the radius such that 9 full sized nodes could fit
	        node.radius =  width * Math.sqrt(node.value) / 25;
	        node.node_id = i;
	    });
        
	    // get the x position for each node
	    computeNodeBreadths(root);
	    // find how deep the tree goes and set the linkDistance to match 
	    max_level = 10;
	    force.linkDistance(.8*width/(max_level + 1));
	    
	    // good starting points
	    root.nodes.forEach(function(d , i) {
	        d.y = d.fixed_y;
	        d.x = d.fixed_x;
	    });

	    // run the call-back function to update positions
	    update(root.nodes, root.links);
    };

    // make the links long nice by using diagonal
    // swap x and y so the curve goes the propper way
    /*
    var diagonal = d3.svg.diagonal()
	    .source(function(d) { return d.flip ? {"x":d.source.y, "y":d.source.x} : {"x":d.source.x, "y":d.source.y}; })
	    .target(function(d) { return d.flip ? {"x":d.target.y, "y":d.target.x} : {"x":d.target.x, "y":d.target.y}; })
	    .projection(function(d) { return d.flip ? [d.y, d.x] : [d.x, d.y]; });
    */
    var diagonal = d3.svg.diagonal()
	    .source(function(d) { return {"x":d.source.x, "y":d.source.y}; })
	    .target(function(d) { return {"x":d.target.x, "y":d.target.y}; })
	    .projection(function(d) { return [d.x, d.y]; });
    
    var first_draw = true;
    var first_size = true;
    // create the update function to draw the tree
    function update(nodes_in, links_in) {
        function wrap(text) {
	        text.each(function() {
		        var text = d3.select(this),
                    width = text.attr("text_width"),
		            words = text.text().split(/\s+/).reverse(),
		            word,
		            line = [],
		            lineNumber = 0,
                    max_height = 1,
		            lineHeight = 1.1, // ems
		            y = text.attr("y"),
		            x = text.attr("x"),
		            dy = parseFloat(text.attr("dy")),
		            tspan = text.text(null).append("tspan").attr("x", x).attr("y", y).attr("dy", dy + "em");
		        while (word = words.pop()) {
		            line.push(word);
		            tspan.text(line.join(" "));
		            if (tspan.node().getComputedTextLength() > width) {
			            line.pop();
			            tspan.text(line.join(" "));
			            line = [word];
			            tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
			            if (lineNumber+1 > max_height) {
			                max_height = lineNumber+1;
			            }
		            }
		        }
                text.attr("text_height", max_height);
                d3.select(text.node().parentNode).attr("text_height", max_height);
                //text.attr("dy", 0.5*max_height+"em")
	        });
	    }
        
        var link = d3.select("#canvus").selectAll(".link");
	    var gnode = d3.select("#canvus").selectAll(".gnode");
	    // Set data as node ids        
	    // add the nodes and links to the tree
	    force
	        .nodes(nodes_in)
	        .links(links_in)
	        .on("tick", tick);

	    // set the data for the links (with unique ids)
	    link = link.data(links_in, function(d) { return d.link_id; });
        
	    // add a path object to each link
	    var lenter = link.enter().insert("path", ".gnode")
            .attr("class", function(d) { return d.dash ? "link dash" : "link"; })
	        .attr("d", diagonal)
            .style("stroke-width",1.5);
        
	    // Exit any old links
	    link.exit().remove();

	    // set the data for the nodes (with unique ids)
	    gnode = gnode.data(nodes_in, function(d) { return d.node_id; });

	    // Exit any old nodes
	    gnode.exit().remove();
        
	    // add a group to the node to translate it
	    var genter = gnode.enter().append("g")
	        .attr("class", function(d) { return d.answer_id ? "gnode" : "gnode metadata-thumbnail"; })
	        .call(force.drag);            
        
	    // add a group to the node to scale it
        // pull on the answer nodes and make them circles
	    var gimage = genter.filter(function(d) { return d.value==0.5; })
            .append("g")
            .attr("class", "gimage")
            .attr("transform", function(d) { return "scale(" + d.radius/50 + ")"; });

	    // add a clipPath for a circle to corp the node image
	    gimage.append("defs")
	        .append("clipPath")
	        .attr("id", function(d) { return "myClip" + d.node_id; })
	        .append("circle")
	        .attr("cx", 0)
	        .attr("cy", 0)
	        .attr("r", 45);

	    // add a black circle in the background
        /*
	    gimage.append("circle")
	        .attr("color", "black")
	        .attr("cx", 0)
	        .attr("cy", 0)
	        .attr("r", 45);
        */

        gimage.append("text")
            .attr("text-anchor", "middle")
            .attr("text_width", 90)
            .attr("dy", "0.3em")
            .attr("y",-60)
            .attr("x",0)
            .text(function(d) { return d.name; })
            .attr("style","font-size:15.8px")

        gimage.selectAll("text")
            .call(wrap);
        
        gimage.insert("g","text")
            .attr("transform", "translate(-50, -95)")
            .append("rect")
            .attr("class", function(d) { return "tier"+d.t; })
            .attr("width", 100)
            .attr("height", 140)
            .attr("rx", 20);
        
	    // add the inital image to the node
	    gimage.insert("image","text")
	        .attr("xlink:href", "./images/workflow_orig.png")
            .attr("id", "node_image")
	        .attr("x", -50)
	        .attr("y", function(d) { return -image_offset[d.answer_id][1]*100-50; })
	        .attr("clip-path", function(d) { return "url(#myClip" + d.node_id + ")"; })
	        .attr("width", 100)
	        .attr("height", 4900);
        
        // draw rect for question nodes
        var gquestion = genter.filter(function(d) { return d.value==0.05; })
            .append("g")
            .attr("class", "gquestion");
        
        gquestion.append("text")
            .attr("x", 0)
            .attr("y", 1)
            .attr("text_width", function(d) { return d.width*kx; })
            .attr("dy", "0.3em")
            .attr("text-anchor", "middle")
            .text(function(d) { return d.question; });

        gquestion.selectAll("text")
            .call(wrap);

        gquestion.insert("g", "text")
            .attr("text_height", function(d) { return d3.select(this.parentNode).attr("text_height"); })
            .attr("transform", function(d) { return 'translate(' + [-0.5*d.width*kx, -10] + ')'; })
            .attr("class", "question_box")
            .append("rect")
            .attr("class", function(d) { return "tier"+d.t; })
            .attr("width", function(d) { return d.width*kx; })
            .attr("height", function(d) { return parseInt(d3.select(this.parentNode).attr("text_height"))*1.4+'em'; })
            .attr("rx", 5)
            .attr("ry", 5);
	    // start the nodes moving
	    force.start();

	    // call-back to set how the nodes will move
	    function tick(e) {
	        // make sure the force gets smaller as the simulation runs
	        var kx = 10 * e.alpha;	        
	        root.nodes.forEach(function(d, i) {
                d.y = d.fixed_y
                d.x = d.fixed_x
	        });
	        // Translate the node group to the new position
	        gnode.attr("transform", function(d) {
                //if (d.width) {
                //    kx = width/18.0;
                //    new_x = d.x - 0.5*d.width*kx;
                //    new_y = d.y - 10;
                //    return 'translate(' + [new_x, new_y] + ')';
                //} else {
		            return 'translate(' + [d.x, d.y] + ')';
                //}
	        });    
	        link.attr("d",diagonal);            
	    };
    };
    // Find the x positions for each node
    function computeNodeBreadths(root) {
        root.nodes.forEach(function(node) {
            node.fixed_y = node.group_y * ky;
            node.fixed_x = node.group_x * kx;
	    });
    };
};
