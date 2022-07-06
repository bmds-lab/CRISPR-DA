import {csv, text} from "https://cdn.skypack.dev/d3-fetch@3";
import * as d3f from "https://cdn.skypack.dev/d3-fetch@3";

var fp_scores = "scores"; //"test.scores.csv"
var fp_species = "species"; //"test.species.csv" // "Euk.tree.annotations2.csv"
var fp_newick = "newick"; //"test.newick.txt"
var fp_guides = "guides";

var y_scale = 1;
var x_scale = 1;
var duration = 500;
var font_size = 14;

var lca; // last common ancenstor

var annotations = null;
var scores = null;
var root = null;

var tree = null;
var svg = null;
var vis = null;

async function update_data() {
    var root_tax_id = document.getElementById('txtRootTaxId').value;
    
    var e = document.getElementById("txtRootTaxId");
    var root_tax_id = e.options[e.selectedIndex].text;
    var root_tax_id = e.options[e.selectedIndex].value;
    
    var guide_seq = document.getElementById('selGuide').value;

    var newick_str = await d3f.text(fp_newick + "/" + root_tax_id);
    annotations = await d3f.csv(fp_species + "/" + root_tax_id);
    scores = await d3f.csv(fp_scores + "/" + root_tax_id + "/" + guide_seq);
    
    emptyTreeAndVis();
    
    root = parseNewick(newick_str);
    
    traverseAndAnnotate(root);
    getTotalLength(root);
    redrawTree();
}

function emptyTreeAndVis() {
    tree = d3.layout.cluster().size([1920, 1080]);

    var g = svg.select("g")
    if (g.size()) {
        svg.select('g').remove();
    }
    vis = svg.append("g")
       .attr("transform", "translate(50, 50)");
}

async function populateGuidesSelect() {
    var guides = await d3f.csv(fp_guides);
    var selGuide = d3.select("#selGuide");
    guides.forEach(function(d) {
        var seq = d.guide;
        selGuide.append("option").attr("value", seq).text(seq)
    });
}

function parseNewick(a) {
    // Copyright 2011 Jason Davies https://github.com/jasondavies/newick.js
    a = a.toString();
    for (var e = [], r = {}, s = a.split(/\s*(;|\(|\)|,|:)\s*/), t = 0; t < s.length; t++) {
        var n = s[t];
        switch (n) {
            case "(":
                var c = {};
                r.branchset = [c], e.push(r), r = c;
                break;
            case ",":
                var c = {};
                e[e.length - 1].branchset.push(c), r = c;
                break;
            case ")":
                r = e.pop();
                break;
            case ":":
                break;
            default:
                var h = s[t - 1];
                ")" == h || "(" == h || "," == h ? r.name = n : ":" == h && (r.length = parseFloat(n))
        }
    }
    return r
}

function findScore(node) {
    var number = node.name;
    scores.forEach(function(d) {
        if (+d.tax_id == +number) {
            node.score = d.score;
        }
    });
};

function findSpecies(node) {
    var number = node.name;
    annotations.forEach(function(d) {
        if (+d.taxId == +number) {
            node.species = d.organismName;
            node.commonName = d.commonName;
            node.taxonomy = d.taxonomy;
            node.taxonomyStr = d.taxonomy;
        }
    });
};

function traverseAndAnnotate(root_node) {
    root_node.children = root_node.branchset;
    var queue = [];
    queue.push(root_node);
    while (queue.length > 0) {
        var node = queue.shift(); // this is O(n)
        findSpecies(node);
        findScore(node);
        if (typeof node.branchset != "undefined") {
            node.branchset.forEach(function (d) {
               node.children = node.branchset;
               queue.push(d); 
            });
        }
    }
};

function getTotalLength(node) {
    //console.log(node);
    if (typeof node.branchset != "undefined") {
        node.branchset.forEach(function(d) {
            d.totalLength = d.length + node.totalLength;
            getTotalLength(d);
        });
    }
}

function adjustLength(n, offset) {
    //console.log(n);
    y_scale = 200 / Math.log(12 - document.getElementById("rangeHorizontalScale").value);
    x_scale = 0.1 / Math.log(52 - document.getElementById("rangeVerticalScale").value);
    
    if (n.length != null) {
        offset += n.length * y_scale;
    }
    
    n.y = offset;
    n.x = n.x * (font_size * x_scale);
    
    if (n.children) {
        n.children.forEach(function(n) {
            adjustLength(n, offset);
        });
    }
}

function rightAngleDiagonal(d, i) {
    return "M" + d.source.y + "," + d.source.x +
        "V" + d.target.x + "H" + d.target.y;
}

function bolded(is_mousedover) {
    return function(d) {
        d3.select(this).classed("link_bold", is_mousedover);
        d3.select(this.source).classed("link_bold", is_mousedover);
    }
}

function click(d) {
    // collapse subtree
    if (d.children) {
        d.children = null;
    } else {
        d.children = d.branchset;
    }
    redrawTree(d);
}

function name(node) {
    return node.commonName;
    //if (node.children) {
    //    return "";
    //} else {
    //    return node.commonName;
    //}
}

function makeAncestors(leaf) {
    var temp = leaf.parent;
    var ancestors = [leaf.id];
    while (temp.parent) {
        ancestors.unshift(temp.id);
        temp = temp.parent;
    }
    leaf.ancestors = ancestors;
}

function traverseAncestors(leaf, set_inpath) {
    // looks through ancestors and sets them as in the path, until get to the LCA
    // we have to clear the path out regardless
    var temp = leaf;
    while (leaf.id != 0) {
        leaf.in_path = false;
        leaf = leaf.parent;
    }
    
    leaf = temp;
    if (set_inpath) {
        while (leaf.id != lca) {
            //console.log(leaf);
            leaf.in_path = true;
            leaf = leaf.parent;
        }
    }
}

function collapseNthSubtree(root_node, n) {
    //console.log(root_node.commonName, n);
    if (root_node.children) {
        root_node.children.forEach(function (d) {
            collapseNthSubtree(d, n - 1);
        });
    }
    if (n <= 0) {
        if (root_node.children) {
            root_node.children = null;
        } else {
            root_node.children = root_node.branchset;
        }
    }    
    redrawTree(root_node);
}

function getLinkColour(score) {
    // https://gist.github.com/mlocati/7210513
    var perc = score;
    var r, g, b = 0;
    if(perc > 50) {
        r = 255;
        g = Math.round(5.1 * perc);
    }
    else {
        g = 255;
        r = Math.round(510 - 5.10 * perc);
    }
    var h = r * 0x10000 + g * 0x100 + b * 0x1;
    return '#' + ('000000' + h.toString(16)).slice(-6);
}

function getLinkColour_old(value){
    // https://stackoverflow.com/a/17268489/12891825
    //value from 0 to 1
    var hue=((1-value)*120).toString(10);
    return ["hsl(",hue,",100%,50%)"].join("");
}



function redrawTree(source) {
    var nodes = tree.nodes(root);
    var links = tree.links(nodes);

    adjustLength(nodes[0], 0);

    vis.selectAll(".link")
        .data(links, function(d) { 
            return d.target.id;
        })
        .enter()
        .append("path")
        .attr("class", "link")
        .attr("fill", "none")
        .attr("stroke", "black")
        .attr("d", rightAngleDiagonal);

    var node = vis.selectAll("g.node")
        .data(nodes, function(d, i) {
            return d.id || (d.id = i);
        });

    var nodeEnter = node.enter()
        .append("g")
        .attr("class", "node")
        .attr("transform", function(d) {
            return "translate(" + d.y + "," + d.x + ")";
        })
        .on("click", click);

    nodeEnter.append("circle")
        .attr("r", 8)
        .attr("fill", "black")
        .on("mouseover", function(d) {
            d3.select("#commonName").text(d.commonName);
            d3.select("#score").text(d.score);
            d3.select("#taxonomy").text(d.name);
            d3.select("#depth").text(d.depth);
        });

    nodeEnter.append("text")
        .attr("dx", 5)
        .attr("dy", 3)
        .style("font-size", font_size + "px")
        .text(function(d) {
            return name(d)
        })
        .on("mouseover", function(d) {
            d3.select("#commonName").text(d.commonName);
            d3.select("#score").text(d.score);
            d3.select("#taxonomy").text(d.name);
            d3.select("#depth").text(d.depth);
        });
     
    // when nodes are updated
    var nodeUpdate = node.transition()
        .duration(duration)
        .attr("transform", function(d) {
            return "translate(" + d.y + "," + d.x + ")";
        });

    nodeUpdate.select("circle")
        .attr("r", function(d) {
            if (d.id == lca) {
                return 5;
            } else {
                return 3;
            }
        });

    nodeUpdate.select("text")
        .style("fill-opacity", 1)
        .style("font-size", function (d) { 
            if (d.children) { 
                return "10px";
            } else {
                return "14px";
            } 
        })
        .text(function(d) {
            return name(d)
        });
        
    
    // when nodes are collapsed
    var nodeExit = node.exit().transition()
        .duration(duration)
        .attr("transform", function(d) {
            return "translate(" + source.y + "," + source.x + ")";
        });

    // make circle really small 
    nodeExit.select("circle")
      .attr("r", 1e-6);

    // make text invisible
    nodeExit.select("text")
      .style("fill-opacity", 1e-6);

    // update the links
    var link = vis.selectAll("path.link")
        .data(links, function(d) {
            return d.target.id;
        });

    link.enter().insert("path", "g")
      .attr("class", "link") // for some reason, this line is very important
      .attr("d", function(d) {
            var o = {
                x: d.source.x,
                y: d.source.y
            };
            return rightAngleDiagonal({
                source: o,
                target: o
            });
        })
        .attr("fill", "none");

    // Transition links to their new position.
    link.transition()
        .duration(duration)
        .attr("d", rightAngleDiagonal);

    // Transition exiting nodes to the parent's new position.
    link.exit().transition()
        .duration(duration)
        .attr("d", function(d) {
            var o = {
                x: source.x,
                y: source.y
            };
            return rightAngleDiagonal({
                source: o,
                target: o
            });
        });
        
    // Stash the old positions for transition.
    nodes.forEach(function(d) {
        d.x0 = d.x;
        d.y0 = d.y;
    });

    link.forEach(function (d) { 
        d.forEach(function (e) {
            e.style.stroke = getLinkColour(e.__data__.source.score); 
        });
    });
    
};

function zoom() {
    vis.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
}

var zoomListener = d3.behavior.zoom().scaleExtent([0.1, 3]).on("zoom", zoom);

d3.select("#btnUpdate").on("click", function() {
    update_data();
});

d3.select("#rangeHorizontalScale").on("change", function() {
    redrawTree(root);
});
d3.select("#rangeVerticalScale").on("change", function() {
    redrawTree(root);
});

d3.select("#zoomReset").on("click", function() {
    zoomListener.translate([0, 0]).scale(1);
    redrawTree(root);
});

svg = d3.select("#panelTree")
    .append("svg")
    .attr("width", '100%')
    .attr("height", '100%')
    .call(zoomListener);
    
await populateGuidesSelect();

await update_data();

window.vis = vis;
window.svg = svg;



redrawTree(root);
//collapseNthSubtree(root, 2);
