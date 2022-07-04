import {csv, text} from "https://cdn.skypack.dev/d3-fetch@3";
import * as d3f from "https://cdn.skypack.dev/d3-fetch@3";

var fp_scores = "scores"; //"test.scores.csv"
var fp_species = "species"; //"test.species.csv" // "Euk.tree.annotations2.csv"
var fp_newick = "newick"; //"test.newick.txt"

var width = 1200;
var height = 900;
var y_scale = 1;
var x_scale = 1;
var treeWidth = width;
var treeHeight = height;
var duration = 500;
var font_size = 14;

var relatedNodes = []; // global variable for relatedness feature.
var lca; // last common ancenstor
var foundSpecies = []; //for searching a taxonomical unit
var nodeArray; //a global array of all the nodes

var annotations = await d3f.csv(fp_species);
var euk = await d3f.text("/newick");
var euk = await d3f.text(fp_newick);
var scores = await d3f.csv(fp_scores);

<!-- Copyright 2011 Jason Davies https://github.com/jasondavies/newick.js -->
function parseNewick(a) {
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

function zoom() {
    vis.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
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



/* Annotate using breadth-first traversal */
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

/*
Gets the total length of each node from the root (totalLength is this total 
length, and length is the distance from the parent node).
*/
function getTotalLength(node) {
    //console.log(node);
    if (typeof node.branchset != "undefined") {
        node.branchset.forEach(function(d) {
            d.totalLength = d.length + node.totalLength;
            getTotalLength(d);
        });
    }
}


/*This changes the y values to be based off of the actual relatedness distance.
Also an offset seemed potentially helpful, so it's there.*/
function adjustLength(n, offset) {
    y_scale = document.getElementById("rangeHorizontalScale").value * 10;
    x_scale = document.getElementById("rangeVerticalScale").value / 10;
    
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

/*
Describes how to handle clicks. If an internal node (a node with children) is clicked on, 
it will collapse the branches following this node. If a leaf node is clicked, the node 
will be highlighted or unhighlighted.
*/
function click(d) {
    // collapse subtree
    if (d.children) {
        d.children = null;
    } else {
        d.children = d.branchset;
    }

    // highlight edges in common
    if (!d.branchset) {
        // clicked on a leaf
        var nodeIndex = relatedNodes.indexOf(d);
        if (nodeIndex == -1) {
            relatedNodes.push(d);
        } else {
            relatedNodes.splice(nodeIndex, 1);
            traverseAncestors(d, false); // remove highlight
        }
    }
    
    if (relatedNodes.length > 1) {
        //colorRelated();
        relatedNodes.forEach(function(leaf) {
            traverseAncestors(leaf, true);
        });
    } else {
        relatedNodes.forEach(function(leaf) {
            traverseAncestors(leaf, false);
        });
    }
    update(d);
}


function name(node) {
    if (node.children) {
        return "";
    } else {
        return node.species;
    }
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

// looks through ancestors and sets them as in the path, until get to the LCA
function traverseAncestors(leaf, set_inpath) {
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

/*Sets the lca value with the id of the least common ancestor node.*/
function colorRelated() {
    relatedNodes.forEach(function(d) {
        makeAncestors(d);
    });

    var prior = 0;
    var one_node_ancestors = relatedNodes[0].ancestors;
    var j = 0;
    while (j < 25) { // arbitrary number, deeper than all leaves
        relatedNodes.slice(1, relatedNodes.length).forEach(function(d) {
            if (d.ancestors[j] != one_node_ancestors[j]) {
                j = 26; // to break out of the while
            }
        });
        if (j < 25) {
            prior = one_node_ancestors[j];
        }
        j++;
    }
    lca = prior;
}


// define the zoomListener which calls the zoom function on the "zoom" event constrained within the scaleExtents
// item 0: zoom out constraint, item 1: zoom in constraint
var zoomListener = d3.behavior.zoom().scaleExtent([0.1, 5]).on("zoom", zoom);

//the svg for the main tree
var vis = d3.select("body")
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("id", "visDiv")
    .call(zoomListener)
    .append("g")
    .attr("transform", "translate(50, 50)");

var tree = d3.layout.cluster()
    .size([treeWidth, treeHeight]);

var root = parseNewick(euk);

root.parent = {
    x: 0,
    y: 0
};

root.totalLength = 0;

traverseAndAnnotate(root);
getTotalLength(root);

window.vis = vis;
window.tree = tree;
window.root = root;

//for the search bar
d3.select("#searchButton").on("click", function() {
    var name = document.getElementById("nameInput").value;

    foundSpecies = [];
    var tax;

    nodeArray.forEach(function(d) { //going through each of the nodes in tree
        tax = d.taxonomy;
        tax.forEach(function(f) { //going through each entry in taxonomy array
            if (f.includes(name)) { //if has that name in it's taxonomy
                foundSpecies.push(d.id); //push it onto the array
            }
        });
    });

    if (foundSpecies.length == 0) {
        alert("No matches found");
    }

    update(root);
});

d3.select("#rangeHorizontalScale").on("change", function() {
    update(root);
});
d3.select("#rangeVerticalScale").on("change", function() {
    update(root);
});


d3.select("#zoomReset").on("click", function() {
    zoomListener.translate([0, 0]).scale(1);
    foundSpecies = [];
    update(root);
});


// https://gist.github.com/mlocati/7210513
function getLinkColour(score) {
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
    return color;
}

// https://stackoverflow.com/a/17268489/12891825
function getLinkColour_old(value){
    //value from 0 to 1
    var hue=((1-value)*120).toString(10);
    return ["hsl(",hue,",100%,50%)"].join("");
}


function update(source) {    
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
        .attr("r", 1e-6)
        .attr("fill", "dodgerblue")
        .on("mouseover", function(d) {
            d3.select("#commonName").text(d.commonName);
            d3.select("#score").text(d.score);
            d3.select("#taxonomy").text(d.name);
            d3.select("#depth").text(d.depth);
        })
        .classed("node", true);


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
                return 2.5;
            }
        });

    nodeUpdate.select("text")
        .text(function(d) {
            return name(d)
        });


    // when nodes are collapsed
    var nodeExit = node.exit().transition()
        .duration(duration)
        .attr("transform", function(d) {
            return "translate(" + source.y + "," + source.x + ")";
        });

    nodeExit.select("circle")
        .attr("r", 1e-15);

    nodeExit.select("text")
        .style("fill-opacity", 1e-6);





    // update the links
    var link = vis.selectAll("path.link")
        .data(links, function(d) {
            return d.target.id;
        });

    link.enter().insert("path", "g")
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
        .attr("stroke-width", function(d) {
            if (d.target.in_path) {
                return "4px";
            } else {
                return "2px";
            }
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
        })

    link.forEach(function (d) { 
        d.forEach(function (e) {
            e.style.stroke = getLinkColour(e.__data__.source.score); 
        });
    });

};
update(root);