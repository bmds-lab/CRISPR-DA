import {csv, text} from "https://cdn.skypack.dev/d3-fetch@3";
import * as d3f from "https://cdn.skypack.dev/d3-fetch@3";

var fp_scores = "scores"; //"test.scores.csv"
var fp_subtree = "subtree"; //"test.newick.txt"
var fp_guides = "guides";

//var fp_species = "species"; //"test.species.csv" // "Euk.tree.annotations2.csv"
//var fp_subspecies = "subspecies"; //"test.species.csv" // "Euk.tree.annotations2.csv"
//var fp_newick = "newick"; //"test.newick.txt"
//var fp_subnewick = "subnewick"; //"test.newick.txt"

var start_tax_id = 11118; //10239;

var y_scale = 1;
var x_scale = 1;
var duration = 500;
var font_size = 14;

//var annotations = null;
//var scores = null;

var root = null;
var tree = null;
var svg = null;
var vis = null;
var count_nodes = 0;

async function emptyAndUpdateDataAndVis(root_tax_id=start_tax_id) {
    var guide_seq = document.getElementById('selGuide').value;

    // clear the SVG
    tree = d3.layout.cluster().size([1920, 1080]);
    var g = svg.select("g")
    if (g.size()) {
        svg.select('g').remove();
    }
    
    vis = svg.append("g");

    var subtree = await getAnnotatedSubTree(root_tax_id);
    root = subtree;

    redrawTree();
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


function getTotalLength(node) {
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
    //return "M" + d.target.y + "," + d.target.x
    //            + "C" + (d.source.y + 50) + "," + d.target.x
    //            + " " + (d.source.y + 150) + "," + d.source.x // 50 and 150 are coordinates of inflexion, play with it to change links shape
    //            + " " + d.source.y + "," + d.source.x;
    //          
    //          
    //return "M" + d.source.y + "," + d.source.x
    //            + "C" + (d.source.y) + " " + d.source.x
    //            + "," + (d.source.y+10) + " " + (d.source.x+10)
    //            // 50 and 150 are coordinates of inflexion, play with it to change links shape
    //            + ", " + d.target.y + " " + d.target.x;

    return "M" + d.source.y + "," + d.source.x +
        "V" + d.target.x + "H" + d.target.y;
}

async function getAnnotatedSubTree(root_tax_id) {
    var subtree = null;
    var guide_seq = document.getElementById('selGuide').value;
    var score_name = document.getElementById('selScore').value;
    await d3f.json(fp_subtree + "/" + root_tax_id + "/" + guide_seq + "/" + score_name)
        .then(function(json) {

            // parse the subtree newick string
            var newick = json['newick'];
            subtree = parseNewick(newick);

            // add the annotations and scores
            var annots = json['annotations']['commonName'];
            var scores = json['scores']
            var queue = [];
            queue.push(subtree);
            while (queue.length > 0) {
                var node = queue.shift(); // this is O(n)

                node.species = annots[node.name]; //d.organismName;
                node.commonName = annots[node.name]; //d.commonName;
                node.taxonomy = annots[node.name]; //d.taxonomy;
                node.taxonomyStr = annots[node.name]; //d.taxonomy;
                node.score = scores[node.name]

                if (typeof node.branchset != "undefined") {
                    node.branchset.forEach(function (node_child) {
                       node_child.parent = node;
                       node.children = node.branchset;
                       queue.push(node_child);
                    });
                }
            }
            //getTotalLength(subtree);
        })
        .catch(function(error) {
            console.log('something went wrong: ');
            console.log(error);
        });

    return subtree;
}

function doClick(node_clicked) {
    // collapse subtree
    if (node_clicked.children) {
        // remove children so the subtree collapses.
        // children are retained in `.branchset`.
        node_clicked.children = null;
        redrawTree(node_clicked);

    } else {
        if (node_clicked.branchset) {
            // the subtree will expand if `.branchset` contains the children
            node_clicked.children = node_clicked.branchset;
            redrawTree(node_clicked);

        } else {
            // `.branchset` was empty. Request the subtree, if it exists.
            var subtree = getAnnotatedSubTree(node_clicked.name)
                .then(function(st) {
                    st.parent = node_clicked;
                    node_clicked.children = st.children;
                    node_clicked.branchset = st.children;
                    redrawTree(st);
                });
        }
    }
    //redrawTree(node_clicked);
}

function name(node) {
    return node.commonName;
}

function getLinkColour(score) {
    // https://gist.github.com/mlocati/7210513

    //return '#ccc';
    

    if (score == -1 || typeof score == "undefined") {
        return '#ccc';
    }
    //if (score == 0) {
    //    return '#0000FF';
    //}

    var perc = score;
    var r, g, b = 0;
    if(perc < 50) {
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

function redrawTree(source) {
    var nodes = tree.nodes(root);
    var links = tree.links(nodes);

    adjustLength(nodes[0], 0);

    vis.selectAll(".link")
        .data(links, function(d) {
            //console.log(d);
            return d.target.id;
        })
        .enter()
        .append("path")
        .attr("class", "link")
        .attr("fill", "none")
        .attr("stroke", "#ccc")
        .attr("d", rightAngleDiagonal);


    var node = vis.selectAll("g.node")
        .data(nodes, function(d, i) {
            //console.log(i, d.species, d.name);
            count_nodes++;
            return d.id || (d.id = (count_nodes));
        });

    var nodeEnter = node.enter()
        .append("g")
        .attr("class", "node")
        .attr("transform", function(d) {
            return "translate(" + d.y + "," + d.x + ")";
        })
        .on("click", doClick);

    nodeEnter.append("circle")
        .attr("r", 5)
        .style("fill", "#69b3a2")
        .attr("stroke", "black")
        .style("stroke-width", 2)
        .on("mouseover", function(d) {
            d3.select("#commonName").text(d.commonName);
            d3.select("#score").text(d.score);
            d3.select("#taxonomy").text(d.name);
            d3.select("#depth").text(d.depth);
        });

    nodeEnter.append("text")
        .attr("dx", 8)
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

    // when nodes are updated/moved/etc.
    var nodeUpdate = node.transition()
        .duration(duration)
        .attr("transform", function(d) {
            return "translate(" + d.y + "," + d.x + ")";
        })
        .attr("opacity", 1)
        .attr("display", "inline");

    nodeUpdate.select("text")
        .style("font-size", function (d) {
            if (d.children) {
                return (font_size-2) + "px";
            } else {
                return font_size + "px";
            }
        })
        .attr('transform', function(d) {
            if (d.children) {
                return 'rotate(-6, 0, 0)';
            } else {
                return 'rotate(0, 0, 0)';
            }
        })
        .attr("dy", function(d) {
            if (d.children) {
                return -5;
            } else {
                return font_size/3;
            }
        });

    // when nodes are collapsed
    var nodeExit = node.exit().transition()
        .duration(duration)
        .attr("transform", function(d) {
            return "translate(" + source.y + "," + source.x + ")";
        })
        .attr("opacity", 0)
        .attr("display", "none");

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
            e.style.stroke = getLinkColour(e.__data__.target.score);
        });
    });

};

function zoom() {
    vis.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
}

var zoomListener = d3.behavior.zoom().scaleExtent([0.5, 5]).on("zoom", zoom);

d3.select("#btnUpdate").on("click", function(t) {
    var e = document.getElementById("txtRootTaxId");
    var root_tax_id = e.options[e.selectedIndex].text;
    var root_tax_id = e.options[e.selectedIndex].value;

    emptyAndUpdateDataAndVis(root_tax_id);
});

d3.select("#rangeHorizontalScale").on("change", function() {
    redrawTree(root);
});
d3.select("#rangeVerticalScale").on("change", function() {
    redrawTree(root);
});

//d3.select("#zoomReset").on("click", function() {
//    zoomListener.translate([0, 0]).scale(1);
//    redrawTree(root);
//});

svg = d3.select("#panelTree")
    .append("svg")
    .attr("width", '100%')
    .attr("height", '100%')
    .call(zoomListener);

populateGuidesSelect()
    .then(emptyAndUpdateDataAndVis);
    //.then(redrawTree);
    //.then(function() {
    //    collapseNthSubtree(root, 2);
    //});

window.vis = vis;
window.svg = svg;




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












// archive
/*
function findSpecies(node) {
    var number = node.name;
    annotations.forEach(function(d) {
        if (+d.taxId == +number) { // + converts to numerical representation
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


function bolded(is_mousedover) {
    return function(d) {
        d3.select(this).classed("link_bold", is_mousedover);
        d3.select(this.source).classed("link_bold", is_mousedover);
    }
}








            //// get annotations
            //d3f.csv(fp_subspecies + "/" + d.name)
            //    .then(function (csv) {
            //        csv.forEach(function (row) {
            //            annotations.push(row);
            //        });
            //    })
            //    .catch(function(error) {
            //        console.log('something went wrong: ');
            //        console.log(error);
            //    });
            //
            //console.log(annotations);
            //
            //// get newick string for subtree
            //d3f.text(fp_subnewick + "/" + d.name)
            //    .then(function (newick) {
            //        var subtree = parseNewick(newick);
            //        traverseAndAnnotate(subtree);
            //        getTotalLength(subtree);
            //
            //        //d.branchset = subtree.children;
            //        //d.children = d.branchset;
            //
            //        console.log(d);
            //        console.log(d.parent);
            //        console.log(subtree);
            //
            //        //d.parent = subtree;
            //        d.children = subtree.children;
            //        redrawTree(d);
            //    })
            //    .catch(function(error) {
            //        console.log('something went wrong: ');
            //        console.log(error);
            //    });






    .attr("y", function(d) {
            if (d.children) {
                console.log('adjusting y of label', d.y, d.y-10);
                return d.y-10;
            } else {
                return d.y
            }
        });

function findScore(node) {
    var number = node.name;
    scores.forEach(function(d) {
        if (+d.tax_id == +number) { // + converts to numerical representation
            node.score = d.score;
        }
    });
};
*/