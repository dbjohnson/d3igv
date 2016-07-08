var bases = ["C", "G", "A", "T"],
    colorScale = d3.scale.category10().domain(bases),
    positionalSort = true,
    lastSortPos = 0,
    svgRefSeq = null,
    svgCoverage = null,
    svgReads = null;

///////////////////////////////////////////////////////////////

function initSVGs(width, margin, refseq, coverage, reads) {
  function makeSVG(track) {
    return d3.select(track.div).append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", track.height + margin.top)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
  }
  svgRefSeq = makeSVG(refseq);
  svgCoverage = makeSVG(coverage);
  svgReads = makeSVG(reads);
}

///////////////////////////////////////////////////////////////

function render(data, gridSizeX) {
  positionalSort = true
  sortReads(data, position=Math.floor(data.reference.length / 2));
  renderRefSeq(data, gridSizeX);
  renderCoverage(data, gridSizeX);
  renderReads(data, gridSizeX);
}

///////////////////////////////////////////////////////////////

function renderRefSeq(data, gridSizeX, trackHeight=10) {
  svgRefSeq.selectAll(".referenceLabel")
      .data(data.reference.split(""))
      .enter().append("text")
        .text(function (d) { return d; })
        .attr("x", function (d, i) { return i * gridSizeX; })
        .attr("y", 0)
        .style("text-anchor", "middle")
        .attr("transform", "translate(" + gridSizeX / 2 + ", "  + trackHeight / 1.5 + ")")
        .style("fill", function(d) { return colorScale(d); });
}

///////////////////////////////////////////////////////////////

function renderCoverage(data, gridSizeX, trackHeight=40) {
  var maxCoverage = d3.max(data.coverage, function(d) {
    return d3.sum(bases.map(function (b) { return d[b]; }));
  });

  var yscale = d3.scale.linear()
     .domain([0, maxCoverage])
    .range([0, trackHeight]);

  var stacked = []
      referenceBases = data.reference.split("");

  for (i in data.coverage) {
    var y = 0,
        cov = data.coverage[i],
        referenceBase = referenceBases[i],
        sortedBases = bases.slice().sort(function (a, b) {
          // sort bases by coverage so that the bases with the highest
          // read count go at the bottom of the stacked bar.  The reference
          // base always goes at the top (not interesting)
          if (b === referenceBase) {
            return -1;
          }
          else if (a === referenceBase) {
            return 1;
          }
          else {
            return cov[a] - cov[b]; }
        })


    for (b in sortedBases) {
      var base = sortedBases[b],
          height = yscale(cov[base])

      stacked.push({
        base: base,
        position: i,
        x: i * gridSizeX,
        y: trackHeight - y - height,
        height: height,
        color: base === referenceBase ? "grey" : colorScale(base),
        stackOrder: sortedBases
      })
      y += height;
    }
  }

  svgCoverage.selectAll("*").remove();
  var cards = svgCoverage.selectAll(".card").data(stacked);
  cards.enter().append("rect")
      .attr("x", function(d) { return d.x })
      .attr("y", function(d) { return d.y; })
      .attr("width", gridSizeX)
      .attr("height", function(d) { return d.height; })
      .attr('fill', function(d, i) { return d.color; } )
      .on("mouseover", function(d, i) {
          tooltip.html(d.stackOrder.map(function (b) {
                  return "<font color=" + colorScale(b) + ">" +
                         b + ": " + data.coverage[d.position][b];
           }).reverse().join("<br>"))
              .style("top", (d3.event.pageY - 28) + "px")
              .style("left", (d3.event.pageX + 10) + "px")
          tooltip.transition()
              .duration(200)
              .style("opacity", .9);
          })
      .on("mouseout", function(d) {
          tooltip.transition()
              .duration(500)
              .style("opacity", 0);
      })
      .on({"click": function(d) {
        if (d.position != lastSortPos) {
          positionalSort = true;
        }
        else {
          positionalSort = !positionalSort;
        }
        sortReads(data, position=positionalSort ? d.position : null);
        renderReads(data, gridSizeX);
      }});
}

///////////////////////////////////////////////////////////////

function renderReads(data, gridSizeX, trackHeight=500) {
  function expandReads(data) {
    var expanded = []
    for (var i = 0; i < data.reads.length; i++) {
      var read = data.reads[i]
          read_bases = read.sequence.split("");
      for (var pos = 0; pos < read_bases.length; pos++) {
        expanded.push({x: read.start + pos,
                       y: i,
                       fwd: read.fwd,
                       base: read_bases[pos]});
      }
    }
    return(expanded);
  }

  var oldCards = svgReads.selectAll("*");
  oldCards.remove()
  var newCards = svgReads.selectAll(".card").data(expandReads(data)),
      gridSizeY = trackHeight / data.reads.length,
      referenceBases = data.reference.split("");

  newCards.enter().append("rect")
      .attr("x", function(d) { return d.x * gridSizeX; })
      .attr("y", function(d) { return d.y * gridSizeY; })
      // .attr("rx", 2)
      // .attr("ry", 2)
      // .attr("class", "bordered")
      .attr("width", gridSizeX)
      .attr("height", gridSizeY)
      .attr("fill", function(d) { return d.base === referenceBases[d.x] ?
                                          (d.fwd ? "purple" : "pink") : colorScale(d.base); })
      ;
}

///////////////////////////////////////////////////////////////

function sortReads(data, position=null) {
  function compare(a, b, byPosition=true) {
    /*
    sort levels:
    1) read covers base
    2) read base coverage
    3) foward reads
    4) by start delay
    */

    if (byPosition && position) {
      lastSortPos = position
      function readAtBase(read) {
        return read.sequence[position - read.start];
      }

      function baseCoverage(base) {
        if (base === data.reference[position]) {
          // put wild-type reads at the bottom
          return 0;
        }
        else {
          return data.coverage[position][base];
        }
      }

      var readAbase = readAtBase(a),
          readBbase = readAtBase(b),
          delta = 0;

      if (readAbase && readBbase) {
        delta = baseCoverage(readBbase) - baseCoverage(readAbase);
        if (delta === 0) {
          if (readAbase === readBbase) {
            return compare(a, b, false)
          }
          delta = readAbase > readBbase ? 1 : -1;
        }
        return delta
      }
      else if (readAbase || readBbase) {
        delta = readAbase ? -10 : 10;
        return delta * data.reference.length + compare(a, b, false);
      }
      else {
        return compare(a, b, false);
      }
    }
    else {
      var s = 0
      if (a.fwd !== b.fwd) {
        s += a.fwd ? -data.reference.length : data.reference.length;
      }
      function startDelay(read) {
        return read.fwd ? read.start : data.reference.length - (read.start + read.sequence.length);
      }

      s += startDelay(a) - startDelay(b);
      return s;
    }
  }

  data.reads.sort(compare);
}
