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
  // reformat coverage to [[{base: "A", count: 13}, ...], [{base: "C", count: 0}, ...], ... ]
  // to play nice with d3.layout.stack
  var coverage = bases.map(function (b) {
    return data.coverage.map(function (d) {
        return {base: b,
                count: d[b]};
      })
  });

  var stack = d3.layout.stack()
    .x(function(d, i) { return i * gridSizeX; })
    .y(function(d) { return d.count; })
  var stacked = stack(coverage);

  var maxY = d3.max(stacked, function(d) {
    return d3.max(d, function(d) {
      return d.y0 + d.y;
    });
  });

  var y = d3.scale.linear()
    .domain([0, maxY])
    .range([trackHeight, 0]);

  var x = d3.scale.linear()
    .domain([0, data.reference.length])
    .range([0, gridSizeX * data.reference.length]);

  // bind a <g> tag for each layer
  var layers = svgCoverage.selectAll('g.layer')
    .data(stacked, function(d) { return d[0].base; })
      .enter()
        .append('g')
          .attr('class', 'layer');

  // bind a <rect> to each value inside the layer
  var referenceBases = data.reference.split("");
  layers.selectAll('rect')
    .data(function(d) { return d; })
    .enter()
      .append('rect')
        .attr('x', function(d, i) { return x(i); })
        .attr('width', gridSizeX)
        .attr('y', function(d) {
          return y(d.y0 + d.y);
        }).attr('height', function(d) {
          return trackHeight - y(d.y)
        })
        .attr('fill', function(d, i) {return d.base === referenceBases[i] ? "grey" : colorScale(d.base)})
        .on("mouseover", function(d, i) {
            tooltip.html(bases.map(function (b) { return b + ": " + data.coverage[i][b]; }).join("<br>"))
                .style("left", (d3.event.pageX + 10) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
            tooltip.transition()
                .duration(200)
                .style("opacity", .9);
            })
        .on("mouseout", function(d) {
            tooltip.transition()
                .duration(500)
                .style("opacity", 0);
        })
        .on({"click": function(d, i) {
          if (i != lastSortPos) {
            positionalSort = true;
          }
          else {
            positionalSort = !positionalSort;
          }
          sortReads(data, position=positionalSort ? i : null);
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

  svgReads.selectAll("*").remove();
  var cards = svgReads.selectAll(".card").data(expandReads(data))
      gridSizeY = trackHeight / data.reads.length;

  cards.append("title");
  cards.enter().append("rect")
      .attr("x", function(d) { return d.x * gridSizeX; })
      .attr("y", function(d) { return d.y * gridSizeY; })
      .attr("rx", 2)
      .attr("ry", 2)
      .attr("class", "bordered")
      .attr("width", gridSizeX)
      .attr("height", gridSizeY);

  var referenceBases = data.reference.split("");
  cards.transition().duration(500)
      .style("fill", function(d) {
        if (d.base === referenceBases[d.x]) {
          return(d.fwd ? "purple" : "pink")
        }
        else {
          return colorScale(d.base);
        }
      });

  cards.select("title").text(function(d) { return d.value; });
  cards.exit().remove();
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

  data.reads.sort(compare)
}
