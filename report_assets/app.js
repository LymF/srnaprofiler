(function () {
  "use strict";

  var D = JSON.parse(document.getElementById("report-data").textContent);
  var S = D.strings;
  var META = D.meta;
  var rows = D.rows;
  var SVGNS = "http://www.w3.org/2000/svg";
  var PAGE = 50;
  var CLASS_KEY = {piRNA: "pi", siRNA: "si", miRNA: "mi", phasiRNA: "ph", phasiRNA24: "ph",
    hcsiRNA: "hc", ambiguous: "am", other: "ot"};
  var LEVEL = ["absent", "weak", "present", "unavailable"];
  var STATUS_ICON = {absent: "-", weak: "!", present: "+", unavailable: "?"};
  var DOT_NAMES = ["Size", "Starts with U", "A at position 10", "Ping-pong", "Phasing"];

  // ---- DOM helpers: every string goes in through text nodes or setAttribute ----
  function h(tag, cls, kids) {
    var e = document.createElement(tag);
    if (cls) { e.className = cls; }
    add(e, kids);
    return e;
  }
  function add(e, kids) {
    if (kids === undefined || kids === null) { return e; }
    if (!Array.isArray(kids)) { kids = [kids]; }
    for (var i = 0; i < kids.length; i++) {
      var k = kids[i];
      if (k === undefined || k === null || k === false) { continue; }
      if (typeof k === "string" || typeof k === "number") {
        e.appendChild(document.createTextNode(String(k)));
      } else {
        e.appendChild(k);
      }
    }
    return e;
  }
  function sv(tag, attrs, kids) {
    var e = document.createElementNS(SVGNS, tag);
    for (var k in attrs) {
      if (Object.prototype.hasOwnProperty.call(attrs, k)) { e.setAttribute(k, String(attrs[k])); }
    }
    add(e, kids);
    return e;
  }
  function clear(e) { while (e.firstChild) { e.removeChild(e.firstChild); } return e; }
  function fmt(n) {
    if (n === null || n === undefined) { return "n/a"; }
    return Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }
  function pct(x) { return x === null || x === undefined ? "n/a" : Math.round(x * 100) + "%"; }
  function ckey(c) { return CLASS_KEY[c] || "ot"; }
  function tag(c) { return h("span", "rpt-srna-tag k-" + ckey(c), c || "unclassified"); }
  function where(r) { return r.h + ":" + fmt(r.s) + "-" + fmt(r.e); }
  function help(text) { return h("p", "rpt-srna-help", [h("b", "", "How to read: "), text]); }
  function card(title, body) {
    return h("div", "rpt-card", [
      h("div", "rpt-card__header", h("div", "rpt-card__title", title)),
      h("div", "rpt-card__body", body)]);
  }
  function stat(label, value, sub, mod) {
    return h("div", "rpt-stat-chip" + (mod ? " rpt-stat-chip--" + mod : ""), [
      h("div", "rpt-stat-chip__label", label),
      h("div", "rpt-stat-chip__value", value),
      sub ? h("div", "rpt-stat-chip__sub", sub) : null]);
  }
  function heading(title, sub) {
    return h("div", "rpt-section-header", h("div", "", [
      h("div", "rpt-section-title", title), sub ? h("div", "rpt-section-sub", sub) : null]));
  }
  function emptyMsg(text) { return h("p", "rpt-srna-empty", text); }

  // ---- tabs ----
  var TABS = [["overview", "Overview"], ["clusters", "Clusters"], ["glossary", "Glossary"],
    ["references", "References"]];
  var tabBtns = {};
  function show(id) {
    var secs = ["overview", "clusters", "detail", "glossary", "references"];
    for (var i = 0; i < secs.length; i++) {
      var s = document.getElementById("sec-" + secs[i]);
      s.className = "rpt-section" + (secs[i] === id ? " active" : "");
    }
    var hi = id === "detail" ? "clusters" : id;
    for (var k in tabBtns) {
      if (Object.prototype.hasOwnProperty.call(tabBtns, k)) {
        tabBtns[k].className = "rpt-nav__tab" + (k === hi ? " active" : "");
      }
    }
    window.scrollTo(0, 0);
  }
  function buildTabs() {
    var box = clear(document.getElementById("tabs"));
    TABS.forEach(function (t) {
      var b = h("button", "rpt-nav__tab", t[1]);
      b.setAttribute("type", "button");
      b.addEventListener("click", function () { show(t[0]); });
      tabBtns[t[0]] = b;
      box.appendChild(b);
    });
    add(clear(document.getElementById("navmeta")), META.title || "srnaprofiler");
  }

  // ---- charts ----
  function donut(counts) {
    var keys = Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; });
    var total = 0;
    keys.forEach(function (k) { total += counts[k]; });
    if (!total) { return emptyMsg("No clusters to draw."); }
    var svg = sv("svg", {viewBox: "0 0 330 170", width: "100%", role: "img"});
    var g = sv("g", {transform: "translate(85,85)"});
    svg.appendChild(g);
    var a = -Math.PI / 2, R = 70, r = 44;
    keys.forEach(function (k) {
      var f = counts[k] / total * 2 * Math.PI;
      if (f >= 2 * Math.PI - 1e-9) {
        g.appendChild(sv("circle", {r: (R + r) / 2, fill: "none", "stroke-width": R - r,
          "class": "f-" + ckey(k) + " s-" + ckey(k)}));
      } else {
        var big = f > Math.PI ? 1 : 0;
        var d = "M" + (R * Math.cos(a)) + "," + (R * Math.sin(a)) +
          "A" + R + "," + R + " 0 " + big + " 1 " + (R * Math.cos(a + f)) + "," + (R * Math.sin(a + f)) +
          "L" + (r * Math.cos(a + f)) + "," + (r * Math.sin(a + f)) +
          "A" + r + "," + r + " 0 " + big + " 0 " + (r * Math.cos(a)) + "," + (r * Math.sin(a)) + "Z";
        g.appendChild(sv("path", {d: d, "class": "f-" + ckey(k)}));
      }
      a += f;
    });
    g.appendChild(sv("text", {"text-anchor": "middle", y: 4, "class": "t-strong"}, fmt(total)));
    keys.slice(0, 6).forEach(function (k, i) {
      svg.appendChild(sv("rect", {x: 190, y: 30 + i * 22, width: 12, height: 12, rx: 3,
        "class": "f-" + ckey(k)}));
      svg.appendChild(sv("text", {x: 208, y: 41 + i * 22},
        k + " " + (counts[k] / total * 100).toFixed(1) + "%"));
    });
    return svg;
  }

  // Stacked bars. layers: [{name, cls, vals}], cats: labels, opts: {hi: index to mark}
  function bars(cats, layers, opts) {
    opts = opts || {};
    var W = 330, H = 170, L = 30, B = 22, T = 8, n = cats.length;
    var svg = sv("svg", {viewBox: "0 0 " + W + " " + H, width: "100%", role: "img"});
    if (!n) { return svg; }
    var mx = 0, i, j;
    for (i = 0; i < n; i++) {
      var s = 0;
      for (j = 0; j < layers.length; j++) { s += layers[j].vals[i] || 0; }
      if (s > mx) { mx = s; }
    }
    if (!(mx > 0)) { mx = 1; }
    var bw = (W - L - 6) / n, ph = H - B - T;
    svg.appendChild(sv("line", {x1: L, x2: W - 4, y1: H - B, y2: H - B, "class": "ax"}));
    svg.appendChild(sv("text", {x: L - 3, y: T + 8, "text-anchor": "end"}, opts.maxLabel ? opts.maxLabel(mx) : fmtShort(mx)));
    svg.appendChild(sv("text", {x: L - 3, y: H - B, "text-anchor": "end"}, "0"));
    var every = Math.max(1, Math.ceil(n / 12));
    for (i = 0; i < n; i++) {
      var y = H - B;
      for (j = 0; j < layers.length; j++) {
        var v = layers[j].vals[i] || 0;
        var hh = v / mx * ph;
        if (hh > 0) {
          y -= hh;
          var cl = "f-" + layers[j].cls + (opts.hi === i ? " hi" : "");
          svg.appendChild(sv("rect", {x: L + i * bw + 1, y: y, width: Math.max(1, bw - 2), height: hh,
            "class": cl}));
        }
      }
      if (i % every === 0 || opts.hi === i) {
        svg.appendChild(sv("text", {x: L + i * bw + bw / 2, y: H - 8, "text-anchor": "middle"}, cats[i]));
      }
    }
    if (opts.legend) { legend(svg, layers); }
    return svg;
  }
  function legend(svg, layers) {
    layers.forEach(function (ly, i) {
      svg.appendChild(sv("rect", {x: 262, y: 10 + i * 14, width: 9, height: 9, "class": "f-" + ly.cls}));
      svg.appendChild(sv("text", {x: 275, y: 18 + i * 14}, ly.name));
    });
  }
  function fmtShort(v) {
    if (v >= 1000) { return Math.round(v / 100) / 10 + "k"; }
    return v >= 10 ? String(Math.round(v)) : String(Math.round(v * 10) / 10);
  }

  function sizeHistogram() {
    var lo = 99, hi = 0, i;
    rows.forEach(function (r) {
      if (typeof r.z === "number" && r.z > 0) { lo = Math.min(lo, r.z); hi = Math.max(hi, r.z); }
    });
    if (lo > hi) { return emptyMsg("No dominant read size is available."); }
    lo = Math.max(16, lo); hi = Math.min(40, Math.max(hi, lo));
    var cats = [], byClass = {};
    for (i = lo; i <= hi; i++) { cats.push(String(i)); }
    rows.forEach(function (r) {
      if (typeof r.z !== "number" || r.z < lo || r.z > hi) { return; }
      if (!byClass[r.c]) { byClass[r.c] = cats.map(function () { return 0; }); }
      byClass[r.c][r.z - lo] += 1;
    });
    var layers = Object.keys(byClass).map(function (k) {
      return {name: k, cls: ckey(k), vals: byClass[k]};
    });
    return bars(cats, layers, {legend: true});
  }

  // ---- overview ----
  function renderOverview() {
    var box = clear(document.getElementById("sec-overview"));
    var lib = D.library || {};
    var amb = D.counts.ambiguous || 0;
    add(box, [
      heading("What is in this library?", "A plain-language summary before any technical detail."),
      h("div", "rpt-srna-say", META.summary)]);
    D.warnings.forEach(function (w) { box.appendChild(h("div", "rpt-srna-warn", w)); });
    add(box, [
      h("h3", "rpt-srna-h", "Main numbers"),
      h("div", "rpt-stats-grid", [
        stat("Clusters analysed", fmt(rows.length), "regions with many reads", "accent"),
        stat("Hidden by the depth filter", fmt(META.n_hidden),
          "of " + fmt(META.n_total) + " (below " + META.floor + " reads)", "warn"),
        stat("Ambiguous", fmt(amb), "not enough evidence for one class"),
        stat("Reads in analysed clusters", fmt(META.reads_shown), "", "success")]),
      h("h3", "rpt-srna-h", "Which kinds of small RNA?"),
      h("div", "rpt-srna-grid", [
        card("Classes found", [donut(D.counts),
          help("each slice is the share of clusters in one class. Ambiguous is not an error: " +
            "the evidence did not point to a single class.")]),
        card("Dominant read size (nt)", [sizeHistogram(),
          help("each cluster is counted once, at its most common read length, coloured by class. " +
            "About 21 to 22 nt suggests siRNA or miRNA; about 24 to 30 nt suggests piRNA.")])]),
      h("h3", "rpt-srna-h", "The library"),
      h("div", "rpt-stats-grid", [
        stat("Library size", fmt(lib.library_reads), "reads, estimated from RPM"),
        stat("Reads mapping to one place", pct(lib.unique_share), "weighted by cluster reads"),
        stat("Reads in all clusters", fmt(lib.cluster_reads), "before the depth filter"),
        stat("Fragments below the pipeline floor", fmt(lib.fragments),
          lib.fragment_reads === null || lib.fragment_reads === undefined ? "" :
            fmt(lib.fragment_reads) + " reads")])]);
  }

  // ---- cluster list ----
  var view = {cls: "all", sort: "conf", page: 0, floor: META.floor};
  function dots(v) {
    var box = h("span", "rpt-srna-dots");
    for (var i = 0; i < 5; i++) {
      var lv = LEVEL[parseInt(v.charAt(i), 10)] || "unavailable";
      var d = h("span", "dot " + lv);
      d.setAttribute("title", DOT_NAMES[i] + ": " + (S["status." + lv] || lv));
      box.appendChild(d);
    }
    return box;
  }
  function filtered() {
    var fl = Math.max(view.floor, META.floor);
    var out = [];
    rows.forEach(function (r, i) {
      if (view.cls !== "all" && r.c !== view.cls) { return; }
      if (r.r !== null && r.r < fl) { return; }
      out.push(i);
    });
    var cmp = {
      conf: function (a, b) { return (rows[b].f === null ? -1 : rows[b].f) - (rows[a].f === null ? -1 : rows[a].f) || rows[b].r - rows[a].r; },
      reads: function (a, b) { return rows[b].r - rows[a].r; },
      name: function (a, b) { return rows[a].n < rows[b].n ? -1 : rows[a].n > rows[b].n ? 1 : 0; }
    }[view.sort];
    out.sort(cmp);
    return out;
  }
  function renderList() {
    var box = clear(document.getElementById("sec-clusters"));
    add(box, heading("Which clusters deserve attention?",
      "Ordered by confidence by default. Click a row to see why."));

    var bar = h("div", "rpt-srna-filters");
    var classes = ["all"].concat(Object.keys(D.counts).sort());
    classes.forEach(function (c) {
      var b = h("button", "rpt-srna-chip" + (view.cls === c ? " on" : ""),
        c === "all" ? "All" : c + " (" + D.counts[c] + ")");
      b.setAttribute("type", "button");
      b.addEventListener("click", function () { view.cls = c; view.page = 0; renderList(); });
      bar.appendChild(b);
    });
    box.appendChild(bar);

    var ctl = h("div", "rpt-srna-filters");
    var sel = h("select", "rpt-srna-input");
    [["conf", "Sort: confidence"], ["reads", "Sort: reads"], ["name", "Sort: name"]].forEach(function (o) {
      var op = h("option", "", o[1]);
      op.setAttribute("value", o[0]);
      if (view.sort === o[0]) { op.setAttribute("selected", "selected"); }
      sel.appendChild(op);
    });
    sel.addEventListener("change", function () { view.sort = sel.value; view.page = 0; renderList(); });
    var inp = h("input", "rpt-srna-input");
    inp.setAttribute("type", "number");
    inp.setAttribute("min", "0");
    inp.setAttribute("value", String(view.floor));
    inp.addEventListener("change", function () {
      var v = parseFloat(inp.value);
      view.floor = isNaN(v) ? META.floor : v;
      view.page = 0;
      renderList();
    });
    add(ctl, [sel, h("label", "rpt-srna-lbl", ["Minimum reads ", inp])]);
    box.appendChild(ctl);
    box.appendChild(h("p", "rpt-srna-note",
      "This box can only narrow the table to the " + fmt(rows.length) + " clusters already in this file " +
      "(at least " + META.floor + " reads). To see clusters with fewer reads, rerun the report with " +
      "--min-reads set lower."));

    var idx = filtered();
    var pages = Math.max(1, Math.ceil(idx.length / PAGE));
    if (view.page >= pages) { view.page = pages - 1; }
    var slice = idx.slice(view.page * PAGE, view.page * PAGE + PAGE);

    var tb = h("tbody");
    if (!slice.length) {
      var td = h("td", "", emptyMsg(rows.length ? "No clusters match these filters." :
        "No clusters passed the depth filter, so there is nothing to list."));
      td.setAttribute("colspan", "6");
      tb.appendChild(h("tr", "", td));
    }
    slice.forEach(function (i) {
      var r = rows[i];
      var tr = h("tr", r.d ? "row" : "");
      var last;
      if (r.d) {
        tr.setAttribute("tabindex", "0");
        tr.addEventListener("click", function () { openCluster(i); });
        tr.addEventListener("keydown", function (ev) { if (ev.key === "Enter") { openCluster(i); } });
      }
      if (r.m && r.m.length) { last = "Missing: " + r.m.join("; "); }
      else { last = r.d ? "" : "no detail page (cap)"; }
      if (r.m && r.m.length && !r.d) { last += " (no detail page, cap)"; }
      var conf = r.f === null ? 0 : Math.max(0, Math.min(1, r.f));
      var fill = h("i");
      fill.setAttribute("style", "width:" + Math.round(conf * 100) + "%");
      add(tr, [
        h("td", "", [h("b", "", r.n), h("div", "rpt-srna-sub", where(r))]),
        h("td", "", tag(r.c)),
        h("td", "", [h("span", "rpt-srna-bar", fill), " " + (r.f === null ? "n/a" : r.f.toFixed(2))]),
        h("td", "", dots(r.v)),
        h("td", "num", fmt(r.r)),
        h("td", "rpt-srna-sub", last)]);
      tb.appendChild(tr);
    });
    var head = h("tr", "", ["Cluster", "Class", "Confidence", "Evidence", "Reads", ""].map(function (t) {
      return h("th", "", t);
    }));
    box.appendChild(h("div", "rpt-card rpt-srna-scroll",
      h("table", "rpt-srna-table", [h("thead", "", head), tb])));

    var pg = h("div", "rpt-srna-filters");
    var prev = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", "Previous");
    var next = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", "Next");
    prev.setAttribute("type", "button"); next.setAttribute("type", "button");
    if (view.page <= 0) { prev.setAttribute("disabled", "disabled"); }
    if (view.page >= pages - 1) { next.setAttribute("disabled", "disabled"); }
    prev.addEventListener("click", function () { view.page -= 1; renderList(); });
    next.addEventListener("click", function () { view.page += 1; renderList(); });
    add(pg, [prev, h("span", "rpt-srna-note", "Page " + (view.page + 1) + " of " + pages +
      " (" + fmt(idx.length) + " clusters)"), next]);
    box.appendChild(pg);

    box.appendChild(h("p", "rpt-srna-note", [
      "Evidence dots, in order: size, starts with U, A at position 10, ping-pong, phasing. ",
      h("span", "rpt-srna-dots", h("span", "dot present")), " present ",
      h("span", "rpt-srna-dots", h("span", "dot weak")), " weak ",
      h("span", "rpt-srna-dots", h("span", "dot absent")), " absent ",
      h("span", "rpt-srna-dots", h("span", "dot unavailable")), " not available in this table."]));
  }

  // ---- cluster page ----
  function lineChart(sig) {
    var W = 330, H = 170, L = 30, mid = 85, i;
    var svg = sv("svg", {viewBox: "0 0 " + W + " " + H, width: "100%", role: "img"});
    var plus = sig.cv[2], minus = sig.cv[3], n = Math.max(plus.length, minus.length);
    var mx = 0;
    for (i = 0; i < n; i++) { mx = Math.max(mx, plus[i] || 0, minus[i] || 0); }
    if (!(mx > 0)) { mx = 1; }
    var bw = n ? (W - L - 4) / n : 1;
    svg.appendChild(sv("line", {x1: L, x2: W - 4, y1: mid, y2: mid, "class": "ax"}));
    for (i = 0; i < n; i++) {
      var p = (plus[i] || 0) / mx * 66, m = (minus[i] || 0) / mx * 66;
      if (p > 0) { svg.appendChild(sv("rect", {x: L + i * bw, y: mid - p, width: Math.max(1, bw), height: p, "class": "f-plus"})); }
      if (m > 0) { svg.appendChild(sv("rect", {x: L + i * bw, y: mid, width: Math.max(1, bw), height: m, "class": "f-minus"})); }
    }
    svg.appendChild(sv("text", {x: L - 3, y: 18, "text-anchor": "end"}, fmtShort(mx)));
    svg.appendChild(sv("text", {x: L - 3, y: H - 40, "text-anchor": "end"}, fmtShort(mx)));
    svg.appendChild(sv("text", {x: L, y: H - 8}, fmt(sig.cv[0])));
    svg.appendChild(sv("text", {x: W - 4, y: H - 8, "text-anchor": "end"}, fmt(sig.cv[0] + n * sig.cv[1])));
    legend(svg, [{name: "+ strand", cls: "plus"}, {name: "- strand", cls: "minus"}]);
    return svg;
  }
  function figures(sig) {
    if (!sig) { return emptyMsg("Figures are not available for this cluster."); }
    var lo = sig.sz[0], cats = [], i;
    for (i = 0; i < sig.sz[1].length; i++) { cats.push(String(lo + i)); }
    var ntCats = [], A = [], C = [], G = [], U = [];
    sig.nt.forEach(function (p, k) {
      ntCats.push(String(k + 1));
      A.push(p[0] / 100); C.push(p[1] / 100); G.push(p[2] / 100); U.push(p[3] / 100);
    });
    var ppCats = [];
    for (i = 1; i <= sig.pp.length; i++) { ppCats.push(String(i)); }
    var ntSvg = bars(ntCats, [{name: "A", cls: "nA", vals: A}, {name: "C", cls: "nC", vals: C},
      {name: "G", cls: "nG", vals: G}, {name: "U", cls: "nU", vals: U}],
      {legend: true, maxLabel: function () { return "1"; }});
    return h("div", "rpt-srna-grid", [
      card("Read sizes by strand", [bars(cats, [{name: "+ strand", cls: "plus", vals: sig.sz[1]},
        {name: "- strand", cls: "minus", vals: sig.sz[2]}], {legend: true}),
        help("each bar is the number of reads of one length. Tight clusters have one tall bar.")]),
      card("Nucleotide at each position", [ntSvg,
        help("reads of the dominant size, base by base. A tall U at position 1 and A at position 10 " +
          "is the 1U/10A signature.")]),
      card("Coverage along the cluster", [lineChart(sig),
        help("reads per bin along the region. Above the line is the + strand, below is the - strand.")]),
      card("Overlap between strands (ping-pong)", [bars(ppCats, [{name: "pairs", cls: "plus", vals: sig.pp}],
        {hi: 9}), help("how often reads on opposite strands overlap by 1 to " + sig.pp.length +
        " bases. A peak at 10 means ping-pong.")])]);
  }
  function openCluster(i) {
    var r = rows[i], d = D.details[String(i)];
    var box = clear(document.getElementById("sec-detail"));
    var back = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", "Back to the list");
    back.setAttribute("type", "button");
    back.addEventListener("click", function () { show("clusters"); });
    box.appendChild(back);
    if (!d) {
      add(box, [heading(r.n, where(r)), emptyMsg("no detail page (cap)")]);
      show("detail");
      return;
    }
    var conf = r.f === null ? 0 : Math.max(0, Math.min(1, r.f));
    var fill = h("i");
    fill.setAttribute("style", "width:" + Math.round(conf * 100) + "%");
    add(box, [
      heading(r.n, where(r) + " | " + fmt(r.r) + " reads" +
        (typeof r.z === "number" ? " | dominant size " + r.z + " nt" : "")),
      h("div", "rpt-card rpt-srna-verdict", [
        tag(r.c),
        h("span", "rpt-srna-big", r.f === null ? "n/a" : r.f.toFixed(2)),
        h("span", "rpt-srna-bar wide", fill),
        h("span", "rpt-srna-note", "confidence, from 0 to 1; it is not a probability of being right")]),
      h("div", "rpt-srna-say", d.verdict),
      h("h3", "rpt-srna-h", "Why this class?")]);
    var list = h("div", "rpt-card rpt-card__body");
    d.checks.forEach(function (c) {
      var lv = LEVEL[parseInt(c[1], 10)] || "unavailable";
      list.appendChild(h("div", "rpt-srna-ev", [
        h("div", "ico " + lv, STATUS_ICON[lv]),
        h("div", "", [h("b", "", (S["check." + c[0] + ".title"] || c[0]) + ": " + (S["status." + lv] || lv)),
          h("div", "rpt-srna-sub", c[2])])]));
    });
    box.appendChild(list);
    add(box, [h("h3", "rpt-srna-h", "The signals, in figures"), figures(d.sig),
      h("h3", "rpt-srna-h", "Why not another class?")]);
    var wt = h("tbody");
    d.why.forEach(function (w) {
      wt.appendChild(h("tr", "", [h("td", "", w[0] ? tag(w[0]) : ""), h("td", "", w[1])]));
    });
    box.appendChild(h("div", "rpt-card", h("table", "rpt-srna-table", wt)));
    box.appendChild(h("h3", "rpt-srna-h", "Cautions"));
    if (d.cautions.length) {
      var ul = h("ul", "rpt-srna-list");
      d.cautions.forEach(function (c) { ul.appendChild(h("li", "", c)); });
      box.appendChild(h("div", "rpt-card rpt-card__body", ul));
    } else {
      box.appendChild(emptyMsg("No specific cautions for this cluster."));
    }
    add(box, [h("h3", "rpt-srna-h", "Dominant sequence"),
      h("div", "rpt-card rpt-card__body", d.seq ? h("code", "rpt-srna-seq", d.seq) :
        "No dominant sequence is recorded for this cluster.")]);
    show("detail");
  }

  // ---- glossary and references ----
  var refNum = {};
  D.references.forEach(function (r, i) { refNum[r.key] = i + 1; });
  function renderGlossary() {
    var box = clear(document.getElementById("sec-glossary"));
    add(box, heading("Glossary", "For readers new to small RNAs. Numbers in brackets point to the References tab."));
    var grid = h("div", "rpt-srna-grid");
    D.glossary.forEach(function (g) {
      var cites = (g.cite || []).map(function (k) { return refNum[k]; }).filter(Boolean);
      grid.appendChild(h("div", "rpt-card rpt-card__body", [
        h("h3", "rpt-srna-term", g.term), g.text,
        cites.length ? h("span", "rpt-srna-cite", " [" + cites.join(", ") + "]") : null]));
    });
    box.appendChild(grid);
  }
  function italic(text) {
    var parts = text.split("*"), out = [];
    parts.forEach(function (p, i) { out.push(i % 2 ? h("em", "", p) : p); });
    return out;
  }
  function renderRefs() {
    var box = clear(document.getElementById("sec-references"));
    add(box, heading("References"));
    var ol = h("ol", "rpt-srna-list");
    D.references.forEach(function (r) { ol.appendChild(h("li", "", italic(r.text))); });
    box.appendChild(h("div", "rpt-card rpt-card__body", ol));
  }

  buildTabs();
  renderOverview();
  renderList();
  renderGlossary();
  renderRefs();
  show("overview");
})();
