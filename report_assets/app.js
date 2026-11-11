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
  var LEVEL = ["absent", "weak", "present", "unavailable", "low", "moderate", "high"];
  var STATUS_ICON = {absent: "-", weak: "!", present: "+", unavailable: "?",
    low: "+", moderate: "!", high: "-"};
  // colour of a level: multimapping reads low as good, high as bad
  var LEVEL_CSS = {absent: "absent", weak: "weak", present: "present", unavailable: "unavailable",
    low: "present", moderate: "weak", high: "absent"};
  var KT = META.kt || {};

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
  // text lookup with {name} placeholders
  function fmt(text, vars) {
    return String(text).replace(/\{(\w+)\}/g, function (m, k) {
      return vars && vars[k] !== undefined ? String(vars[k]) : m;
    });
  }
  function T(key, vars) { return S[key] === undefined ? key : fmt(S[key], vars); }
  function num(n) {
    if (n === null || n === undefined) { return T("ui.na"); }
    return Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }
  function pct(x) { return x === null || x === undefined ? T("ui.na") : Math.round(x * 100) + "%"; }
  function ckey(c) { return CLASS_KEY[c] || "ot"; }
  function tag(c) { return h("span", "rpt-srna-tag k-" + ckey(c), c || T("ui.unclassified")); }
  // 1-based coordinates (the TSV start is 0-based)
  function where(r) {
    return T("ui.where", {chrom: r.h, start: r.s === null ? T("ui.na") : num(r.s + 1),
      end: num(r.e)});
  }
  // default names are "chrom:start-end" (0-based): do not repeat them
  function autoName(r) { return r.n === r.h + ":" + r.s + "-" + r.e; }
  function confText(r) {
    if (r.c === "ambiguous") { return r.nc ? T("ui.conf.amb", {cls: r.nc}) : T("ui.conf.amb.none"); }
    return r.f === null ? T("ui.na") : r.f.toFixed(2);
  }
  function help(text) { return h("p", "rpt-srna-help", [h("b", "", T("ui.howto")), text]); }
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
  var TABS = [["overview", "ui.tab.overview"], ["clusters", "ui.tab.clusters"],
    ["glossary", "ui.tab.glossary"], ["references", "ui.tab.references"]];
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
      var b = h("button", "rpt-nav__tab", T(t[1]));
      b.setAttribute("type", "button");
      b.addEventListener("click", function () { show(t[0]); });
      tabBtns[t[0]] = b;
      box.appendChild(b);
    });
    add(clear(document.getElementById("navmeta")), META.title || T("ui.brand"));
  }

  // ---- charts ----
  function donut(counts) {
    var keys = Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; });
    var total = 0;
    keys.forEach(function (k) { total += counts[k]; });
    if (!total) { return emptyMsg(T("ui.donut.empty")); }
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
    g.appendChild(sv("text", {"text-anchor": "middle", y: 4, "class": "t-strong"}, num(total)));
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
    if (lo > hi) { return emptyMsg(T("ui.sizes.empty")); }
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
      heading(T("ui.ov.title"), T("ui.ov.sub")),
      h("div", "rpt-srna-say", META.summary)]);
    D.warnings.forEach(function (w) { box.appendChild(h("div", "rpt-srna-warn", w)); });
    add(box, [
      h("h3", "rpt-srna-h", T("ui.ov.main")),
      h("div", "rpt-stats-grid", [
        stat(T("ui.stat.analysed"), num(rows.length), T("ui.stat.analysed.sub"), "accent"),
        stat(T("ui.stat.hidden"), num(META.n_hidden),
          T("ui.stat.hidden.sub", {total: num(META.n_total), floor: META.floor}), "warn"),
        stat(T("ui.stat.amb"), num(amb), T("ui.stat.amb.sub")),
        stat(T("ui.stat.reads"), num(META.reads_shown), "", "success")]),
      h("h3", "rpt-srna-h", T("ui.ov.kinds")),
      h("div", "rpt-srna-grid", [
        card(T("ui.card.classes"), [donut(D.counts), help(T("ui.help.classes"))]),
        card(T("ui.card.sizes"), [sizeHistogram(), help(KT.sizes)])]),
      h("h3", "rpt-srna-h", T("ui.ov.library")),
      h("div", "rpt-stats-grid", [
        stat(T("ui.lib.size"), num(lib.library_reads), T("ui.lib.size.sub")),
        stat(T("ui.lib.unique"), pct(lib.unique_share), T("ui.lib.unique.sub")),
        stat(T("ui.lib.reads"), num(lib.cluster_reads), T("ui.lib.reads.sub")),
        stat(T("ui.lib.frag"), num(lib.fragments),
          lib.fragment_reads === null || lib.fragment_reads === undefined ? "" :
            T("ui.lib.frag.sub", {n: num(lib.fragment_reads)}))])]);
  }

  // ---- cluster list ----
  var view = {cls: "all", sort: "conf", page: 0, floor: META.floor};
  function dots(v) {
    var box = h("span", "rpt-srna-dots");
    for (var i = 0; i < 5; i++) {
      var lv = LEVEL[parseInt(v.charAt(i), 10)] || "unavailable";
      var d = h("span", "dot " + lv);
      d.setAttribute("title", (META.dot_names || [])[i] + ": " + T("status." + lv));
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
      // ambiguous rows after resolved ones, then confidence, then reads
      conf: function (a, b) {
        return (rows[a].c === "ambiguous") - (rows[b].c === "ambiguous") ||
          (rows[b].f === null ? -1 : rows[b].f) - (rows[a].f === null ? -1 : rows[a].f) ||
          (rows[b].r || 0) - (rows[a].r || 0);
      },
      reads: function (a, b) { return (rows[b].r || 0) - (rows[a].r || 0); },
      name: function (a, b) { return rows[a].n < rows[b].n ? -1 : rows[a].n > rows[b].n ? 1 : 0; }
    }[view.sort];
    out.sort(cmp);
    return out;
  }
  function renderList() {
    var box = clear(document.getElementById("sec-clusters"));
    add(box, heading(T("ui.list.title"), T("ui.list.sub")));

    var bar = h("div", "rpt-srna-filters");
    var classes = ["all"].concat(Object.keys(D.counts).sort());
    classes.forEach(function (c) {
      var b = h("button", "rpt-srna-chip" + (view.cls === c ? " on" : ""),
        c === "all" ? T("ui.all") : c + " (" + D.counts[c] + ")");
      b.setAttribute("type", "button");
      b.addEventListener("click", function () { view.cls = c; view.page = 0; renderList(); });
      bar.appendChild(b);
    });
    box.appendChild(bar);

    var ctl = h("div", "rpt-srna-filters");
    var sel = h("select", "rpt-srna-input");
    [["conf", T("ui.sort.conf")], ["reads", T("ui.sort.reads")], ["name", T("ui.sort.name")]].forEach(function (o) {
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
    add(ctl, [sel, h("label", "rpt-srna-lbl", [T("ui.floor.label"), inp])]);
    box.appendChild(ctl);
    box.appendChild(h("p", "rpt-srna-note",
      T("ui.floor.note", {n: num(rows.length), floor: META.floor})));

    var idx = filtered();
    var pages = Math.max(1, Math.ceil(idx.length / PAGE));
    if (view.page >= pages) { view.page = pages - 1; }
    var slice = idx.slice(view.page * PAGE, view.page * PAGE + PAGE);
    box.appendChild(h("p", "rpt-srna-note", T("ui.shown", {n: num(idx.length), m: num(rows.length)})));

    var tb = h("tbody");
    if (!slice.length) {
      var td = h("td", "", emptyMsg(rows.length ? T("ui.list.none.filter") : T("ui.list.none.rows")));
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
      if (r.u) {
        last = T(r.d ? "ui.row.unresolved" : "ui.row.unresolved.nodetail", {reason: r.u});
      } else {
        last = r.d ? "" : T("ui.row.nodetail");
      }
      var conf = r.f === null ? 0 : Math.max(0, Math.min(1, r.f));
      var fill = h("i");
      fill.setAttribute("style", "width:" + Math.round(conf * 100) + "%");
      add(tr, [
        h("td", "", [h("b", "", r.n), autoName(r) ? null : h("div", "rpt-srna-sub", where(r))]),
        h("td", "", tag(r.c)),
        h("td", "", [h("span", "rpt-srna-bar", fill), " " + confText(r)]),
        h("td", "", dots(r.v)),
        h("td", "num", num(r.r)),
        h("td", "rpt-srna-sub", last)]);
      tb.appendChild(tr);
    });
    var head = h("tr", "", [T("ui.th.cluster"), T("ui.th.class"), T("ui.th.conf"), T("ui.th.evidence"), T("ui.th.reads"), ""].map(function (t) {
      return h("th", "", t);
    }));
    box.appendChild(h("div", "rpt-card rpt-srna-scroll",
      h("table", "rpt-srna-table", [h("thead", "", head), tb])));

    var pg = h("div", "rpt-srna-filters");
    var prev = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", T("ui.page.prev"));
    var next = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", T("ui.page.next"));
    prev.setAttribute("type", "button"); next.setAttribute("type", "button");
    if (view.page <= 0) { prev.setAttribute("disabled", "disabled"); }
    if (view.page >= pages - 1) { next.setAttribute("disabled", "disabled"); }
    prev.addEventListener("click", function () { view.page -= 1; renderList(); });
    next.addEventListener("click", function () { view.page += 1; renderList(); });
    add(pg, [prev, h("span", "rpt-srna-note", T("ui.page.info",
      {p: view.page + 1, pages: pages, n: num(idx.length)})), next]);
    box.appendChild(pg);

    box.appendChild(h("p", "rpt-srna-note", [
      META.dot_legend + " ",
      h("span", "rpt-srna-dots", h("span", "dot present")), " " + T("status.present") + " ",
      h("span", "rpt-srna-dots", h("span", "dot weak")), " " + T("status.weak") + " ",
      h("span", "rpt-srna-dots", h("span", "dot absent")), " " + T("status.absent") + " ",
      h("span", "rpt-srna-dots", h("span", "dot unavailable")), " " + T("status.unavailable")]));
  }

  // ---- cluster page ----
  // bars are summed per-base depth per bin: divide by the bin for mean depth per base
  function lineChart(sig, end) {
    var W = 330, H = 170, L = 30, mid = 85, i;
    var svg = sv("svg", {viewBox: "0 0 " + W + " " + H, width: "100%", role: "img"});
    var plus = sig.cv[2], minus = sig.cv[3], n = Math.max(plus.length, minus.length);
    var mx = 0;
    for (i = 0; i < n; i++) { mx = Math.max(mx, plus[i] || 0, minus[i] || 0); }
    if (!(mx > 0)) { mx = 1; }
    var perBase = mx / (sig.cv[1] || 1);
    var bw = n ? (W - L - 4) / n : 1;
    svg.appendChild(sv("line", {x1: L, x2: W - 4, y1: mid, y2: mid, "class": "ax"}));
    for (i = 0; i < n; i++) {
      var p = (plus[i] || 0) / mx * 66, m = (minus[i] || 0) / mx * 66;
      if (p > 0) { svg.appendChild(sv("rect", {x: L + i * bw, y: mid - p, width: Math.max(1, bw), height: p, "class": "f-plus"})); }
      if (m > 0) { svg.appendChild(sv("rect", {x: L + i * bw, y: mid, width: Math.max(1, bw), height: m, "class": "f-minus"})); }
    }
    svg.appendChild(sv("text", {x: L - 3, y: 18, "text-anchor": "end"}, fmtShort(perBase)));
    svg.appendChild(sv("text", {x: L - 3, y: H - 40, "text-anchor": "end"}, fmtShort(perBase)));
    svg.appendChild(sv("text", {x: L, y: H - 8}, num(sig.cv[0] + 1)));
    svg.appendChild(sv("text", {x: W - 4, y: H - 8, "text-anchor": "end"},
      num(typeof end === "number" ? end : sig.cv[0] + n * sig.cv[1])));
    legend(svg, [{name: T("ui.strand.plus"), cls: "plus"}, {name: T("ui.strand.minus"), cls: "minus"}]);
    return svg;
  }
  function figures(sig, r) {
    if (!sig) { return emptyMsg(T("ui.fig.na")); }
    if (sig.err) { return emptyMsg(T("ui.fig.err")); }
    if (sig.none) { return emptyMsg(T("ui.fig.noreads")); }
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
      card(T("ui.fig.sizes"), [bars(cats, [{name: T("ui.strand.plus"), cls: "plus", vals: sig.sz[1]},
        {name: T("ui.strand.minus"), cls: "minus", vals: sig.sz[2]}], {legend: true}),
        help(T("ui.fig.sizes.help"))]),
      card(T("ui.fig.nt"), [ntSvg, help(KT.nt)]),
      card(T("ui.fig.cov"), [lineChart(sig, r.e), help(T("ui.fig.cov.help"))]),
      card(KT.pp_title, [bars(ppCats, [{name: T("ui.pairs"), cls: "plus", vals: sig.pp}],
        META.kingdom === "plant" ? {} : {hi: 9}),
        help(fmt(KT.pp_help, {n: sig.pp.length}))])]);
  }
  function openCluster(i) {
    var r = rows[i], d = D.details[String(i)];
    var box = clear(document.getElementById("sec-detail"));
    var back = h("button", "rpt-btn rpt-btn--ghost rpt-btn--sm", T("ui.detail.back"));
    back.setAttribute("type", "button");
    back.addEventListener("click", function () { show("clusters"); });
    box.appendChild(back);
    if (!d) {
      add(box, [heading(r.n, where(r)), emptyMsg(T("ui.row.nodetail"))]);
      show("detail");
      return;
    }
    var conf = r.f === null ? 0 : Math.max(0, Math.min(1, r.f));
    var fill = h("i");
    fill.setAttribute("style", "width:" + Math.round(conf * 100) + "%");
    var confNote = r.c === "ambiguous" ? T("ui.detail.conf.amb") :
      r.c === "miRNA" ? T("ui.detail.conf.mirna") : META.mode === "ml" ?
      T("ui.detail.conf.ml") : T("ui.detail.conf", {min: META.conf_min, k: META.conf_k});
    add(box, [
      heading(r.n, T("ui.detail.sub", {where: where(r), reads: num(r.r)}) +
        (typeof r.z === "number" ? T("ui.detail.size", {z: r.z}) : "")),
      h("div", "rpt-card rpt-srna-verdict", [
        tag(r.c),
        h("span", "rpt-srna-big", confText(r)),
        r.c === "ambiguous" ? null : h("span", "rpt-srna-bar wide", fill),
        h("span", "rpt-srna-note", confNote)]),
      h("div", "rpt-srna-say", d.verdict),
      h("h3", "rpt-srna-h", T("ui.detail.why"))]);
    var list = h("div", "rpt-card rpt-card__body");
    if (d.note) { list.appendChild(h("p", "rpt-srna-note", d.note)); }
    d.checks.forEach(function (c) {
      var lv = LEVEL[parseInt(c[1], 10)] || "unavailable";
      // checks that do not define the class are grey context, not red/green
      var css = c[3] ? LEVEL_CSS[lv] : "context";
      var line = T("ui.check.line", {title: T("check." + c[0] + ".title"), status: T("status." + lv)});
      list.appendChild(h("div", "rpt-srna-ev", [
        h("div", "ico " + css, STATUS_ICON[lv]),
        h("div", "", [h("b", "", c[3] ? line : line + " (" + T("check.context") + ")"),
          h("div", "rpt-srna-sub", c[2])])]));
    });
    box.appendChild(list);
    add(box, [h("h3", "rpt-srna-h", T("ui.detail.figs")), figures(d.sig, r),
      h("h3", "rpt-srna-h", T("ui.detail.whynot"))]);
    var wt = h("tbody");
    d.why.forEach(function (w) {
      wt.appendChild(h("tr", "", [h("td", "", w[0] ? tag(w[0]) : ""), h("td", "", w[1])]));
    });
    box.appendChild(h("div", "rpt-card", h("table", "rpt-srna-table", wt)));
    box.appendChild(h("h3", "rpt-srna-h", T("ui.detail.cautions")));
    if (d.cautions.length) {
      var ul = h("ul", "rpt-srna-list");
      d.cautions.forEach(function (c) { ul.appendChild(h("li", "", c)); });
      box.appendChild(h("div", "rpt-card rpt-card__body", ul));
    } else {
      box.appendChild(emptyMsg(T("ui.detail.nocautions")));
    }
    add(box, [h("h3", "rpt-srna-h", T("ui.detail.seq")),
      h("div", "rpt-card rpt-card__body", d.seq ? h("code", "rpt-srna-seq", d.seq) :
        T("ui.detail.noseq"))]);
    show("detail");
  }

  // ---- glossary and references ----
  var refNum = {};
  D.references.forEach(function (r, i) { refNum[r.key] = i + 1; });
  function renderGlossary() {
    var box = clear(document.getElementById("sec-glossary"));
    add(box, heading(T("ui.gl.title"), T("ui.gl.sub")));
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
    add(box, heading(T("ui.ref.title")));
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
