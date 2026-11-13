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
  // colour of a level: multimapping reads low as good, high as bad
  var LEVEL_CSS = {absent: "absent", weak: "weak", present: "present", unavailable: "unavailable",
    low: "present", moderate: "weak", high: "absent"};
  var KT = META.kt || {};

  // ---- figure constants: the lab's R colours ----
  var BASES = ["A", "C", "G", "U"];
  var BASE_FILL = {A: "#108041", C: "#3B54A2", G: "#F9A31B", U: "#EA2627"};
  // letter ink inside a segment, by the fill's luminance
  var BASE_INK = {A: "#ffffff", C: "#ffffff", G: "#1e293b", U: "#ffffff"};
  var STRANDS = ["+", "-"];
  var STRAND_FILL = {"+": "#996631", "-": "#0099CC"};
  var PP_RED = "#e41a1c", GREY60 = "#999999", PERIOD_RED = "#e41a1c";
  var INK = "#1e293b", INK2 = "#475569", GRID = "#e2e8f0", AXIS = "#94a3b8", SURF = "#ffffff";
  var ACCENT = "#2f86d6", OTHER_BAR = "#94a3b8";
  var FS = 13;
  var FONT = "IBM Plex Sans, system-ui, sans-serif";
  var PERIODS = {phasiRNA: 21, phasiRNA24: 24};
  var BAR_MAX = 24, GAP = 2, RAD = 4, LEFT = 64, RIGHT = 14;

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
      if (Object.prototype.hasOwnProperty.call(attrs, k) && attrs[k] !== null && attrs[k] !== undefined) {
        e.setAttribute(k, String(attrs[k]));
      }
    }
    add(e, kids);
    return e;
  }
  function clear(e) { while (e.firstChild) { e.removeChild(e.firstChild); } return e; }
  function button(label, cls) {
    var b = h("button", cls || "rpt-btn rpt-btn--ghost rpt-btn--sm", label);
    b.setAttribute("type", "button");
    return b;
  }
  // text lookup with {name} placeholders
  function fmt(text, vars) {
    return String(text).replace(/\{(\w+)\}/g, function (m, k) {
      return vars && vars[k] !== undefined ? String(vars[k]) : m;
    });
  }
  function T(key, vars) { return S[key] === undefined ? key : fmt(S[key], vars); }

  // ---- numbers ----
  function commas(s) { return s.replace(/\B(?=(\d{3})+(?!\d))/g, ","); }
  function num(n) {
    if (n === null || n === undefined || !isFinite(n)) { return T("ui.na"); }
    return (n < 0 ? "-" : "") + commas(String(Math.round(Math.abs(n))));
  }
  // exact value with up to 3 decimals (weighted reads)
  function exact(n) {
    if (n === null || n === undefined) { return T("ui.na"); }
    var p = String(Math.round(n * 1000) / 1000).split(".");
    return commas(p[0]) + (p[1] ? "." + p[1] : "");
  }
  // a dequantised value: as many digits as the quantisation carries
  function val(v) {
    if (v === null || v === undefined || !isFinite(v)) { return T("ui.na"); }
    var a = Math.abs(v);
    if (a === 0) { return "0"; }
    if (a >= 100) { return num(v); }
    if (a >= 10) { return v.toFixed(1); }
    if (a >= 1) { return v.toFixed(2); }
    return String(Number(v.toPrecision(2)));
  }
  function pct(x) {
    return x === null || x === undefined ? T("ui.na") : T("ui.unit.pct", {v: Math.round(x * 100)});
  }
  function pval(p) { return p < 0.001 ? p.toExponential(1) : String(Number(p.toPrecision(2))); }
  function f2(v) { return Math.round(v * 100) / 100; }

  // ---- shared bits ----
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
  function card(title, body) {
    return h("div", "rpt-card", [
      h("div", "rpt-card__header", h("div", "rpt-card__title", title)),
      h("div", "rpt-card__body", body)]);
  }
  // text cut by an ellipsis stays readable in full on hover
  function titled(e, v) {
    if (typeof v === "string" || typeof v === "number") { e.setAttribute("title", String(v)); }
    return e;
  }
  function stat(label, value, sub, mod) {
    return h("div", "rpt-stat-chip" + (mod ? " rpt-stat-chip--" + mod : ""), [
      h("div", "rpt-stat-chip__label", label),
      titled(h("div", "rpt-stat-chip__value", value), value),
      sub ? titled(h("div", "rpt-stat-chip__sub", sub), sub) : null]);
  }
  function heading(title, sub, cls) {
    return h("div", "rpt-section-header" + (cls ? " " + cls : ""), h("div", "", [
      h("h2", "rpt-section-title", title), sub ? h("div", "rpt-section-sub", sub) : null]));
  }
  function emptyMsg(text) { return h("p", "rpt-srna-empty", text); }
  function note(text) { return h("p", "rpt-srna-note", text); }
  // the only annotation verdicts that may reach a class name
  var AN_OK = ["supports", "conflicts", "neutral"];
  // per-class table of annotation verdicts over the listed loci (meta.an_sum); text nodes only
  function annotationSummary(sum) {
    var head = h("tr", "", ["ui.an.th.class", "ui.an.th.supports", "ui.an.th.conflicts",
      "ui.an.th.neutral"].map(function (k) { return h("th", "", T(k)); }));
    var body = h("tbody", "", Object.keys(sum).map(function (c) {
      return h("tr", "", [h("td", "", c)].concat(sum[c].map(function (n) { return h("td", "", num(n)); })));
    }));
    return [h("h3", "rpt-srna-h", T("ui.an.title")), note(T("ui.an.note")),
      h("div", "rpt-srna-scroll", h("table", "rpt-srna-table", [h("thead", "", head), body]))];
  }

  // ---- tooltip: one element, filled with text nodes ----
  var tipEl = document.getElementById("tip");
  function showTip(t, x, y) {
    clear(tipEl);
    tipEl.appendChild(h("div", "hd", t.head));
    t.rows.forEach(function (rw) {
      var key = null;
      if (rw.key) {
        key = h("span", "key");
        key.style.borderTopColor = rw.key;
      }
      tipEl.appendChild(h("div", "rw", [key, rw.value === undefined ? null : h("strong", "", rw.value),
        h("span", "", rw.label)]));
    });
    tipEl.removeAttribute("hidden");
    var w = tipEl.offsetWidth, ht = tipEl.offsetHeight;
    var left = x + 14, top = y + 14;
    if (left + w > window.innerWidth - 8) { left = x - 14 - w; }
    if (top + ht > window.innerHeight - 8) { top = y - 14 - ht; }
    tipEl.style.transform = "translate(" + Math.round(Math.max(8, left)) + "px," +
      Math.round(Math.max(8, top)) + "px)";
  }
  function hideTip() { tipEl.setAttribute("hidden", ""); }
  function tipText(t) {
    return [t.head].concat(t.rows.map(function (rw) {
      return (rw.value === undefined ? "" : rw.value + " ") + rw.label;
    })).join("; ");
  }

  // ---- tabs ----
  var TABS = [["about", "ui.about.title"], ["overview", "ui.tab.overview"], ["clusters", "ui.tab.clusters"]];
  var SECTIONS = ["about", "overview", "clusters", "detail"];
  var tabBtns = {};
  var indicator = h("span", "rpt-srna-ind still");
  indicator.setAttribute("aria-hidden", "true");
  var current = "overview";
  function placeIndicator() {
    var b = tabBtns[current === "detail" ? "clusters" : current];
    if (!b) { return; }
    indicator.style.transform = "translateX(" + b.offsetLeft + "px) scaleX(" + (b.offsetWidth / 100) + ")";
  }
  function show(id) {
    current = id;
    for (var i = 0; i < SECTIONS.length; i++) {
      var s = document.getElementById("sec-" + SECTIONS[i]);
      s.className = "rpt-section" + (SECTIONS[i] === id ? " active" : "");
    }
    var hi = id === "detail" ? "clusters" : id;
    for (var k in tabBtns) {
      if (Object.prototype.hasOwnProperty.call(tabBtns, k)) {
        tabBtns[k].className = "rpt-nav__tab" + (k === "about" ? " rpt-nav__tab--about" : "") +
          (k === hi ? " active" : "");
        tabBtns[k].setAttribute("aria-selected", k === hi ? "true" : "false");
        tabBtns[k].setAttribute("tabindex", k === hi ? "0" : "-1");
      }
    }
    placeIndicator();
    hideTip();
    window.scrollTo(0, 0);
  }
  function buildTabs() {
    var box = clear(document.getElementById("tabs"));
    box.setAttribute("role", "tablist");
    TABS.forEach(function (t, i) {
      var b = button(T(t[1]), "rpt-nav__tab");
      b.setAttribute("role", "tab");
      b.setAttribute("id", "tab-" + t[0]);
      b.setAttribute("aria-controls", "sec-" + t[0]);
      b.addEventListener("click", function () { show(t[0]); });
      // arrows, Home and End move focus and selection along the tabs
      b.addEventListener("keydown", function (ev) {
        var n = TABS.length, j = {ArrowRight: i + 1, ArrowLeft: i - 1 + n, Home: 0, End: n - 1}[ev.key];
        if (j === undefined) { return; }
        if (ev.preventDefault) { ev.preventDefault(); }
        var id = TABS[j % n][0];
        show(id);
        tabBtns[id].focus();
      });
      tabBtns[t[0]] = b;
      box.appendChild(b);
    });
    SECTIONS.forEach(function (id) {
      var sec = document.getElementById("sec-" + id);
      sec.setAttribute("role", "tabpanel");
      sec.setAttribute("aria-labelledby", "tab-" + (id === "detail" ? "clusters" : id));
    });
    box.appendChild(indicator);
    add(clear(document.getElementById("navmeta")), META.title || T("ui.brand"));
  }

  // ---- About ----
  function icon(d, size, box) {
    var s = sv("svg", {width: size, height: size, viewBox: box, fill: "currentColor",
      "aria-hidden": "true", focusable: "false"});
    s.appendChild(sv("path", {d: d}));
    return s;
  }
  var GH_PATH = "M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49" +
    "-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82" +
    ".72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15" +
    "-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82" +
    " 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48" +
    " 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z";
  var ARROW_PATH = "M5 3l5 5-5 5-1.4-1.4L7.2 8 3.6 4.4z";
  function renderAbout() {
    var box = clear(document.getElementById("sec-about"));
    var links = h("div", "rpt-about__links");
    var gh = META.github || "";
    if (/^https:\/\//.test(gh)) {
      var a = h("a", "rpt-about__card", [
        h("span", "rpt-about__card-icon rpt-about__card-icon--github", icon(GH_PATH, 26, "0 0 16 16")),
        h("span", "rpt-about__card-body", [
          h("div", "rpt-about__card-title", T("ui.about.gh.title")),
          h("div", "rpt-about__card-sub", T("ui.about.gh.sub"))]),
        h("span", "rpt-about__card-arrow", icon(ARROW_PATH, 16, "0 0 16 16"))]);
      a.setAttribute("href", gh);
      a.setAttribute("target", "_blank");
      a.setAttribute("rel", "noopener noreferrer");
      links.appendChild(a);
    }
    add(box, h("div", "rpt-about", [
      h("div", "rpt-about__hero", [
        h("h1", "rpt-about__title", T("ui.about.name")),
        h("p", "rpt-about__tagline", META.tagline || ""),
        META.version ? h("span", "rpt-about__version", T("ui.about.version", {v: META.version})) : null]),
      links,
      h("div", "rpt-about__block", [
        h("h3", "rpt-srna-h", T("ui.about.install")),
        h("pre", "rpt-about__pre", h("code", "", (META.install || []).join("\n")))])]));
  }

  // ---- overview charts, in real pixels like the cluster figures ----
  var ovLive = [];
  function sized(height, label, draw) {
    var plot = h("div", "rpt-srna-fig__plot");
    plot.style.height = height + "px";
    var lastW = 0;
    function render() {
      var cw = Math.floor(plot.clientWidth || 0);
      if (!cw) {
        if (lastW) { return; }
        cw = 330;
      }
      if (cw === lastW) { return; }
      lastW = cw;
      clear(plot);
      var svg = svgRoot(cw, height, label);
      svg.setAttribute("role", "img");
      draw(svg, cw, height);
      plot.appendChild(svg);
    }
    ovLive.push({plot: plot, render: render});
    return plot;
  }
  // sorted horizontal bars, class name written on the row, so colour is never the only cue
  function donut(counts) {
    var keys = Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; });
    var total = 0;
    keys.forEach(function (k) { total += counts[k]; });
    if (!total) { return emptyMsg(T("ui.donut.empty")); }
    return sized(keys.length * 26 + 8, T("ui.card.classes"), function (svg, cw) {
      donutDraw(svg, cw, keys, counts, total);
    });
  }
  function donutDraw(svg, cw, keys, counts, total) {
    var labelW = 96, valW = 112, x0 = labelW + 8, maxW = Math.max(40, cw - x0 - valW);
    var top = counts[keys[0]];
    keys.forEach(function (k, i) {
      var y = 4 + i * 26, w = Math.max(2, Math.round(counts[k] / top * maxW));
      svg.appendChild(sv("text", {x: labelW, y: y + 15, "text-anchor": "end"}, k));
      svg.appendChild(sv("rect", {x: x0, y: y + 2, width: w, height: 18, rx: 3, "class": "f-" + ckey(k)}));
      svg.appendChild(sv("text", {x: x0 + w + 8, y: y + 15},
        (counts[k] / total * 100).toFixed(1) + "% \u00b7 " + num(counts[k])));
    });
  }
  // stacked bars of the dominant size per class
  function sizeHistogram() {
    var lo = 99, hi = 0, i;
    rows.forEach(function (r) {
      if (typeof r.z === "number" && r.z > 0) { lo = Math.min(lo, r.z); hi = Math.max(hi, r.z); }
    });
    if (lo > hi) { return emptyMsg(T("ui.sizes.empty")); }
    lo = Math.max(16, lo); hi = Math.min(40, Math.max(hi, lo));
    var n = hi - lo + 1, byClass = {}, tot = [];
    for (i = 0; i < n; i++) { tot.push(0); }
    rows.forEach(function (r) {
      if (typeof r.z !== "number" || r.z < lo || r.z > hi) { return; }
      if (!byClass[r.c]) { byClass[r.c] = tot.map(function () { return 0; }); }
      byClass[r.c][r.z - lo] += 1;
      tot[r.z - lo] += 1;
    });
    return sized(170, T("ui.card.sizes"), function (svg, W, H) { histDraw(svg, W, H, lo, n, byClass, tot); });
  }
  function histDraw(svg, W, H, lo, n, byClass, tot) {
    var L = 44, B = 24, TOP = 8, mx = Math.max.apply(null, tot) || 1, i;
    var slot = (W - L - 100) / n, bw = Math.min(BAR_MAX, slot - GAP), ph = H - B - TOP;
    svg.appendChild(sv("line", {x1: L, x2: W - 100, y1: H - B, y2: H - B, "class": "ax"}));
    svg.appendChild(sv("text", {x: L - 6, y: TOP + 10, "text-anchor": "end", "class": "tk-y"}, num(mx)));
    var every = Math.max(1, Math.ceil(26 / slot)), keys = Object.keys(byClass).sort();
    for (i = 0; i < n; i++) {
      var y = H - B;
      keys.forEach(function (k) {
        var hh = byClass[k][i] / mx * ph;
        if (hh > 0) {
          y -= hh;
          svg.appendChild(sv("rect", {x: L + i * slot + (slot - bw) / 2, y: y, width: bw, height: hh,
            "class": "f-" + ckey(k)}));
        }
      });
      if (i % every === 0) {
        svg.appendChild(sv("text", {x: L + i * slot + slot / 2, y: H - 6, "text-anchor": "middle",
          "class": "tk-x"}, String(lo + i)));
      }
    }
    keys.slice(0, 6).forEach(function (k, j) {
      svg.appendChild(sv("rect", {x: W - 92, y: 10 + j * 20, width: 12, height: 12, rx: 3, "class": "f-" + ckey(k)}));
      svg.appendChild(sv("text", {x: W - 74, y: 21 + j * 20}, k));
    });
  }

  // ---- library quality: read-length flags, trimming numbers, whole-library length figure ----
  function pct1(x) { return T("ui.unit.pct", {v: (Math.round(x * 1000) / 10).toFixed(1)}); }
  function lqLabel(i, n) { return i === 0 ? "<16" : i === n - 1 ? ">40" : String(15 + i); }
  function lqFlags(lq) {
    var box = h("div", "rpt-srna-strip");
    lq.fl.forEach(function (f) {
      var v = f[0] === "dominant" ? T("ui.unit.nt", {v: f[1]}) : pct1(f[1]);
      box.appendChild(h("div", "", [h("div", "lb", T("lq." + f[0] + ".label")), h("div", "vl", v),
        h("div", "lv " + (f[2] ? "present" : "absent"), T("lq." + f[0] + "." + f[2]))]));
    });
    return box;
  }
  function lqStats(lq, total) {
    var cards = [stat(T("lq.mapped"), num(total), T("lq.mapped.sub"))], t = lq.tr;
    if (lq.u) { cards.push(stat(T("lq.unmapped"), num(lq.u), T("lq.unmapped.sub"))); }
    if (t && t.b !== undefined) { cards.push(stat(T("lq.tr.before"), num(t.b))); }
    if (t && t.a !== undefined) {
      cards.push(stat(T("lq.tr.after"), num(t.a),
        t.b ? T("lq.tr.after.sub", {pct: pct1(t.a / t.b)}) : ""));
    }
    if (t && t.ad !== undefined) {
      cards.push(stat(T("lq.tr.adapter"), num(t.ad), t.b ? T("lq.tr.adapter.sub", {pct: pct1(t.ad / t.b)}) : ""));
    }
    return h("div", "rpt-stats-grid", cards);
  }
  function lqFigure(lq, total) {
    var c = lq.c, n = c.length, top = 0, i;
    for (i = 1; i < n; i++) { if (c[i] > c[top]) { top = i; } }
    var inBand = 0;
    for (i = 5; i <= 10; i++) { inBand += c[i]; }
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, y0 = 34, y1 = H - 44, len = y1 - y0;
      var tk = ticks(maxOf(c), 3), sc = len / tk.top;
      var cols = slots(x0, x1, n), slot = (x1 - x0) / n, bw = Math.min(BAR_MAX, slot - GAP);
      // the 20 to 25 nt band sits behind the bars and is named on the figure
      svg.appendChild(sv("rect", {x: f2(cols[5][0]), y: y0, width: f2(cols[10][1] - cols[5][0]), height: len,
        fill: GRID, "class": "band"}));
      svg.appendChild(text((cols[5][0] + cols[10][1]) / 2, 14, T("lq.band.label", {pct: pct1(inBand / total)}),
        {anchor: "middle", fill: INK, weight: 500, cls: "dl"}));
      yGrid(svg, x0, x1, y1, len, -1, tk);
      c.forEach(function (v, j) {
        if (v * sc >= 0.5) {
          svg.appendChild(tagMark(barPath(cols[j][0] + (slot - bw) / 2, y1, bw, v * sc, -1, "mk-bar", ACCENT), j));
        }
      });
      zeroLine(svg, x0, x1, y1);
      svg.appendChild(text(cols[top][0] + slot / 2, y1 - c[top] * sc - 6,
        T("ui.unit.nt", {v: lqLabel(top, n)}), {anchor: "middle", fill: INK, weight: 500, cls: "dl"}));
      var k = every(slot, 3), lab = [];
      for (i = 0; i < n; i++) {
        if (i % k === 0 || i === n - 1) { lab.push([cols[i][0] + slot / 2, lqLabel(i, n)]); }
      }
      xTicks(svg, lab, y1 + 18);
      xTitle(svg, T("ui.axis.length"), W, H - 6);
      yTitle(svg, T("ui.axis.reads"), (y0 + y1) / 2);
      hits(svg, cols, y0, y1, function (j) {
        return {head: T("ui.unit.nt", {v: lqLabel(j, n)}),
          rows: [{key: ACCENT, value: exact(c[j]), label: T("ui.axis.reads")}, {label: pct1(c[j] / total)}]};
      });
    }
    return figure({letter: "", title: T("ui.fig.lq"), n: T("ui.fig.n", {n: num(total)}), tip: T("ui.fig.lq.tip"),
      legend: [{kind: "sw", color: GRID, label: T("lq.band")}], height: 300, draw: draw,
      file: T("ui.brand") + "_lengths", live: ovLive, wide: true,
      table: function () {
        return {head: [T("ui.axis.length"), T("ui.axis.reads"), T("ui.axis.share")], exact: true,
          rows: c.map(function (v, j) { return [lqLabel(j, n), exact(v), pct1(v / total)]; })};
      }});
  }
  function libraryQuality(lq) {
    var total = 0;
    lq.c.forEach(function (v) { total += v; });
    return [h("h3", "rpt-srna-h", T("ui.ov.lq")), note(T("ui.ov.lq.note")), lqFlags(lq), lqStats(lq, total),
      lqFigure(lq, total)];
  }

  function renderOverview() {
    var box = clear(document.getElementById("sec-overview"));
    var lib = D.library || {};
    var amb = D.counts.ambiguous || 0;
    add(box, [
      heading(T("ui.ov.title"), T("ui.ov.sub"), "rpt-srna-enter"),
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
        card(T("ui.card.classes"), [donut(D.counts), note(T("ui.help.classes"))]),
        card(T("ui.card.sizes"), [sizeHistogram(), note(KT.sizes)])]),
      h("h3", "rpt-srna-h", T("ui.ov.library")),
      h("div", "rpt-stats-grid", [
        stat(T("ui.lib.size"), num(lib.library_reads), T("ui.lib.size.sub")),
        stat(T("ui.lib.unique"), pct(lib.unique_share), T("ui.lib.unique.sub")),
        stat(T("ui.lib.reads"), num(lib.cluster_reads), T("ui.lib.reads.sub")),
        stat(T("ui.lib.frag"), num(lib.fragments),
          lib.fragment_reads === null || lib.fragment_reads === undefined ? "" :
            T("ui.lib.frag.sub", {n: num(lib.fragment_reads)}))])]);
    if (META.an_sum) { add(box, annotationSummary(META.an_sum)); }
    if (lib.qc) { add(box, libraryQuality(lib.qc)); }
  }

  // ---- cluster list ----
  var view = {cls: "all", an: "all", sort: "conf", page: 0, floor: META.floor};
  var listH2 = null;
  var listFocus = null;  // data-ctl key of the control to refocus after a re-render
  function dots(v) {
    var box = h("span", "rpt-srna-dots");
    for (var i = 0; i < 5; i++) {
      var lv = LEVEL[parseInt(v.charAt(i), 10)] || "unavailable";
      var d = h("span", "dot " + lv);
      var lab = (META.dot_names || [])[i] + ": " + T("status." + lv);
      d.setAttribute("title", lab);
      d.setAttribute("role", "img");
      d.setAttribute("aria-label", lab);
      box.appendChild(d);
    }
    return box;
  }
  function evScore(v) {
    var s = 0;
    for (var i = 0; i < 5; i++) {
      var k = parseInt(v.charAt(i), 10);
      s += (k === 2 || k === 5 || k === 6) ? 10 : (k === 1 || k === 4) ? 1 : 0;
    }
    return s;
  }
  function filtered() {
    var fl = Math.max(view.floor, META.floor);
    var out = [];
    rows.forEach(function (r, i) {
      if (view.cls !== "all" && r.c !== view.cls) { return; }
      if (r.r !== null && r.r < fl) { return; }
      if (META.an && view.an !== "all" && r.an !== view.an) { return; }
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
      // number of supporting evidence checks (present/moderate/high), weak ones break ties
      evidence: function (a, b) {
        return evScore(rows[b].v) - evScore(rows[a].v) || (rows[b].r || 0) - (rows[a].r || 0);
      },
      name: function (a, b) { return rows[a].n < rows[b].n ? -1 : rows[a].n > rows[b].n ? 1 : 0; }
    }[view.sort];
    out.sort(cmp);
    return out;
  }
  function renderList() {
    var box = clear(document.getElementById("sec-clusters"));
    var refocus = null;
    function keep(e, key) {
      e.setAttribute("data-ctl", key);
      if (listFocus === key) { refocus = e; }
      return e;
    }
    var lhead = heading(T("ui.list.title"), T("ui.list.sub"));
    listH2 = lhead.firstChild.firstChild;
    add(box, lhead);

    var bar = h("div", "rpt-srna-filters");
    var classes = ["all"].concat(Object.keys(D.counts).sort());
    classes.forEach(function (c) {
      var b = button(c === "all" ? T("ui.all") : c + " (" + D.counts[c] + ")",
        "rpt-srna-chip" + (view.cls === c ? " on" : ""));
      b.setAttribute("aria-pressed", view.cls === c ? "true" : "false");
      keep(b, "cls:" + c);
      b.addEventListener("click", function () { view.cls = c; view.page = 0; listFocus = "cls:" + c; renderList(); });
      bar.appendChild(b);
    });
    box.appendChild(bar);

    var ctl = h("div", "rpt-srna-filters");
    var sel = h("select", "rpt-srna-input");
    [["conf", T("ui.sort.conf")], ["evidence", T("ui.sort.evidence")], ["reads", T("ui.sort.reads")], ["name", T("ui.sort.name")]].forEach(function (o) {
      var op = h("option", "", o[1]);
      op.setAttribute("value", o[0]);
      if (view.sort === o[0]) { op.setAttribute("selected", "selected"); }
      sel.appendChild(op);
    });
    sel.setAttribute("aria-label", T("ui.sort.aria"));
    keep(sel, "sort");
    sel.addEventListener("change", function () { view.sort = sel.value; view.page = 0; listFocus = "sort"; renderList(); });
    var inp = h("input", "rpt-srna-input");
    inp.setAttribute("type", "number");
    inp.setAttribute("min", "0");
    inp.setAttribute("autocomplete", "off");
    inp.setAttribute("aria-label", T("ui.floor.aria"));
    keep(inp, "floor");
    inp.setAttribute("value", String(view.floor));
    inp.addEventListener("change", function () {
      var v = parseFloat(inp.value);
      view.floor = isNaN(v) ? META.floor : v;
      view.page = 0;
      listFocus = "floor";
      renderList();
    });
    var anSel = null;
    if (META.an) {
      anSel = h("select", "rpt-srna-input");
      [["all", T("ui.an.all")], ["supports", T("ui.an.th.supports")],
        ["conflicts", T("ui.an.th.conflicts")], ["neutral", T("ui.an.th.neutral")]].forEach(function (o) {
        var op = h("option", "", o[1]);
        op.setAttribute("value", o[0]);
        if (view.an === o[0]) { op.setAttribute("selected", "selected"); }
        anSel.appendChild(op);
      });
      anSel.setAttribute("aria-label", T("ui.an.filter"));
      keep(anSel, "an");
      anSel.addEventListener("change", function () {
        view.an = AN_OK.indexOf(anSel.value) >= 0 ? anSel.value : "all";
        view.page = 0; listFocus = "an"; renderList();
      });
    }
    add(ctl, [sel, anSel ? h("label", "rpt-srna-lbl", [T("ui.an.filter") + " ", anSel]) : null,
      h("label", "rpt-srna-lbl", [T("ui.floor.label"), inp])]);
    box.appendChild(ctl);
    box.appendChild(note(T("ui.floor.note", {n: num(rows.length), floor: META.floor})));

    var idx = filtered();
    var pages = Math.max(1, Math.ceil(idx.length / PAGE));
    if (view.page >= pages) { view.page = pages - 1; }
    var slice = idx.slice(view.page * PAGE, view.page * PAGE + PAGE);
    var shown = note(T("ui.shown", {n: num(idx.length), m: num(rows.length)}));
    shown.setAttribute("aria-live", "polite");
    box.appendChild(shown);

    var tb = h("tbody");
    if (!slice.length) {
      var td = h("td", "", emptyMsg(rows.length ? T("ui.list.none.filter") : T("ui.list.none.rows")));
      td.setAttribute("colspan", META.an ? "7" : "6");
      tb.appendChild(h("tr", "", td));
    }
    slice.forEach(function (i) {
      var r = rows[i];
      var tr = h("tr", r.d ? "row" : "");
      var last;
      if (r.d) {
        tr.setAttribute("tabindex", "0");
        tr.addEventListener("click", function () { openCluster(i); });
        tr.addEventListener("keydown", function (ev) {
          if (ev.key === "Enter" || ev.key === " ") {
            if (ev.preventDefault) { ev.preventDefault(); }
            openCluster(i);
          }
        });
      }
      if (r.u) {
        last = T(r.d ? "ui.row.unresolved" : "ui.row.unresolved.nodetail", {reason: r.u});
      } else {
        last = r.d ? "" : T("ui.row.nodetail");
      }
      var conf = r.f === null ? 0 : Math.max(0, Math.min(1, r.f));
      var fill = h("i");
      fill.style.width = Math.round(conf * 100) + "%";
      add(tr, [
        h("td", "", [h("strong", "rpt-srna-name", r.n), autoName(r) ? null : h("div", "rpt-srna-sub", where(r))]),
        h("td", "", tag(r.c)),
        h("td", "", [h("span", "rpt-srna-bar", fill), " " + confText(r)]),
        h("td", "", dots(r.v)),
        META.an ? h("td", "", AN_OK.indexOf(r.an) >= 0 ? h("span", "rpt-srna-an " + r.an, r.an) : null) : null,
        h("td", "num", num(r.r)),
        h("td", "rpt-srna-sub", last)]);
      tb.appendChild(tr);
    });
    var heads = [T("ui.th.cluster"), T("ui.th.class"), T("ui.th.conf"), T("ui.th.evidence")];
    if (META.an) { heads.push(T("ui.an.filter")); }
    heads.push(T("ui.th.reads"), "");
    var readsAt = heads.length - 2;
    var head = h("tr", "", heads.map(function (t, j) {
      return h("th", j === readsAt ? "num" : "", t);
    }));
    box.appendChild(h("div", "rpt-card rpt-srna-scroll",
      h("table", "rpt-srna-table", [h("thead", "", head), tb])));

    var pg = h("div", "rpt-srna-filters");
    var prev = button(T("ui.page.prev"));
    var next = button(T("ui.page.next"));
    if (view.page <= 0) { prev.setAttribute("disabled", "disabled"); }
    if (view.page >= pages - 1) { next.setAttribute("disabled", "disabled"); }
    keep(prev, "prev");
    keep(next, "next");
    prev.addEventListener("click", function () { view.page -= 1; listFocus = "prev"; renderList(); });
    next.addEventListener("click", function () { view.page += 1; listFocus = "next"; renderList(); });
    add(pg, [prev, h("span", "rpt-srna-note", T("ui.page.info",
      {p: view.page + 1, pages: pages, n: num(idx.length)})), next]);
    box.appendChild(pg);

    box.appendChild(h("p", "rpt-srna-note", [
      META.dot_legend + " ",
      h("span", "rpt-srna-dots", h("span", "dot present")), " " + T("status.present") + " ",
      h("span", "rpt-srna-dots", h("span", "dot weak")), " " + T("status.weak") + " ",
      h("span", "rpt-srna-dots", h("span", "dot absent")), " " + T("status.absent") + " ",
      h("span", "rpt-srna-dots", h("span", "dot unavailable")), " " + T("status.unavailable")]));
    // a disabled page button cannot take focus; the list heading keeps the place instead
    if (refocus && refocus.hasAttribute && refocus.hasAttribute("disabled")) { refocus = null; }
    if (refocus && refocus.focus) { refocus.focus(); }
    listFocus = null;
  }

  // ==== figure helpers: SVG in real CSS pixels ====

  // [max, ints] back to values; top 99 (or 999 for the shared sizes max)
  function deq(pair, n, top) {
    var m = pair ? pair[0] : 0, q = pair ? pair[1] : [], out = [];
    for (var i = 0; i < n; i++) { out.push(i < q.length ? q[i] * m / (top || 99) : 0); }
    return out;
  }
  function sum(a) { var s = 0; for (var i = 0; i < a.length; i++) { s += a[i]; } return s; }
  function maxOf(a) { var m = 0; for (var i = 0; i < a.length; i++) { if (a[i] > m) { m = a[i]; } } return m; }

  function niceStep(max, count) {
    var raw = max / Math.max(1, count);
    var p = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
    var m = raw / p;
    return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * p;
  }
  function ticks(max, count) {
    if (!(max > 0)) { return {top: 1, step: 1, list: [0, 1], compact: false}; }
    var st = niceStep(max, count), top = Math.ceil(max / st - 1e-9) * st, list = [];
    for (var v = 0; v <= top + st / 2; v += st) { list.push(Math.round(v / st) * st); }
    return {top: top, step: st, list: list, compact: top >= 10000};
  }
  function tickLabel(v, tk) {
    if (v === 0) { return "0"; }
    if (tk.compact) {
      return v >= 1e6 ? T("ui.unit.m", {v: Math.round(v / 1e5) / 10}) : T("ui.unit.k", {v: Math.round(v / 100) / 10});
    }
    if (tk.step >= 1) { return num(v); }
    return v.toFixed(Math.min(4, Math.ceil(-Math.log(tk.step) / Math.LN10 - 1e-9)));
  }

  function svgRoot(W, H, label) {
    var s = sv("svg", {width: W, height: H, viewBox: "0 0 " + W + " " + H, role: "group",
      "aria-label": label, "font-family": FONT, "font-size": FS, fill: INK2});
    s.appendChild(sv("rect", {x: 0, y: 0, width: W, height: H, fill: SURF, "class": "bg"}));
    return s;
  }
  function text(x, y, s, o) {
    o = o || {};
    return sv("text", {x: f2(x), y: f2(y), "text-anchor": o.anchor || "start", fill: o.fill || INK2,
      "font-weight": o.weight || null, "data-ink": o.fill ? "1" : null, "class": o.cls || null, transform: o.transform || null,
      "pointer-events": "none"}, s);
  }
  function hline(svg, x0, x1, y, color, cls) {
    svg.appendChild(sv("line", {x1: f2(x0), x2: f2(x1), y1: f2(y), y2: f2(y), stroke: color,
      "stroke-width": 1, "shape-rendering": "crispEdges", "class": cls || null}));
  }
  function yTitle(svg, label, yMid) {
    svg.appendChild(text(0, 0, label, {anchor: "middle", cls: "yl",
      transform: "translate(14," + f2(yMid) + ") rotate(-90)"}));
  }
  function xTitle(svg, label, W, y) {
    svg.appendChild(text(LEFT + (W - LEFT - RIGHT) / 2, y, label, {anchor: "middle", cls: "xl"}));
  }
  // grid and tick labels of one half: dir -1 grows up from base, +1 grows down
  function yGrid(svg, x0, x1, base, len, dir, tk) {
    tk.list.forEach(function (t) {
      if (t === 0) { return; }
      var y = base + dir * t / tk.top * len;
      hline(svg, x0, x1, y, GRID, "grid");
      svg.appendChild(text(x0 - 8, y + 4.5, tickLabel(t, tk), {anchor: "end", cls: "tk-y"}));
    });
  }
  function zeroLine(svg, x0, x1, y) {
    hline(svg, x0, x1, y, AXIS, "zero");
    svg.appendChild(text(x0 - 8, y + 4.5, "0", {anchor: "end", cls: "tk-y"}));
  }
  // x tick labels at chosen positions
  function xTicks(svg, list, y) {
    var W = parseFloat(svg.getAttribute("width"));
    list.forEach(function (t) {
      // a label at the plot edge is pulled inside the figure instead of being cut
      var half = String(t[1]).length * 4, x = t[0], anchor = "middle";
      if (x + half > W - 2) { x = W - 2; anchor = "end"; } else if (x - half < 2) { x = 2; anchor = "start"; }
      svg.appendChild(text(x, y, t[1], {anchor: anchor, cls: "tk-x"}));
    });
  }
  // a bar with its data end rounded and its baseline square (absolute M/L/Q only)
  function barPath(x, base, w, len, dir, cls, fill) {
    var r = Math.min(RAD, w / 2, len), far = base + dir * len, k = -dir;
    var d = "M" + f2(x) + "," + f2(base) + "L" + f2(x) + "," + f2(far + k * r) +
      "Q" + f2(x) + "," + f2(far) + " " + f2(x + r) + "," + f2(far) +
      "L" + f2(x + w - r) + "," + f2(far) + "Q" + f2(x + w) + "," + f2(far) + " " + f2(x + w) + "," + f2(far + k * r) +
      "L" + f2(x + w) + "," + f2(base) + "Z";
    return sv("path", {d: d, fill: fill, "class": cls});
  }
  function segRect(x, a, b, w, cls, fill) {
    return sv("rect", {x: f2(x), y: f2(Math.min(a, b)), width: f2(w), height: f2(Math.abs(b - a)),
      fill: fill, "class": cls});
  }
  // stacked column: parts [{v, fill, letter}] from the baseline outwards, 2px surface gaps
  function stack(svg, x, w, base, dir, parts, scale, key) {
    var cum = 0, last = -1, i;
    for (i = 0; i < parts.length; i++) { if (parts[i].v > 0) { last = i; } }
    for (i = 0; i <= last; i++) {
      var p = parts[i];
      if (!(p.v > 0)) { continue; }
      var len = p.v * scale, a = cum + (cum > 0 ? GAP : 0), b = cum + len;
      cum = b;
      if (b - a < 0.5) { continue; }
      var cls = "mk-seg " + (dir < 0 ? "up" : "dn");
      svg.appendChild(tagMark(i === last ? barPath(x, base + dir * a, w, b - a, dir, cls, p.fill) :
        segRect(x, base + dir * a, base + dir * b, w, cls, p.fill), key));
      // letter only where it fits
      if (p.letter && b - a >= 16 && w >= 12) {
        svg.appendChild(tagMark(text(x + w / 2, base + dir * (a + b) / 2 + 4.5, p.letter,
          {anchor: "middle", fill: p.ink, weight: 500, cls: "lt"}), key));
      }
    }
    return cum;
  }
  function emptyNote(svg, x0, x1, y, key) {
    var t = text((x0 + x1) / 2, y + 4.5, T(key), {anchor: "middle", cls: "nostrand"});
    // a surface halo keeps the note readable over a gridline
    t.setAttribute("stroke", SURF);
    t.setAttribute("stroke-width", "4");
    t.setAttribute("paint-order", "stroke");
    svg.appendChild(t);
  }

  // Per-render state of the figure being drawn: marks by column key, focus restorers.
  var drawCtx = {marks: {}, focus: {}, sets: 0};
  // a bar, segment or letter belongs to column `key`; the hovered column lifts, the rest dim
  function tagMark(e, key) {
    if (key === undefined || key === null) { return e; }
    key = String(key);
    e.setAttribute("data-i", key);
    (drawCtx.marks[key] = drawCtx.marks[key] || []).push(e);
    return e;
  }
  // Contiguous hit columns: each is a tooltip and focus target (pointer, focus, arrow keys).
  // keyOf(i) names the marks of column i (default: i).
  function hits(svg, cols, y0, y1, tipOf, keyOf) {
    var els = [], on = -1, cur = [], marks = drawCtx.marks, setNo = drawCtx.sets++;
    keyOf = keyOf || function (i) { return String(i); };
    function mark(i) {
      if (on >= 0 && els[on]) { els[on].classList.remove("on"); }
      cur.forEach(function (m) { m.classList.remove("cur"); });
      on = i;
      cur = i >= 0 ? marks[keyOf(i)] || [] : [];
      if (i >= 0) { els[i].classList.add("on"); }
      cur.forEach(function (m) { m.classList.add("cur"); });
      svg.classList.toggle("lift", cur.length > 0);
    }
    function focusAt(j) {
      j = Math.max(0, Math.min(els.length - 1, j));
      els.forEach(function (e, k) { e.setAttribute("tabindex", k === j ? "0" : "-1"); });
      els[j].focus();
    }
    cols.forEach(function (c, i) {
      var e = sv("rect", {x: f2(c[0]), y: f2(y0), width: f2(Math.max(0.5, c[1] - c[0])), height: f2(y1 - y0),
        fill: "transparent", "class": "rpt-srna-hit", tabindex: i === 0 ? 0 : -1, role: "img",
        "data-h": setNo + ":" + i});
      drawCtx.focus[setNo + ":" + i] = function () { focusAt(i); };
      e.addEventListener("pointermove", function (ev) { mark(i); showTip(tipOf(i), ev.clientX, ev.clientY); });
      e.addEventListener("pointerdown", function (ev) { mark(i); showTip(tipOf(i), ev.clientX, ev.clientY); });
      // a lifted finger keeps the tooltip; the next touch elsewhere replaces it
      e.addEventListener("pointerleave", function (ev) {
        if (ev && ev.pointerType === "touch") { return; }
        mark(-1); hideTip();
      });
      e.addEventListener("focus", function () {
        var t = tipOf(i), bb = e.getBoundingClientRect();
        e.setAttribute("aria-label", tipText(t));
        mark(i);
        showTip(t, bb.right, bb.top);
      });
      e.addEventListener("blur", function () { mark(-1); hideTip(); });
      e.addEventListener("keydown", function (ev) {
        var step = {ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1, PageDown: 10, PageUp: -10}[ev.key];
        if (ev.key === "Home") { step = -i; } else if (ev.key === "End") { step = els.length - 1 - i; }
        if (ev.key === "Escape") { mark(-1); hideTip(); return; }
        if (step === undefined) { return; }
        if (ev.preventDefault) { ev.preventDefault(); }
        focusAt(i + step);
      });
      els.push(e);
      svg.appendChild(e);
    });
    return els;
  }
  function slots(x0, x1, n) {
    var w = (x1 - x0) / n, out = [];
    for (var i = 0; i < n; i++) { out.push([x0 + i * w, x0 + (i + 1) * w]); }
    return out;
  }
  // label every k-th category so labels never collide (13px digits are ~8px wide)
  function every(slotW, chars) { return Math.max(1, Math.ceil((chars * 8 + 10) / slotW)); }

  // ---- the figure block: panel letter, title with n, legend, table view, SVG download ----
  var live = [];
  var observer = null;
  function figure(o) {
    var plot = h("div", "rpt-srna-fig__plot" + (o.width ? " rpt-srna-fig__plot--scroll" : ""));
    plot.style.height = o.height + "px";
    var tbl = h("div", "rpt-srna-fig__table");
    tbl.setAttribute("hidden", "");
    tbl.setAttribute("tabindex", "0");
    tbl.setAttribute("role", "region");
    var lab = (o.letter ? o.letter + " " : "") + o.title;
    tbl.setAttribute("aria-label", lab + " " + T("ui.fig.table"));
    var leg = legend(o.legend || []);
    var svg = null, lastW = 0, tableOn = false;
    function render(force) {
      var cw = Math.floor(plot.clientWidth || 0);
      if (!cw) {
        if (lastW) { return; }
        cw = 640;
      }
      if (cw === lastW && !force) { return; }
      lastW = cw;
      var ae = document.activeElement;
      var fk = svg && ae && ae.parentNode === svg && ae.getAttribute ? ae.getAttribute("data-h") : null;
      clear(plot);
      // a figure wider than its box (o.width) scrolls inside it
      svg = svgRoot(o.width ? Math.max(cw, o.width) : cw, o.height, lab);
      drawCtx = {marks: {}, focus: {}, sets: 0};
      var drawn = o.draw(svg, cw, o.height);
      // a draw that sizes itself (the genome tracks) returns the height it used
      if (typeof drawn === "number") { plot.style.height = drawn + "px"; }
      var restore = drawCtx.focus;
      plot.appendChild(svg);
      // a redraw under keyboard focus keeps the focused column
      if (fk && restore[fk]) { restore[fk](); }
    }
    var tools = h("div", "rpt-srna-fig__tools");
    if (o.table) {
      var tb = button(T("ui.fig.table"));
      tb.setAttribute("aria-pressed", "false");
      tb.addEventListener("click", function () {
        tableOn = !tableOn;
        if (tableOn && !tbl.firstChild) {
          var t = o.table();
          add(tbl, [t.plain ? null : note(t.exact ? T("ui.fig.exact") : T("ui.fig.precision")), dataTable(t.head, t.rows, t.left)].concat(t.more || []));
        }
        if (tableOn) { plot.setAttribute("hidden", ""); tbl.removeAttribute("hidden"); } else {
          tbl.setAttribute("hidden", ""); plot.removeAttribute("hidden"); render(true);
        }
        tb.setAttribute("aria-pressed", tableOn ? "true" : "false");
        clear(tb);
        add(tb, T(tableOn ? "ui.fig.chart" : "ui.fig.table"));
        hideTip();
      });
      var dl = button(T("ui.fig.svg"));
      dl.addEventListener("click", function () {
        if (!svg) { render(true); }
        download(svg, o.file + ".svg");
      });
      add(tools, [tb, dl]);
    }
    var title = h("h3", "rpt-srna-fig__title", [o.letter ? h("span", "rpt-srna-fig__letter", o.letter) : null,
      h("span", "", o.title), o.n ? h("span", "rpt-srna-fig__n", o.n) : null]);
    if (o.tip) { title.setAttribute("title", o.tip); }
    var el = h("section", "rpt-srna-fig" + (o.wide ? " rpt-srna-fig--wide" : ""), [
      h("div", "rpt-srna-fig__head", [title, tools]), leg,
      o.body || plot, o.body ? null : tbl]);
    el.setAttribute("data-fig", o.letter || o.id || o.file);
    if (!o.body) { (o.live || live).push({plot: plot, render: render}); }
    return el;
  }
  function legend(items) {
    if (!items.length) { return null; }
    var ul = h("ul", "rpt-srna-legend");
    items.forEach(function (it) {
      var key = null;
      if (it.kind === "sw") {
        key = h("span", "sw");
        key.style.background = it.color;
      } else if (it.kind === "ln" || it.kind === "dash") {
        key = h("span", it.kind === "dash" ? "ln-dash" : "ln");
        key.style.borderTopColor = it.color;
      }
      if (key) { key.setAttribute("aria-hidden", "true"); }
      ul.appendChild(h("li", it.head ? "hd" : "", [key, it.label]));
    });
    return ul;
  }
  // `left` columns (default 1) are text; the rest are numbers
  function dataTable(head, body, left) {
    var nl = left || 1;
    var tr = h("tr", "", head.map(function (t, j) { return h("th", j >= nl ? "num" : "", t); }));
    var tb = h("tbody");
    body.forEach(function (r) {
      tb.appendChild(h("tr", "", r.map(function (c, j) { return h("td", j >= nl ? "num" : "", c); })));
    });
    return h("table", "rpt-srna-table", [h("thead", "", tr), tb]);
  }
  // the exported figure carries no hit layer and no hover state
  function exportCopy(svg) {
    var c = svg.cloneNode(true);
    function strip(el) {
      for (var i = el.children.length - 1; i >= 0; i--) {
        var k = el.children[i];
        if (!k.classList) { continue; }
        if (k.classList.contains("rpt-srna-hit")) { el.removeChild(k); continue; }
        k.classList.remove("cur");
        strip(k);
      }
    }
    c.classList.remove("lift");
    strip(c);
    return c;
  }
  function download(svg, name) {
    var s = new window.XMLSerializer().serializeToString(exportCopy(svg));
    var url = window.URL.createObjectURL(new window.Blob([s], {type: "image/svg+xml"}));
    var a = h("a");
    a.setAttribute("href", url);
    a.setAttribute("download", name);
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.setTimeout(function () { window.URL.revokeObjectURL(url); }, 1000);
  }
  // one observer for the overview charts and the open cluster page
  function watch() {
    if (observer) { observer.disconnect(); }
    observer = null;
    var all = ovLive.concat(live);
    if (window.ResizeObserver) {
      observer = new window.ResizeObserver(function (entries) {
        entries.forEach(function (en) {
          all.forEach(function (f) { if (f.plot === en.target) { f.render(); } });
        });
      });
      all.forEach(function (f) { observer.observe(f.plot); });
    }
  }
  window.addEventListener("resize", function () {
    placeIndicator();
    if (!window.ResizeObserver) { ovLive.concat(live).forEach(function (f) { f.render(); }); }
  });

  function strandLabel(s) { return T(s === "+" ? "ui.strand.sense" : "ui.strand.antisense"); }
  function strandLegend() {
    return [{kind: "sw", color: STRAND_FILL["+"], label: strandLabel("+")},
      {kind: "sw", color: STRAND_FILL["-"], label: strandLabel("-")}];
  }
  function baseLegend() {
    return BASES.map(function (b) { return {kind: "sw", color: BASE_FILL[b], label: b}; });
  }
  // two stacked panels, each titled above its top tick, one x axis title at the bottom
  function twoPanels(H) {
    var ph = (H - 132) / 2, p1 = 28, p2 = p1 + ph + 58;
    return [{sub: p1 - 14, y0: p1, y1: p1 + ph, tk: p1 + ph + 18},
      {sub: p2 - 14, y0: p2, y1: p2 + ph, tk: p2 + ph + 18, xl: H - 6}];
  }

  // ---- decoded signals of one cluster ----
  function decode(sig) {
    var lo = sig.sz.lo, nL = 35 - lo + 1, sizes = {}, z = {};
    BASES.forEach(function (b) {
      sizes[b] = {"+": deq([sig.sz.max, sig.sz[b][0]], nL, 999), "-": deq([sig.sz.max, sig.sz[b][1]], nL, 999)};
    });
    STRANDS.forEach(function (s) {
      z[s] = {ph: deq(sig.ph[s], 100), cv: deq(sig.cv[s], sig.nb), en: deq(sig.en[s], sig.nb)};
      z[s].reads = (sig.en[s] && sig.en[s][0] > 0) || (sig.cv[s] && sig.cv[s][0] > 0);
    });
    return {lo: lo, nL: nL, sizes: sizes, s: z, pp: deq(sig.pp, 30)};
  }
  function strandNote(dd, s) { return dd.s[s].reads ? "ui.fig.nodata" : "ui.fig.nostrand"; }

  // A: size x strand x 5' base, stacked, sense up and antisense down
  function figSizes(sig, dd, r) {
    var tot = {"+": [], "-": []}, i;
    for (i = 0; i < dd.nL; i++) {
      STRANDS.forEach(function (s) {
        tot[s].push(sum(BASES.map(function (b) { return dd.sizes[b][s][i]; })));
      });
    }
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, y0 = 24, y1 = H - 44, mid = (y0 + y1) / 2, half = mid - y0;
      var tk = ticks(Math.max(maxOf(tot["+"]), maxOf(tot["-"])), 3), sc = half / tk.top;
      yGrid(svg, x0, x1, mid, half, -1, tk);
      yGrid(svg, x0, x1, mid, half, 1, tk);
      var cols = slots(x0, x1, dd.nL), slot = (x1 - x0) / dd.nL, bw = Math.min(BAR_MAX, slot * 0.72);
      var tips = [];
      for (var j = 0; j < dd.nL; j++) {
        var bx = cols[j][0] + (slot - bw) / 2;
        STRANDS.forEach(function (s) {
          // U nearest the baseline, A at the tip (calcDistributionPerBase_publication.R)
          var parts = ["U", "G", "C", "A"].map(function (b) {
            return {v: dd.sizes[b][s][j], fill: BASE_FILL[b], letter: b, ink: BASE_INK[b]};
          });
          stack(svg, bx, bw, mid, s === "+" ? -1 : 1, parts, sc, j);
        });
      }
      zeroLine(svg, x0, x1, mid);
      STRANDS.forEach(function (s) {
        if (!(maxOf(tot[s]) > 0)) { emptyNote(svg, x0, x1, s === "+" ? mid - half / 2 : mid + half / 2, strandNote(dd, s)); }
      });
      var k = every(slot, 2), lab = [];
      for (j = 0; j < dd.nL; j++) {
        if ((dd.lo + j) % k === 0) { lab.push([cols[j][0] + slot / 2, String(dd.lo + j)]); }
      }
      xTicks(svg, lab, y1 + 18);
      // the dominant size, labelled once
      var di = sig.dom - dd.lo;
      if (typeof sig.dom === "number" && di >= 0 && di < dd.nL) {
        var cx = cols[di][0] + slot / 2, up = tot["+"][di] * sc, dn = tot["-"][di] * sc;
        var y = up > 0 || dn > half - 22 ? mid - up - 6 : mid + dn + 16;
        svg.appendChild(text(cx, y, T("ui.unit.nt", {v: sig.dom}), {anchor: "middle", fill: INK, weight: 500, cls: "dl"}));
      }
      xTitle(svg, T("ui.axis.length"), W, H - 6);
      yTitle(svg, T("ui.axis.reads"), mid);
      hits(svg, cols, y0, y1, function (j) { return tips[j] || (tips[j] = sizeTip(j)); });
    }
    function sizeTip(j) {
      var rowsT = [];
      STRANDS.forEach(function (s) {
        rowsT.push({key: null, value: val(tot[s][j]), label: strandLabel(s)});
        rowsT.push({label: BASES.map(function (b) { return b + " " + val(dd.sizes[b][s][j]); }).join("  ")});
      });
      return {head: T("ui.unit.nt", {v: dd.lo + j}), rows: rowsT};
    }
    return figure({letter: "A", title: T("ui.fig.sizes"), n: T("ui.fig.n", {n: num(sig.n)}),
      tip: T("ui.fig.sizes.tip"), legend: baseLegend().concat([{kind: "none", label: T("ui.fig.mirror")}]),
      height: 380, draw: draw, file: fileName(r, "A"),
      table: function () {
        var out = [];
        for (var j = 0; j < dd.nL; j++) {
          STRANDS.forEach(function (s) {
            out.push([String(dd.lo + j), s].concat(BASES.map(function (b) { return val(dd.sizes[b][s][j]); }),
              [val(tot[s][j])]));
          });
        }
        return {head: [T("ui.axis.length"), T("ui.th.strand")].concat(BASES, [T("ui.th.total")]), rows: out};
      }});
  }

  // B: phasing, one panel per strand, red dashed lines at multiples of the period
  function figPhasing(sig, dd, r) {
    var period = PERIODS[r.c] || (typeof sig.dom === "number" ? sig.dom : null);
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, P = twoPanels(H);
      var tk = ticks(Math.max(maxOf(dd.s["+"].ph), maxOf(dd.s["-"].ph)), 3);
      function xa(d) { return x0 + (d - 1) / 99 * (x1 - x0); }
      STRANDS.forEach(function (s, k) {
        var p = P[k], len = p.y1 - p.y0;
        svg.appendChild(text(x0, p.sub, strandLabel(s), {fill: INK, weight: 500, cls: "sub"}));
        yGrid(svg, x0, x1, p.y1, len, -1, tk);
        zeroLine(svg, x0, x1, p.y1);
        if (period) {
          for (var m = period; m <= 100; m += period) {
            svg.appendChild(sv("line", {x1: f2(xa(m)), x2: f2(xa(m)), y1: f2(p.y0), y2: f2(p.y1), stroke: PERIOD_RED,
              "stroke-width": 1.5, "stroke-dasharray": "4 4", "class": "period"}));
          }
        }
        var v = dd.s[s].ph;
        if (dd.s[s].reads) {
          var d = "";
          for (var i = 0; i < 100; i++) { d += (i ? "L" : "M") + f2(xa(i + 1)) + "," + f2(p.y1 - v[i] / tk.top * len); }
          svg.appendChild(sv("path", {d: d, fill: "none", stroke: STRAND_FILL[s], "stroke-width": 2,
            "stroke-linejoin": "round", "stroke-linecap": "round", "class": "mk-line"}));
        } else {
          emptyNote(svg, x0, x1, (p.y0 + p.y1) / 2, "ui.fig.nostrand");
        }
        xTicks(svg, [1, 20, 40, 60, 80, 100].map(function (t) { return [xa(t), String(t)]; }), p.tk);
      });
      xTitle(svg, T("ui.axis.distance"), W, P[1].xl);
      P.forEach(function (p) { yTitle(svg, T("ui.axis.pairs"), (p.y0 + p.y1) / 2); });
      var w = (x1 - x0) / 99, cols = [];
      for (var c = 1; c <= 100; c++) { cols.push([xa(c) - w / 2, xa(c) + w / 2]); }
      hits(svg, cols, P[0].y0, P[1].y1, function (i) {
        return {head: T("ui.tip.dist", {x: i + 1}), rows: STRANDS.map(function (s) {
          return {key: STRAND_FILL[s], value: val(dd.s[s].ph[i]), label: strandLabel(s)};
        })};
      });
    }
    var leg = [{kind: "ln", color: STRAND_FILL["+"], label: strandLabel("+")},
      {kind: "ln", color: STRAND_FILL["-"], label: strandLabel("-")}];
    if (period) { leg.push({kind: "dash", color: PERIOD_RED, label: T("ui.fig.period", {p: period})}); }
    return figure({letter: "B", title: T("ui.fig.phasing"), n: T("ui.fig.n", {n: num(sig.n)}),
      tip: T("ui.fig.phasing.tip"), legend: leg, height: 380, draw: draw, file: fileName(r, "B"),
      table: function () {
        var out = [];
        for (var i = 0; i < 100; i++) { out.push([String(i + 1), val(dd.s["+"].ph[i]), val(dd.s["-"].ph[i])]); }
        return {head: [T("ui.axis.distance"), strandLabel("+"), strandLabel("-")], rows: out};
      }});
  }

  // C and E: a binned track along the whole cluster, sense up and antisense down
  function trackFigure(o, sig, dd, r) {
    var nb = sig.nb, b = o.bin, s0 = o.start;
    // the cluster end clips the last bin; the bins win when the two disagree
    var full = s0 + nb * b;
    var e0 = typeof r.e === "number" && r.e > s0 + (nb - 1) * b && r.e < full ? r.e : full;
    function binEnd(i) { return Math.min(s0 + (i + 1) * b, e0); }
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, y0 = 12, y1 = H - 44, mid = (y0 + y1) / 2, half = mid - y0;
      var tk = ticks(Math.max(maxOf(dd.s["+"][o.key]), maxOf(dd.s["-"][o.key])), 2);
      function xa(p) { return x0 + (p - s0) / Math.max(1, e0 - s0) * (x1 - x0); }
      yGrid(svg, x0, x1, mid, half, -1, tk);
      yGrid(svg, x0, x1, mid, half, 1, tk);
      var cols = [];
      for (var i = 0; i < nb; i++) { cols.push([xa(s0 + i * b), xa(binEnd(i))]); }
      STRANDS.forEach(function (s) {
        var v = dd.s[s][o.key], dir = s === "+" ? -1 : 1, d = "";
        for (var j = 0; j < nb; j++) {
          if (!(v[j] > 0)) { continue; }
          var xa0 = cols[j][0], xb = cols[j][1], gap = o.gap && xb - xa0 >= 3 ? 1 : 0;
          var yv = mid + dir * v[j] / tk.top * half;
          d += "M" + f2(xa0) + "," + f2(mid) + "L" + f2(xa0) + "," + f2(yv) + "L" + f2(xb - gap) + "," + f2(yv) +
            "L" + f2(xb - gap) + "," + f2(mid) + "Z";
        }
        if (d) {
          svg.appendChild(sv("path", {d: d, fill: STRAND_FILL[s], "class": "mk-area " + (dir < 0 ? "up" : "dn")}));
        } else {
          emptyNote(svg, x0, x1, mid + dir * half / 2, strandNote(dd, s));
        }
      });
      zeroLine(svg, x0, x1, mid);
      // genomic ticks, 1-based, spaced so labels never touch
      var room = Math.max(2, Math.floor((x1 - x0) / (num(e0).length * 8 + 16)));
      var st = niceStep(Math.max(1, e0 - s0), room), lab = [];
      for (var p = Math.ceil((s0 + 1) / st) * st; p <= e0; p += st) { lab.push([xa(p - 1), num(p)]); }
      xTicks(svg, lab, y1 + 18);
      xTitle(svg, T("ui.axis.position", {chrom: r.h}), W, H - 6);
      yTitle(svg, o.ylab, mid);
      hits(svg, cols, y0, y1, function (j) {
        return {head: T("ui.tip.pos", {start: num(s0 + j * b + 1), end: num(binEnd(j))}),
          rows: STRANDS.map(function (s) {
            return {key: STRAND_FILL[s], value: val(dd.s[s][o.key][j]), label: strandLabel(s)};
          })};
      });
    }
    return figure({letter: o.letter, title: o.title, n: T("ui.fig.n", {n: num(sig.n)}), tip: o.tip,
      legend: strandLegend().concat([{kind: "none", label: T("ui.fig.bin", {b: num(b)})}]),
      height: o.height, wide: o.wide, draw: draw, file: fileName(r, o.letter),
      table: function () {
        var out = [];
        for (var j = 0; j < nb; j++) {
          out.push([T("ui.tip.pos", {start: num(s0 + j * b + 1), end: num(binEnd(j))}),
            val(dd.s["+"][o.key][j]), val(dd.s["-"][o.key][j])]);
        }
        return {head: [T("ui.th.bin"), strandLabel("+"), strandLabel("-")], rows: out};
      }});
  }

  // D: ping-pong overlap 1..30, the 10 nt bar in red (plotPingPong.R)
  function figPingPong(sig, dd, r) {
    var animal = META.kingdom !== "plant";
    function fillOf(i) { return animal && i === 9 ? PP_RED : GREY60; }
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, y0 = 24, y1 = H - 44, len = y1 - y0;
      var tk = ticks(maxOf(dd.pp), 3), sc = len / tk.top;
      yGrid(svg, x0, x1, y1, len, -1, tk);
      var cols = slots(x0, x1, 30), slot = (x1 - x0) / 30, bw = Math.min(BAR_MAX, slot - GAP), peak = -1;
      dd.pp.forEach(function (v, i) {
        if (v > 0 && (peak < 0 || v > dd.pp[peak])) { peak = i; }
        if (v * sc >= 0.5) {
          svg.appendChild(tagMark(barPath(cols[i][0] + (slot - bw) / 2, y1, bw, v * sc, -1, "mk-bar", fillOf(i)), i));
        }
      });
      zeroLine(svg, x0, x1, y1);
      if (peak >= 0) {
        svg.appendChild(text(cols[peak][0] + slot / 2, y1 - dd.pp[peak] * sc - 6, val(dd.pp[peak]),
          {anchor: "middle", fill: INK, weight: 500, cls: "dl"}));
      } else {
        emptyNote(svg, x0, x1, (y0 + y1) / 2, "ui.fig.nopairs");
      }
      xTicks(svg, [1, 5, 10, 15, 20, 25, 30].map(function (t) { return [cols[t - 1][0] + slot / 2, String(t)]; }), y1 + 18);
      xTitle(svg, T("ui.axis.overlap"), W, H - 6);
      yTitle(svg, T("ui.axis.pairs"), (y0 + y1) / 2);
      hits(svg, cols, y0, y1, function (i) {
        return {head: T("ui.tip.overlap", {x: i + 1}), rows: [{key: fillOf(i), value: val(dd.pp[i]), label: T("ui.axis.pairs")}]};
      });
    }
    return figure({letter: "D", title: KT.pp_title || T("ui.fig.pp"), n: T("ui.fig.n", {n: num(sig.n)}),
      tip: T("ui.fig.pp.tip"),
      legend: animal ? [{kind: "sw", color: PP_RED, label: T("ui.fig.pp10")},
        {kind: "sw", color: GREY60, label: T("ui.fig.ppother")}] : [],
      height: 280, draw: draw, file: fileName(r, "D"),
      table: function () {
        return {head: [T("ui.axis.overlap"), T("ui.axis.pairs")],
          rows: dd.pp.map(function (v, i) { return [String(i + 1), val(v)]; })};
      }});
  }

  // F: base logo of the dominant size, one panel per strand, largest base on top
  function figLogo(sig, dd, r) {
    var dom = sig.dom;
    function frac(s, i) {
      var p = (sig.lg[s] || [])[i] || [0, 0, 0, 0];
      return BASES.map(function (b, k) { return (p[k] || 0) / 99; });
    }
    function empty(s) {
      for (var i = 0; i < 26; i++) { if (sum(frac(s, i)) > 0) { return false; } }
      return true;
    }
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT, P = twoPanels(H), cols = slots(x0, x1, 26), slot = (x1 - x0) / 26;
      var bw = Math.min(BAR_MAX, slot - GAP), tk = {top: 1, step: 0.5, list: [0, 0.5, 1]};
      STRANDS.forEach(function (s, k) {
        var p = P[k], len = p.y1 - p.y0;
        svg.appendChild(text(x0, p.sub, strandLabel(s), {fill: INK, weight: 500, cls: "sub"}));
        [0.5, 1].forEach(function (t) {
          var y = p.y1 - t * len;
          hline(svg, x0, x1, y, GRID, "grid");
          svg.appendChild(text(x0 - 8, y + 4.5, T("ui.unit.pct", {v: t * 100}), {anchor: "end", cls: "tk-y"}));
        });
        zeroLine(svg, x0, x1, p.y1);
        if (empty(s)) {
          emptyNote(svg, x0, x1, (p.y0 + p.y1) / 2, strandNote(dd, s));
        } else {
          for (var i = 0; i < 26; i++) {
            var f = frac(s, i);
            var parts = BASES.map(function (b, j) { return {v: f[j], fill: BASE_FILL[b], letter: b, ink: BASE_INK[b]}; });
            parts.sort(function (a, c) { return a.v - c.v; });
            stack(svg, cols[i][0] + (slot - bw) / 2, bw, p.y1, -1, parts, len / tk.top, i);
          }
        }
        var kk = every(slot, 2), lab = [];
        for (var q = 0; q < 26; q++) { if (q === 0 || (q + 1) % Math.max(5, kk) === 0) { lab.push([cols[q][0] + slot / 2, String(q + 1)]); } }
        xTicks(svg, lab, p.tk);
      });
      xTitle(svg, T("ui.axis.logo"), W, P[1].xl);
      P.forEach(function (p) { yTitle(svg, T("ui.axis.share"), (p.y0 + p.y1) / 2); });
      hits(svg, cols, P[0].y0, P[1].y1, function (i) {
        return {head: T("ui.tip.logo", {x: i + 1}), rows: STRANDS.map(function (s) {
          var f = frac(s, i);
          return {label: strandLabel(s) + "  " + BASES.map(function (b, j) {
            return b + " " + T("ui.unit.pct", {v: Math.round(f[j] * 100)});
          }).join("  ")};
        })};
      });
    }
    return figure({letter: "F", title: T("ui.fig.logo", {z: typeof dom === "number" ? dom : T("ui.na")}),
      n: T("ui.fig.n", {n: num(sig.nd)}), tip: T("ui.fig.logo.tip"), legend: baseLegend(),
      height: 340, wide: true, draw: draw, file: fileName(r, "F"),
      table: function () {
        var out = [];
        STRANDS.forEach(function (s) {
          for (var i = 0; i < 26; i++) {
            out.push([String(i + 1), s].concat(frac(s, i).map(function (x) { return T("ui.unit.pct", {v: Math.round(x * 100)}); })));
          }
        });
        return {head: [T("ui.axis.logo"), T("ui.th.strand")].concat(BASES), rows: out};
      }});
  }

  // G: the most abundant distinct sequences, exact
  function figSeqs(sig, r) {
    var sq = (sig.sq || []).slice(0, 10), body;
    if (!sq.length) {
      body = emptyMsg(T("ui.fig.noseqs"));
    } else {
      var head = h("tr", "", [T("ui.th.seq"), T("ui.th.length"), T("ui.th.strand"), T("ui.th.reads"), T("ui.th.base")].map(
        function (t, j) { return h("th", j === 1 || j === 3 ? "num" : "", t); }));
      var tb = h("tbody");
      sq.forEach(function (q) {
        var sw = h("span", "rpt-srna-sw");
        sw.style.background = BASE_FILL[q[4]] || GREY60;
        sw.setAttribute("aria-hidden", "true");
        tb.appendChild(h("tr", "", [h("td", "rpt-srna-seq", q[0]), h("td", "num", String(q[1])),
          h("td", "rpt-srna-mono", q[2] === "-" ? "−" : "+"), h("td", "num", exact(q[3])),
          h("td", "rpt-srna-mono", [sw, q[4]])]));
      });
      body = h("div", "rpt-srna-scroll", h("table", "rpt-srna-table rpt-srna-seqs", [h("thead", "", head), tb]));
    }
    return figure({letter: "G", title: T("ui.fig.seqs"), n: T("ui.fig.nseq", {n: sq.length}),
      tip: T("ui.fig.seqs.tip"), wide: true, body: body});
  }

  // H: where the cluster sits among the shown clusters (reads; confidence within its class)
  function figContext(d, r) {
    var ctx = d.ctx || [null, null];
    function hist(vals, lo, hi, n) {
      var c = [], w = (hi - lo) / n || 1;
      for (var i = 0; i < n; i++) { c.push(0); }
      vals.forEach(function (v) { c[Math.max(0, Math.min(n - 1, Math.floor((v - lo) / w)))] += 1; });
      return {c: c, lo: lo, w: w, at: function (v) { return Math.max(0, Math.min(n - 1, Math.floor((v - lo) / w))); }};
    }
    var lr = rows.filter(function (x) { return x.r > 0; }).map(function (x) { return Math.log(x.r) / Math.LN10; });
    var mn = Infinity, mx = -Infinity;
    for (var li = 0; li < lr.length; li++) {
      if (lr[li] < mn) { mn = lr[li]; }
      if (lr[li] > mx) { mx = lr[li]; }
    }
    var rlo = lr.length ? Math.floor(mn * 4) / 4 : 0;
    var rhi = lr.length ? Math.ceil(mx * 4 + 1e-9) / 4 : 1;
    if (rhi <= rlo) { rhi = rlo + 0.25; }
    var HR = hist(lr, rlo, rhi, Math.min(24, Math.max(4, Math.round((rhi - rlo) * 8))));
    var confs = r.c === "ambiguous" ? [] : rows.filter(function (x) { return x.c === r.c && x.f !== null; })
      .map(function (x) { return x.f; });
    var HC = hist(confs, 0, 1, 20);
    var me = [r.r > 0 ? HR.at(Math.log(r.r) / Math.LN10) : -1, confs.length > 1 && r.f !== null ? HC.at(r.f) : -1];
    var top = ctx[0] === null ? null : Math.max(1, 100 - ctx[0]);
    function rlab(i) {
      return T("ui.tip.range", {lo: num(Math.pow(10, HR.lo + i * HR.w)), hi: num(Math.pow(10, HR.lo + (i + 1) * HR.w))});
    }
    function clab(i) {
      return T("ui.tip.range", {lo: (HC.lo + i * HC.w).toFixed(2), hi: (HC.lo + (i + 1) * HC.w).toFixed(2)});
    }
    // short panel titles; the full sentences live in the strip's hover tips
    var subs = [T("ui.ctx.sub", {lb: T("ui.ctx.reads"), v: top === null ? T("ui.na") : T("ui.ctx.top", {pct: top})}),
      T("ui.ctx.sub", {lb: T("ui.ctx.conf"), v: ctx[1] === null || ctx[1] === undefined ? T("ui.na") :
        T("ui.ctx.pct", {pct: ctx[1]})})];
    function draw(svg, W, H) {
      var x0 = LEFT, x1 = W - RIGHT;
      var P = [{sub: 14, y0: 28, y1: 118, tk: 136, xl: 158}, {sub: 184, y0: 198, y1: 288, tk: 306, xl: 328}];
      [[HR, me[0], function (i) { return rlab(i); }, T("ui.axis.logreads")],
        [HC, me[1], function (i) { return clab(i); }, T("ui.axis.conf")]].forEach(function (pn, k) {
        var Hh = pn[0], p = P[k], len = p.y1 - p.y0, n = Hh.c.length;
        svg.appendChild(text(x0, p.sub, subs[k], {fill: INK, weight: 500, cls: "sub"}));
        var tk = ticks(maxOf(Hh.c), 2);
        yGrid(svg, x0, x1, p.y1, len, -1, tk);
        zeroLine(svg, x0, x1, p.y1);
        var cols = slots(x0, x1, n), slot = (x1 - x0) / n, bw = Math.min(BAR_MAX, slot - GAP);
        if (k === 1 && confs.length < 2) {
          emptyNote(svg, x0, x1, (p.y0 + p.y1) / 2, "ui.ctx.conf.na");
        } else {
          Hh.c.forEach(function (v, i) {
            if (v > 0) {
              svg.appendChild(tagMark(barPath(cols[i][0] + (slot - bw) / 2, p.y1, bw, Math.max(1, v / tk.top * len), -1,
                "mk-bar", i === pn[1] ? ACCENT : OTHER_BAR), k + ":" + i));
            }
          });
        }
        var lab = [];
        if (k === 0) {
          lab = logTicks(HR.lo, HR.lo + n * HR.w, x0, x1);
        } else {
          lab = [0, 0.25, 0.5, 0.75, 1].map(function (t) { return [x0 + t * (x1 - x0), String(t)]; });
        }
        xTicks(svg, lab, p.tk);
        svg.appendChild(text(x0 + (x1 - x0) / 2, p.xl, pn[3], {anchor: "middle", cls: "xl"}));
        hits(svg, cols, p.y0, p.y1, function (i) {
          var rw = [{key: i === pn[1] ? ACCENT : OTHER_BAR, value: num(Hh.c[i]), label: T("ui.axis.clusters")}];
          if (i === pn[1]) { rw.push({label: T("ui.tip.this")}); }
          return {head: pn[2](i), rows: rw};
        }, function (i) { return k + ":" + i; });
      });
      P.forEach(function (p) { yTitle(svg, T("ui.axis.clusters"), (p.y0 + p.y1) / 2); });
    }
    return figure({letter: "H", title: T("ui.fig.ctx"), n: T("ui.fig.nclusters", {n: num(rows.length)}),
      tip: T("ui.fig.ctx.tip"),
      legend: [{kind: "sw", color: ACCENT, label: T("ui.tip.this")}, {kind: "sw", color: OTHER_BAR, label: T("ui.fig.others")}],
      height: 340, wide: true, draw: draw, file: fileName(r, "H"),
      table: function () {
        var out = [];
        HR.c.forEach(function (v, i) { out.push([T("ui.axis.logreads"), rlab(i), num(v)]); });
        if (confs.length > 1) { HC.c.forEach(function (v, i) { out.push([T("ui.axis.conf"), clab(i), num(v)]); }); }
        return {head: ["", T("ui.th.bin"), T("ui.axis.clusters")], rows: out, exact: true};
      }});
  }

  // 1-2-5 ticks on a log10 axis, thinned so labels keep 64 px apart
  function logTicks(lo, hi, x0, x1) {
    var out = [], last = -1e9;
    for (var e = Math.floor(lo); e <= Math.ceil(hi); e++) {
      [1, 2, 5].forEach(function (m) {
        var v = Math.log(m) / Math.LN10 + e;
        if (v < lo - 1e-9 || v > hi + 1e-9) { return; }
        var x = x0 + (v - lo) / (hi - lo) * (x1 - x0);
        if (x - last >= 64) { out.push([x, num(m * Math.pow(10, e))]); last = x; }
      });
    }
    return out;
  }
  // ---- genomic context: coverage and annotation tracks on one coordinate axis ----
  var GM_GUT = 112, GM_RIGHT = 16;
  var GM_TRACKS = ["g", "t", "r", "p"];
  var GM_LAYER = {g: "gene", t: "te", r: "ncrna", p: "phas"};
  var GM_TITLE = {cv: "ui.gm.track.cv", ends: "ui.gm.track.ends", g: "ui.gm.track.g", t: "ui.gm.track.t",
    r: "ui.gm.track.r", p: "ui.gm.track.p"};
  var GM_KIND = {g: "ui.gm.kind.g"};
  var GM_CAP = {g: 10, t: 14, r: 6, p: 6};
  var GM_LANE = {g: 26, t: 12, r: 20, p: 20};
  var GM_SIZE_FILL = {"18-20": "#94a3b8", "21": "#12a3a3", "22": "#e0921c", "23": "#1f9d6b",
    "24": "#2c5f8d", "25-35": "#7c4dff"};
  var GM_TE = ["LTR", "LINE", "SINE", "DNA", "RC"];
  var GM_TE_FILL = {LTR: "#c46a2a", LINE: "#6b4d99", SINE: "#a8413a", DNA: "#2c5f8d", RC: "#137d6f", other: "#8493a3"};
  var GM_FILL = {r: "#a0527c", p: "#1f9d6b"};
  var GM_GENE = "#25527d", GM_GENE_LIGHT = "#8eaecb", GM_PREC = "#15314b", GM_BAND = "#eaf3fb";
  var GM_CHAR = 8;

  // the part of a family before "/", one of the five classes or "other"
  function gmTeClass(family) {
    var c = String(family || "").split("/")[0].toUpperCase();
    return GM_TE.indexOf(c) >= 0 ? c : "other";
  }
  // a feature span cut to the window; l / r say which end was cut
  function gmClip(s, e, ws, we) {
    return {s: Math.max(s, ws), e: Math.min(e, we), l: s < ws, r: e > we};
  }
  // greedy lane packing: items [{a, b}] (any order); lane index per item, boxes in a lane
  // keep `gap` apart. `cap` limits the lanes: items that fit in none get lane -1.
  function gmLanes(items, gap, cap, first) {
    var order = items.map(function (it, i) { return i; });
    order.sort(function (p, q) {
      var fp = first && first(items[p]) ? 0 : 1, fq = first && first(items[q]) ? 0 : 1;
      return cap && fp !== fq ? fp - fq : items[p].a - items[q].a || items[p].b - items[q].b;
    });
    var ends = [], lane = items.map(function () { return -1; });
    order.forEach(function (i) {
      var k = 0;
      while (k < ends.length && items[i].a < ends[k] + gap) { k++; }
      if (cap && k >= cap) { return; }
      ends[k] = items[i].b;
      lane[i] = k;
    });
    return {lane: lane, n: ends.length};
  }
  // rows of one track as objects
  function gmRows(kind, list) {
    return (list || []).map(function (x) {
      if (kind === "g") { return {kind: kind, label: x[0], strand: x[1], s: x[2], e: x[3], hl: !!x[4], blocks: x[5] || [], type: T(GM_KIND.g)}; }
      if (kind === "t") { return {kind: kind, label: x[0], type: x[1], strand: x[2], s: x[3], e: x[4], hl: !!x[5], cls: gmTeClass(x[1])}; }
      return {kind: kind, label: x[0], type: x[1], strand: x[2], s: x[3], e: x[4], hl: !!x[5]};
    });
  }
  function gmStranded(s) { return s === "+" || s === "-"; }
  // widest stretch of a gene between its exon boxes (pixels, clipped to the window): [x, width]
  function gmGap(rw, a, b, xa, ws, we) {
    var bx = [], best = [(a + b) / 2, 0], at = a;
    rw.blocks.forEach(function (k) {
      var c = gmClip(k[1], k[2], ws, we);
      if (c.e > c.s) { bx.push([xa(c.s), Math.max(xa(c.e), xa(c.s) + 1.5)]); }
    });
    bx.sort(function (p, q) { return p[0] - q[0]; });
    bx.concat([[b, b]]).forEach(function (k) {
      if (k[0] - at > best[1]) { best = [(at + k[0]) / 2, k[0] - at]; }
      at = Math.max(at, k[1]);
    });
    return best;
  }
  // pixel boxes, label places and lanes of one track (pure: xa maps a position to x)
  function gmPlace(kind, rows, xa, x0, x1, ws, we) {
    var items = rows.map(function (rw) {
      var c = gmClip(rw.s, rw.e, ws, we), a = xa(c.s), b = xa(c.e), cx = (a + b) / 2;
      if (b - a < 2) { a = cx - 1; b = cx + 1; }
      var it = {row: rw, clip: c, a: a, b: b, label: null, ea: a, eb: b,
        // the arrow only where at least a pixel of the feature is cut off
        cl: c.l && (ws - rw.s) * (x1 - x0) / Math.max(1, we - ws) >= 1,
        cr: c.r && (rw.e - we) * (x1 - x0) / Math.max(1, we - ws) >= 1};
      var tw = (rw.label ? String(rw.label).length : 0) * GM_CHAR;
      if (kind !== "t" && tw) {
        if (b + 6 + tw <= x1) { it.label = {x: b + 6, anchor: "start", inside: false}; it.eb = b + 6 + tw; }
        else if (a - 6 - tw >= x0) { it.label = {x: a - 6, anchor: "end", inside: false}; it.ea = a - 6 - tw; }
        else if (kind === "g") {
          var gp = gmGap(rw, a, b, xa, ws, we);
          if (gp[1] >= tw + 14) { it.label = {x: gp[0], anchor: "middle", inside: true}; }
        } else if (b - a >= tw + 16) { it.label = {x: (a + b) / 2, anchor: "middle", inside: true}; }
      }
      return it;
    });
    var cap = GM_CAP[kind];
    var ext = items.map(function (it) { return {a: it.ea, b: it.eb, hl: it.row.hl}; });
    var pk = gmLanes(ext, 8, 0), dropped = 0;
    if (pk.n > cap) { pk = gmLanes(ext, 8, cap, function (e) { return e.hl; }); }
    items.forEach(function (it, i) { it.lane = pk.lane[i]; if (it.lane < 0) { dropped++; } });
    return {items: items, lanes: Math.min(pk.n, cap), dropped: dropped};
  }
  // coverage: every series [max, ints 0..99] to values in one unit; the stacks share one scale
  function gmCoverage(cv) {
    var nb = cv.nb || 0, vals = {"+": [], "-": []}, tot = {"+": [], "-": []}, shared = 0, i;
    STRANDS.forEach(function (s) {
      for (i = 0; i < nb; i++) { tot[s].push(0); }
      (cv.k || []).forEach(function (k) {
        var pr = cv[s] ? cv[s][k] : null, v = deq(pr, nb, 99);
        if (pr && pr[0] > shared) { shared = pr[0]; }
        vals[s].push(v);
        for (i = 0; i < nb; i++) { tot[s][i] += v[i]; }
      });
    });
    var ends = {"+": deq(cv.e && cv.e["+"], nb, 99), "-": deq(cv.e && cv.e["-"], nb, 99)};
    var eShared = Math.max(cv.e && cv.e["+"] ? cv.e["+"][0] : 0, cv.e && cv.e["-"] ? cv.e["-"][0] : 0);
    return {vals: vals, tot: tot, ends: ends, shared: shared, eShared: eShared,
      top: Math.max(maxOf(tot["+"]), maxOf(tot["-"])), eTop: Math.max(maxOf(ends["+"]), maxOf(ends["-"]))};
  }
  // words onto lines of at most `max` characters (a longer word stays whole)
  function gmWrap(s, max) {
    var out = [], cur = "";
    String(s || "").split(" ").forEach(function (w) {
      if (!w) { return; }
      if (cur && (cur + " " + w).length > max) { out.push(cur); cur = w; } else { cur = cur ? cur + " " + w : w; }
    });
    if (cur) { out.push(cur); }
    return out;
  }
  function gmStrandText(s) { return gmStranded(s) ? (s === "+" ? "+" : "−") : T("ui.gm.tip.nostrand"); }
  function gmCoords(s, e) { return T("ui.tip.pos", {start: num(s + 1), end: num(e)}); }
  function gmTip(rw) {
    var out = [{value: T("ui.gm.tip.type"), label: rw.type || ""},
      gmStranded(rw.strand) ? {value: T("ui.gm.tip.strand"), label: gmStrandText(rw.strand)} : {label: T("ui.gm.tip.nostrand")},
      {value: T("ui.gm.tip.pos"), label: gmCoords(rw.s, rw.e)}];
    if (rw.hl) { out.push({label: T("ui.gm.tip.overlap")}); }
    return {head: rw.label || "", rows: out};
  }
  function gmChev(x, y, dir, sz, color, cls) {
    return sv("polyline", {points: f2(x - dir * sz / 2) + "," + f2(y - sz) + " " + f2(x + dir * sz / 2) + "," + f2(y) + " " +
      f2(x - dir * sz / 2) + "," + f2(y + sz), fill: "none", stroke: color, "stroke-width": 1.2,
      "stroke-linecap": "round", "stroke-linejoin": "round", "class": cls || "rpt-gm-chev", "pointer-events": "none"});
  }
  // a small triangle at a window edge says the feature goes on
  function gmCont(x, y, dir, color) {
    return sv("path", {d: "M" + f2(x) + "," + f2(y) + "L" + f2(x + dir * 7) + "," + f2(y - 5) + "L" + f2(x + dir * 7) + "," + f2(y + 5) + "Z",
      fill: SURF, stroke: color, "stroke-width": 1.2, "stroke-linejoin": "round", "class": "rpt-gm-cont", "pointer-events": "none"});
  }
  function gmHalo(t) {
    t.setAttribute("stroke", SURF);
    t.setAttribute("stroke-width", "4");
    t.setAttribute("paint-order", "stroke");
    t.setAttribute("stroke-linejoin", "round");
    return t;
  }
  function gmFeatLabel(it, cy) {
    var l = it.label, t = text(l.x, cy + 4.5, it.row.label, {anchor: l.anchor, fill: INK, cls: "rpt-gm-lb"});
    return gmHalo(t);
  }
  // one gene model: exon boxes (UTR and plain exons low, CDS tall), introns as a line with chevrons
  function gmGene(it, cy, x0, x1, ws, we, xa) {
    var rw = it.row, g = sv("g", {"class": "rpt-gm-feat rpt-gm-gene" + (rw.hl ? " hl" : ""), "data-track": "g",
      "data-clip": (it.clip.l ? "l" : "") + (it.clip.r ? "r" : ""), opacity: rw.hl ? 1 : 0.82});
    var coding = rw.blocks.some(function (b) { return b[0] === "C"; });
    var H = {E: coding ? 8 : 12, U: 8, C: 16}, FILL = {E: coding ? GM_GENE_LIGHT : GM_GENE, U: GM_GENE_LIGHT, C: GM_GENE};
    var ink = rw.hl ? INK : "none", boxes = [];
    g.appendChild(sv("line", {x1: f2(it.a), x2: f2(it.b), y1: f2(cy), y2: f2(cy), stroke: GM_GENE, "stroke-width": 1,
      "class": "rpt-gm-intron"}));
    // exons and UTRs first, the CDS on top of them
    rw.blocks.slice().sort(function (p, q) { return (p[0] === "C") - (q[0] === "C"); }).forEach(function (b) {
      var c = gmClip(b[1], b[2], ws, we);
      if (c.e <= c.s) { return; }
      var a = xa(c.s), w = Math.max(1.5, xa(c.e) - a), h = H[b[0]] || H.E;
      boxes.push([a, a + w]);
      g.appendChild(sv("rect", {x: f2(a), y: f2(cy - h / 2), width: f2(w), height: h, rx: 1.5, fill: FILL[b[0]] || FILL.E,
        stroke: ink, "stroke-width": rw.hl ? 1.2 : 0, "class": "rpt-gm-exon", "data-b": b[0]}));
    });
    if (!rw.blocks.length) {
      g.appendChild(sv("rect", {x: f2(it.a), y: f2(cy - 5), width: f2(it.b - it.a), height: 10, rx: 1.5, fill: GM_GENE_LIGHT,
        stroke: ink, "stroke-width": rw.hl ? 1.2 : 0, "class": "rpt-gm-exon"}));
      boxes.push([it.a, it.b]);
    }
    if (gmStranded(rw.strand)) {
      var dir = rw.strand === "+" ? 1 : -1;
      for (var x = it.a + 8; x < it.b - 4; x += 16) {
        if (!boxes.some(function (bx) { return x > bx[0] - 4 && x < bx[1] + 4; })) { g.appendChild(gmChev(x, cy, dir, 3.5, GM_GENE)); }
      }
    }
    if (it.cl) { g.appendChild(gmCont(x0 + 1, cy, 1, GM_GENE)); }
    if (it.cr) { g.appendChild(gmCont(x1 - 1, cy, -1, GM_GENE)); }
    if (it.label) { g.appendChild(gmFeatLabel(it, cy)); }
    return g;
  }
  // a transposon, structural RNA or PHAS locus: a thin rounded bar
  function gmBar(it, cy, x0, x1, fill, bh) {
    var rw = it.row, g = sv("g", {"class": "rpt-gm-feat rpt-gm-bar" + (rw.hl ? " hl" : ""), "data-track": rw.kind,
      "data-clip": (it.clip.l ? "l" : "") + (it.clip.r ? "r" : ""), opacity: rw.hl ? 1 : 0.82});
    var w = it.b - it.a;
    g.appendChild(sv("rect", {x: f2(it.a), y: f2(cy - bh / 2), width: f2(w), height: bh, rx: Math.min(3, bh / 2, w / 2), fill: fill,
      stroke: rw.hl ? INK : "none", "stroke-width": rw.hl ? 1.2 : 0}));
    if (gmStranded(rw.strand) && w >= 18) {
      var dir = rw.strand === "+" ? 1 : -1;
      for (var x = it.a + 9; x < it.b - 5; x += 12) { g.appendChild(gmChev(x, cy, dir, 2.2, "#ffffff")); }
    }
    if (it.cl) { g.appendChild(gmCont(x0 + 1, cy, 1, fill)); }
    if (it.cr) { g.appendChild(gmCont(x1 - 1, cy, -1, fill)); }
    if (it.label) { g.appendChild(gmFeatLabel(it, cy)); }
    return g;
  }
  // pointer and keyboard targets over the features (arrow keys walk them left to right)
  function gmFeatHits(svg, list) {
    var els = [], on = -1;
    function mark(i) {
      if (on >= 0 && els[on]) { els[on].classList.remove("on"); }
      on = i;
      if (i >= 0) { els[i].classList.add("on"); }
    }
    function focusAt(j) {
      j = Math.max(0, Math.min(els.length - 1, j));
      els.forEach(function (e, k) { e.setAttribute("tabindex", k === j ? "0" : "-1"); });
      els[j].focus();
    }
    list.forEach(function (f, i) {
      var w = Math.max(8, f.w), e = sv("rect", {x: f2(f.x - (w - f.w) / 2), y: f2(f.y), width: f2(w), height: f2(f.h),
        fill: "transparent", "class": "rpt-srna-hit rpt-gm-hit", tabindex: i === 0 ? 0 : -1, role: "img", "data-track": f.track});
      e.setAttribute("aria-label", tipText(f.tip));
      e.addEventListener("pointermove", function (ev) { mark(i); showTip(f.tip, ev.clientX, ev.clientY); });
      e.addEventListener("pointerdown", function (ev) { mark(i); showTip(f.tip, ev.clientX, ev.clientY); });
      e.addEventListener("pointerleave", function (ev) {
        if (ev && ev.pointerType === "touch") { return; }
        mark(-1); hideTip();
      });
      e.addEventListener("focus", function () {
        var bb = e.getBoundingClientRect();
        mark(i); showTip(f.tip, bb.right, bb.top);
      });
      e.addEventListener("blur", function () { mark(-1); hideTip(); });
      e.addEventListener("keydown", function (ev) {
        var step = {ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1}[ev.key];
        if (ev.key === "Escape") { mark(-1); hideTip(); return; }
        if (step === undefined) { return; }
        if (ev.preventDefault) { ev.preventDefault(); }
        focusAt(i + step);
      });
      els.push(e);
      svg.appendChild(e);
    });
    return els;
  }

  function figGenome(d, r) {
    var gm = d.gm, ws = gm.w[0], we = gm.w[1], loc = gm.L || [ws, we], ly = gm.ly || [], mo = gm.mo || {};
    var cv = gm.cv || {nb: 0, k: [], s: ws, b: 1}, cov = gmCoverage(cv), nb = cv.nb || 0, bsz = cv.b || 1;
    var mi = d.mi && !d.mi.nowin && typeof d.mi.lo === "number" && d.mi.L > 0 ? d.mi : null;
    var tracks = GM_TRACKS.map(function (k) { return {kind: k, rows: gmRows(k, gm[k])}; });
    var present = {};
    tracks.forEach(function (t) { t.rows.forEach(function (rw) { if (rw.cls) { present[rw.cls] = true; } }); });
    var teKeys = GM_TE.concat(["other"]).filter(function (k) { return present[k]; });
    function teName(k) { return k === "other" ? T("ui.gm.te.other") : k; }
    function colorOf(rw) { return rw.kind === "t" ? GM_TE_FILL[rw.cls] : GM_FILL[rw.kind]; }

    function draw(svg, cw) {
      var W = parseFloat(svg.getAttribute("width")), x0 = GM_GUT, x1 = W - GM_RIGHT, span = Math.max(1, we - ws);
      function xa(p) { return x0 + (Math.min(Math.max(p, ws), we) - ws) / span * (x1 - x0); }
      // ---- plan: heights first, so the locus band can sit behind everything ----
      var y = 0, plan = [], RULER = 58, AXY = 28, HALF = 66, EH = 24, PREC = mi ? 20 : 0;
      plan.push({id: "ruler", y: 0, h: RULER}); y = RULER;
      plan.push({id: "cv", y: y, h: PREC + 8 + 2 * HALF + 8}); y += plan[1].h;
      plan.push({id: "ends", y: y, h: 8 + 2 * EH + 8}); y += plan[2].h;
      tracks.forEach(function (t) {
        var given = ly.indexOf(GM_LAYER[t.kind]) >= 0;
        if (!given && !t.rows.length) { t.state = "nolayer"; } else if (!t.rows.length) { t.state = "empty"; }
        t.place = gmPlace(t.kind, t.rows, xa, x0, x1, ws, we);
        var more = (mo[t.kind] || 0) + t.place.dropped;
        t.more = more;
        var lh = GM_LANE[t.kind], body = t.state ? 28 : t.place.lanes * lh + 12;
        if (more) { body = Math.max(body, 46); }
        t.y = y; t.h = body; t.lh = lh;
        plan.push({id: t.kind, y: y, h: body}); y += body;
      });
      var H = y + 4;
      // ---- locus band, vertical grid, separators ----
      var lx0 = xa(loc[0]), lx1 = Math.max(xa(loc[1]), lx0 + 2);
      svg.appendChild(sv("rect", {x: f2(lx0), y: AXY, width: f2(lx1 - lx0), height: H - AXY, fill: GM_BAND,
        "class": "rpt-gm-band"}));
      var room = Math.max(2, Math.floor((x1 - x0) / (num(we).length * 8 + 24)));
      var st = niceStep(Math.max(1, we - ws), room), lab = [];
      for (var p = Math.ceil((ws + 1) / st) * st; p <= we; p += st) { lab.push([xa(p - 1), num(p)]); }
      lab.forEach(function (t) {
        svg.appendChild(sv("line", {x1: f2(t[0]), x2: f2(t[0]), y1: AXY, y2: H - 4, stroke: GRID, "stroke-width": 1,
          "shape-rendering": "crispEdges", "class": "rpt-gm-vgrid", opacity: 0.7}));
      });
      plan.slice(1).forEach(function (pl) { hline(svg, 0, W, pl.y, GRID, "rpt-gm-sep"); });
      // ---- ruler ----
      hline(svg, x0, x1, AXY, AXIS, "rpt-gm-axis");
      lab.forEach(function (t) {
        svg.appendChild(sv("line", {x1: f2(t[0]), x2: f2(t[0]), y1: AXY - 5, y2: AXY, stroke: AXIS, "stroke-width": 1,
          "shape-rendering": "crispEdges"}));
      });
      xTicks(svg, lab, AXY - 11);
      svg.appendChild(text(8, AXY - 11, r.h, {cls: "rpt-gm-chrom", fill: INK, weight: 500}));
      [lx0, lx1].forEach(function (x) {
        svg.appendChild(sv("line", {x1: f2(x), x2: f2(x), y1: AXY, y2: H - 4, stroke: ACCENT, "stroke-width": 1,
          "shape-rendering": "crispEdges", "class": "rpt-gm-edge", opacity: 0.85}));
      });
      svg.appendChild(sv("rect", {x: f2(lx0), y: AXY - 1, width: f2(lx1 - lx0), height: 4, fill: ACCENT, "class": "rpt-gm-locus"}));
      var lw = T("ui.gm.locus").length * GM_CHAR, lt = text(lx1 + 6, AXY + 19, T("ui.gm.locus"), {fill: "#1a6cb8", cls: "rpt-gm-loclabel"});
      if (lx1 + 6 + lw > x1) { lt.setAttribute("x", f2(lx0 - 6)); lt.setAttribute("text-anchor", "end"); }
      svg.appendChild(lt);
      function tlabel(pl, title, sub) {
        svg.appendChild(text(8, pl.y + 21, title, {fill: INK, weight: 500, cls: "rpt-gm-tl"}));
        // the unit wraps into the gutter so it never meets the tick labels
        gmWrap(sub, 9).forEach(function (ln, i) { svg.appendChild(text(8, pl.y + 41 + i * 17, ln, {cls: "rpt-gm-ts"})); });
      }
      // ---- coverage ----
      var pc = plan[1], base = pc.y + PREC + 8 + HALF, tk = ticks(cov.top, 2), sc = HALF / tk.top;
      tlabel(pc, T(GM_TITLE.cv), T("ui.axis.depth"));
      yGrid(svg, x0, x1, base, HALF, -1, tk);
      yGrid(svg, x0, x1, base, HALF, 1, tk);
      function colX(i) { return [xa(cv.s + i * bsz), xa(Math.min(cv.s + (i + 1) * bsz, we))]; }
      STRANDS.forEach(function (s) {
        var dir = s === "+" ? -1 : 1;
        (cv.k || []).forEach(function (k, ki) {
          var dpath = "";
          for (var i = 0; i < nb; i++) {
            var v = cov.vals[s][ki][i];
            if (!(v > 0)) { continue; }
            var below = 0;
            for (var j = 0; j < ki; j++) { below += cov.vals[s][j][i]; }
            var c = colX(i), gp = c[1] - c[0] >= 3 ? 0.5 : 0, ya = base + dir * below * sc, yb = base + dir * (below + v) * sc;
            dpath += "M" + f2(c[0] + gp) + "," + f2(ya) + "L" + f2(c[0] + gp) + "," + f2(yb) + "L" + f2(c[1] - gp) + "," + f2(yb) +
              "L" + f2(c[1] - gp) + "," + f2(ya) + "Z";
          }
          if (dpath) { svg.appendChild(sv("path", {d: dpath, fill: GM_SIZE_FILL[k] || OTHER_BAR, "class": "rpt-gm-cov", "data-k": k, "data-s": s})); }
        });
        if (!(maxOf(cov.tot[s]) > 0)) { emptyNote(svg, x0, x1, base + dir * HALF / 2, "ui.fig.nostrand"); }
      });
      zeroLine(svg, x0, x1, base);
      if (nb) {
        var cols = [];
        for (var bi = 0; bi < nb; bi++) { cols.push(colX(bi)); }
        hits(svg, cols, base - HALF, base + HALF, function (i) {
          var rowsT = [];
          STRANDS.forEach(function (s) {
            rowsT.push({value: val(cov.tot[s][i]), label: strandLabel(s)});
            (cv.k || []).forEach(function (k, ki) {
              var v = cov.vals[s][ki][i];
              if (v > 0) { rowsT.push({key: GM_SIZE_FILL[k], value: val(v), label: T("ui.unit.nt", {v: k})}); }
            });
          });
          return {head: gmCoords(cv.s + i * bsz, Math.min(cv.s + (i + 1) * bsz, we)), rows: rowsT};
        });
      }
      if (mi) {
        var px0 = xa(mi.lo), px1 = Math.max(xa(mi.lo + mi.L), px0 + 3), py = pc.y + 11;
        var pg = sv("g", {"class": "rpt-gm-prec"});
        pg.appendChild(sv("line", {x1: f2(px0), x2: f2(px1), y1: py, y2: py, stroke: GM_PREC, "stroke-width": 3, "stroke-linecap": "round"}));
        [px0, px1].forEach(function (x) { pg.appendChild(sv("line", {x1: f2(x), x2: f2(x), y1: py - 4, y2: py + 4, stroke: GM_PREC, "stroke-width": 2})); });
        var pt = text(px1 + 8, py + 4.5, T("mi.title"), {fill: INK, cls: "rpt-gm-lb"});
        if (px1 + 8 + T("mi.title").length * GM_CHAR > x1) { pt.setAttribute("x", f2(px0 - 8)); pt.setAttribute("text-anchor", "end"); }
        pg.appendChild(pt);
        svg.appendChild(pg);
        gmFeatHits(svg, [{x: px0, y: py - 8, w: px1 - px0, h: 16, track: "mi",
          tip: {head: T("mi.title"), rows: [{value: T("ui.gm.tip.pos"), label: gmCoords(mi.lo, mi.lo + mi.L)}]}}]);
      }
      // ---- 5' ends ----
      var pe = plan[2], eb = pe.y + 8 + EH, ek = ticks(cov.eTop, 1), esc = EH / ek.top;
      tlabel(pe, T(GM_TITLE.ends));
      yGrid(svg, x0, x1, eb, EH, -1, ek);
      yGrid(svg, x0, x1, eb, EH, 1, ek);
      STRANDS.forEach(function (s) {
        var dir = s === "+" ? -1 : 1, dp = "";
        for (var i = 0; i < nb; i++) {
          var v = cov.ends[s][i];
          if (!(v > 0)) { continue; }
          var c = colX(i), gp = c[1] - c[0] >= 3 ? 0.5 : 0, yv = eb + dir * Math.max(1, v * esc);
          dp += "M" + f2(c[0] + gp) + "," + f2(eb) + "L" + f2(c[0] + gp) + "," + f2(yv) + "L" + f2(c[1] - gp) + "," + f2(yv) +
            "L" + f2(c[1] - gp) + "," + f2(eb) + "Z";
        }
        if (dp) { svg.appendChild(sv("path", {d: dp, fill: STRAND_FILL[s], "class": "rpt-gm-end", "data-s": s})); }
      });
      zeroLine(svg, x0, x1, eb);
      if (nb) {
        var ecols = [];
        for (var ei = 0; ei < nb; ei++) { ecols.push(colX(ei)); }
        hits(svg, ecols, eb - EH, eb + EH, function (i) {
          return {head: gmCoords(cv.s + i * bsz, Math.min(cv.s + (i + 1) * bsz, we)), rows: STRANDS.map(function (s) {
            return {key: STRAND_FILL[s], value: val(cov.ends[s][i]), label: strandLabel(s)};
          })};
        });
      }
      // ---- annotation tracks ----
      var fh = [];
      tracks.forEach(function (t) {
        tlabel(t, T(GM_TITLE[t.kind]));
        if (t.state) {
          svg.appendChild(gmHalo(text(x0 + 10, t.y + 21, t.state === "empty" ? T("ui.gm.empty") : T("ui.gm.nolayer", {kind: GM_LAYER[t.kind]}),
            {cls: "rpt-gm-empty"})));
        }
        t.place.items.forEach(function (it) {
          if (it.lane < 0) { return; }
          var cy = t.y + 6 + it.lane * t.lh + t.lh / 2, rw = it.row;
          svg.appendChild(t.kind === "g" ? gmGene(it, cy, x0, x1, ws, we, xa) :
            gmBar(it, cy, x0, x1, colorOf(rw), t.kind === "t" ? 8 : 12));
          fh.push({x: it.a, y: cy - t.lh / 2 + 1, w: it.b - it.a, h: t.lh - 2, track: t.kind, tip: gmTip(rw)});
        });
        if (t.more) { svg.appendChild(text(8, t.y + 41, T("ui.gm.more", {n: num(t.more)}), {cls: "rpt-gm-more"})); }
      });
      fh.sort(function (a, b) { return a.x - b.x; });
      if (fh.length) { gmFeatHits(svg, fh); }
      // the sized SVG
      svg.setAttribute("height", H);
      svg.setAttribute("viewBox", "0 0 " + W + " " + H);
      svg.firstChild.setAttribute("height", H);
      return H;
    }

    var legendItems = [{kind: "none", label: T("ui.fig.mirror")}];
    legendItems = legendItems.concat([{kind: "none", label: T("ui.gm.legend.size"), head: 1}],
      (cv.k || []).map(function (k) { return {kind: "sw", color: GM_SIZE_FILL[k] || OTHER_BAR, label: T("ui.unit.nt", {v: k})}; }),
      [{kind: "none", label: T("ui.gm.track.ends"), head: 1}], strandLegend());
    if (teKeys.length) {
      legendItems = legendItems.concat([{kind: "none", label: T("ui.gm.legend.te"), head: 1}],
        teKeys.map(function (k) { return {kind: "sw", color: GM_TE_FILL[k], label: teName(k)}; }));
    }
    return figure({letter: "I", title: T("ui.gm.title"),
      n: T("ui.gm.range", {chrom: r.h, start: num(loc[0] + 1), end: num(loc[1]), n: num(loc[1] - loc[0])}),
      legend: legendItems, height: 420, wide: true, width: 560, draw: draw, file: fileName(r, "I"),
      table: function () {
        var out = [];
        tracks.forEach(function (t) {
          t.rows.forEach(function (rw) {
            out.push([T(GM_TITLE[t.kind]), rw.label || "", rw.type || "", gmStrandText(rw.strand), num(rw.s + 1), num(rw.e),
              rw.hl ? "✓" : ""]);
          });
        });
        return {head: [T("ui.gm.th.track"), T("ui.gm.th.name"), T("ui.gm.tip.type"), T("ui.th.strand"), T("ui.gm.th.start"),
          T("ui.gm.th.end"), T("ui.gm.tip.overlap")], rows: out, plain: true, left: 4};
      }});
  }

  function fileName(r, letter) {
    var base = String(r.n).replace(/[^\w.-]+/g, "_").replace(/^_+|_+$/g, "").slice(0, 80);
    return (base || T("ui.brand")) + "_" + letter;
  }

  function figures(d, r) {
    var sig = d.sig;
    // the genomic context leads the page when the report was built with an annotation
    var gmFig = d.gm ? figGenome(d, r) : null;
    function lead(msg) { return gmFig ? h("div", "rpt-srna-figs", [gmFig, msg]) : msg; }
    if (!sig) { return lead(emptyMsg(T("ui.fig.na"))); }
    if (sig.err) { return lead(emptyMsg(T("ui.fig.err"))); }
    if (sig.none) { return lead(emptyMsg(T("ui.fig.noreads"))); }
    var dd = decode(sig);
    return h("div", "rpt-srna-figs", [
      gmFig,
      figSizes(sig, dd, r),
      figPhasing(sig, dd, r),
      trackFigure({letter: "C", key: "cv", bin: sig.cv.b, start: sig.cv.s, title: T("ui.fig.cov"),
        tip: T("ui.fig.cov.help"), ylab: T("ui.axis.depth"), height: 260, wide: true}, sig, dd, r),
      figPingPong(sig, dd, r),
      trackFigure({letter: "E", key: "en", bin: sig.en.b, start: sig.cv.s, title: T("ui.fig.ends"),
        tip: T("ui.fig.ends.tip"), ylab: T("ui.axis.ends"), height: 280, gap: true}, sig, dd, r),
      figLogo(sig, dd, r),
      figSeqs(sig, r),
      figContext(d, r)]);
  }

  // ---- cluster page: header, verdict, measured-value strip, figures ----
  function stripValue(k, v) {
    if (v === null || v === undefined) { return T("ui.na"); }
    if (k === "size") { return T("ui.unit.nt", {v: v}); }
    if (k === "strand") { return T("ui.unit.plus", {v: Math.round(v * 100)}); }
    if (k === "pingpong") { return T("ui.unit.z", {v: v.toFixed(1)}); }
    if (k === "phasing") { return T("ui.unit.p", {v: pval(v)}); }
    return pct(v);
  }
  function strip(d, r) {
    var box = h("div", "rpt-srna-strip");
    var def = (META.defining || {})[r.c] || [];
    (d.strip || []).forEach(function (it) {
      var k = it[0], lv = LEVEL[parseInt(it[1], 10)] || "unavailable";
      var defining = def.indexOf(k) >= 0;
      // checks that do not define the class are grey context; a value without a score is neutral
      var css = lv === "unavailable" ? "unavailable" : defining ? LEVEL_CSS[lv] : "context";
      var cell = h("div", "", [h("div", "lb", T("check." + k + ".title")),
        h("div", "vl", stripValue(k, it[2])), h("div", "lv " + css, T("status." + lv))]);
      if (!defining) { cell.setAttribute("title", T("check.context")); }
      box.appendChild(cell);
    });
    var ctx = d.ctx || [null, null];
    var rc = h("div", "", [h("div", "lb", T("ui.ctx.reads")),
      h("div", "vl", ctx[0] === null ? T("ui.na") : T("ui.ctx.top", {pct: Math.max(1, 100 - ctx[0])}))]);
    if (ctx[0] !== null) { rc.setAttribute("title", T("ui.tip.ctx", {pct: Math.max(1, 100 - ctx[0])})); }
    var cc = h("div", "", [h("div", "lb", T("ui.ctx.conf")),
      h("div", "vl", ctx[1] === null || ctx[1] === undefined ? T("ui.na") : T("ui.ctx.pct", {pct: ctx[1]}))]);
    cc.setAttribute("title", ctx[1] === null || ctx[1] === undefined ? T("ui.ctx.conf.na") :
      T("ui.tip.conf", {pct: ctx[1], cls: r.c}));
    add(box, [rc, cc]);
    return box;
  }
  // ---- class checklist: why the class was assigned (or withheld), one row per check ----
  // 16x16 glyphs, so the result never rests on colour alone: tick, cross, dash, ring
  var CK_PATH = {1: "M6.5 11.2L3.3 8l-1.1 1.1 4.3 4.3 8-8-1.1-1.1z",
    0: "M4 3L3 4l4 4-4 4 1 1 4-4 4 4 1-1-4-4 4-4-1-1-4 4z",
    2: "M3 7h10v2H3z",
    3: "M8 3a5 5 0 100 10A5 5 0 008 3zm0 1.5a3.5 3.5 0 110 7 3.5 3.5 0 010-7z"};
  function checklist(d, r) {
    var ck = d.ck, plan = ck && (META.cl || {})[ck[0]];
    if (r.c === "miRNA") {
      return h("section", "rpt-srna-ck", [h("h3", "rpt-srna-h", T("ck.title", {cls: r.c})), note(T("ck.mirna"))]);
    }
    if (!ck || !plan) { return null; }
    var head = h("tr", "", ["ck.th.check", "ck.th.observed", "ck.th.rule", "ck.th.status"].map(function (k) {
      return h("th", "", T(k));
    }));
    var tb = h("tbody");
    plan.forEach(function (p, i) {
      var row = ck[i + 1] || [2, ""], st = row[0];
      var name = h("th", "", T(p[0]));
      name.setAttribute("scope", "row");
      var gl = icon(CK_PATH[st] || CK_PATH[2], 14, "0 0 16 16");
      tb.appendChild(h("tr", "", [name, h("td", "rpt-srna-mono", row[1]), h("td", "", p[1]),
        h("td", "rpt-srna-ck__st s" + st, [gl, h("span", "", T("ck.st." + st))])]));
    });
    return h("section", "rpt-srna-ck", [
      h("h3", "rpt-srna-h", T(r.c === "ambiguous" ? "ck.title.amb" : "ck.title", {cls: ck[0]})),
      h("div", "rpt-srna-scroll", h("table", "rpt-srna-table", [h("thead", "", head), tb])),
      note(T("ck.note")), META.mode === "ml" ? note(T("ck.ml")) : null]);
  }
  // ---- miRNA hairpin: precursor, fold, reads and numbers (a miRNA cluster matched to --mirna) ----
  // Segment colours follow miRDeep2 (mature red, loop yellow, star purple), taken from the palette
  // already in the page: the ping-pong red, the warning amber and the piRNA purple. The solid
  // colour is only a band; letters sit on a light tint in dark ink and the segments are named in text.
  var MIR_FILL = {m: "#e41a1c", lp: "#e0921c", s: "#7c4dff"};
  var MIR_TINT = {m: "#fbd0d0", lp: "#fbe9c6", s: "#e4dcff"};
  var MIR_NAME = {m: "mi.seg.mature", lp: "mi.seg.loop", s: "mi.seg.star"};
  var MONO = "IBM Plex Mono, ui-monospace, monospace";
  var HP = {cell: 8, x0: 48, right: 196, row: 18, cov: 64};
  function hpX(p) { return HP.x0 + p * HP.cell; }
  function miSeg(mi, p) {
    var keys = ["m", "s", "lp"];
    for (var i = 0; i < keys.length; i++) {
      var g = mi[keys[i]];
      if (g && p >= g[0] && p < g[1]) { return keys[i]; }
    }
    return "";
  }
  // a run of letters locked to the cell grid whatever the font's own advance
  function mono(a, y, s, fill, weight) {
    return sv("text", {x: f2(hpX(a)), y: f2(y), textLength: s.length * HP.cell, lengthAdjust: "spacing",
      "font-family": MONO, fill: fill, "data-ink": "1", "font-weight": weight || null, "pointer-events": "none"}, s);
  }
  // partner index (0-based) of each base from the dot-bracket, -1 when unpaired
  function miPairs(db) {
    var p = [], st = [], i;
    for (i = 0; i < db.length; i++) { p.push(-1); }
    for (i = 0; i < db.length; i++) {
      var c = db.charAt(i);
      if (c === "(") { st.push(i); } else if (c === ")" && st.length) { var j = st.pop(); p[i] = j; p[j] = i; }
    }
    return p;
  }
  var D2 = {pad: 26, maxH: 420, maxScale: 2.4, minGap: 13};
  function miD2Box(mi) {
    if (!mi.db || !mi.xy || mi.xy.length !== 2 * mi.L) { return null; }
    var w = 0, h = 0, i;
    for (i = 0; i < mi.L; i++) { w = Math.max(w, mi.xy[2 * i]); h = Math.max(h, mi.xy[2 * i + 1]); }
    return {w: Math.max(w, 1), h: Math.max(h, 1)};
  }
  function miLayout(mi) {
    var L = {reads: !!(mi.rd && mi.rd.length)}, y = 6;
    var box = miD2Box(mi);
    if (box) {
      L.d2 = {top: y, h: Math.min(D2.maxH, Math.ceil(box.h * D2.maxScale) + 2 * D2.pad), box: box};
      y += L.d2.h + 12;
    }
    if (L.reads) { L.covTitle = y + 10; L.covTop = y + 18; L.covBase = L.covTop + HP.cov; y = L.covBase + 12; }
    L.lab = y + 12;
    L.band = y + 18;
    L.seq = y + 26;
    y = L.seq + 18;
    if (mi.db) { L.db = y + 15; y = L.db + 4; }
    L.tick = y + 15;
    y = L.tick + 8;
    L.rd = y;
    y += (L.reads ? mi.rd.length : 0) * HP.row;
    if (L.reads && mi.more) { L.more = y + 16; y += 24; }
    L.H = y + 8;
    return L;
  }
  // the classic 2D drawing: backbone, base-pair bonds, one circle per base. Scaled to the
  // card width (W) and the panel height, aspect kept. Letters only where the bases are >= 13 px apart.
  function miDrawD2(svg, W, mi, L) {
    var n = mi.L, b = L.d2.box, xy = mi.xy, pr = miPairs(mi.db), i;
    var sc = Math.min((W - 2 * D2.pad) / b.w, (L.d2.h - 2 * D2.pad) / b.h, D2.maxScale);
    var ox = (W - b.w * sc) / 2, oy = L.d2.top + (L.d2.h - b.h * sc) / 2;
    var px = [], py = [], gap = Infinity;
    for (i = 0; i < n; i++) { px.push(ox + xy[2 * i] * sc); py.push(oy + xy[2 * i + 1] * sc); }
    for (i = 1; i < n; i++) { gap = Math.min(gap, Math.hypot(px[i] - px[i - 1], py[i] - py[i - 1])); }
    var rad = Math.max(1.5, Math.min(9, gap * 0.42)), letters = gap >= D2.minGap;
    var pts = [];
    for (i = 0; i < n; i++) { pts.push(f2(px[i]) + "," + f2(py[i])); }
    svg.appendChild(sv("polyline", {points: pts.join(" "), fill: "none", stroke: AXIS, "stroke-width": 1,
      "class": "mk-backbone"}));
    var cx = 0, cy = 0;
    for (i = 0; i < n; i++) {
      cx += px[i] / n; cy += py[i] / n;
      if (pr[i] > i) {
        svg.appendChild(sv("line", {x1: f2(px[i]), y1: f2(py[i]), x2: f2(px[pr[i]]), y2: f2(py[pr[i]]),
          stroke: INK2, "stroke-width": 1.2, "class": "mk-bond"}));
      }
    }
    for (i = 0; i < n; i++) {
      var k = miSeg(mi, i), base = mi.seq.charAt(i);
      var tip = T("mi.d2.base", {pos: i + 1, base: base}) + (pr[i] >= 0 ? ", " + T("mi.d2.pair", {pos: pr[i] + 1}) : "") +
        (k ? ", " + T(MIR_NAME[k]) : "");
      svg.appendChild(sv("circle", {cx: f2(px[i]), cy: f2(py[i]), r: f2(rad), fill: k ? MIR_TINT[k] : SURF,
        stroke: k ? MIR_FILL[k] : INK2, "stroke-width": 1, "class": "mk-base"}, sv("title", {}, tip)));
      if (letters) {
        var fs = Math.min(12, rad * 1.45);
        svg.appendChild(sv("text", {x: f2(px[i]), y: f2(py[i] + fs * 0.36), "text-anchor": "middle", "font-size": f2(fs),
          "font-family": MONO, fill: INK, "data-ink": "1", "pointer-events": "none", "class": "mk-letter"}, base));
      }
    }
    // outward direction of base i, to put a label beside it
    function out(i) {
      var dx, dy;
      if (pr[i] >= 0) { dx = px[i] - px[pr[i]]; dy = py[i] - py[pr[i]]; } else {
        dx = px[i] - (px[Math.max(i - 1, 0)] + px[Math.min(i + 1, n - 1)]) / 2;
        dy = py[i] - (py[Math.max(i - 1, 0)] + py[Math.min(i + 1, n - 1)]) / 2;
        if (Math.hypot(dx, dy) < 1) { dx = px[i] - cx; dy = py[i] - cy; }
      }
      var m = Math.hypot(dx, dy) || 1;
      return [dx / m, dy / m];
    }
    for (i = 9; i < n; i += 10) {
      var d = out(i), o = rad + 8;
      svg.appendChild(text(px[i] + d[0] * o, py[i] + d[1] * o + 3, String(i + 1), {anchor: "middle", cls: "tk-x"}));
    }
    var e5 = [px[0] - px[1], py[0] - py[1]], e3 = [px[n - 1] - px[n - 2], py[n - 1] - py[n - 2]];
    [[e5, 0, "5′"], [e3, n - 1, "3′"]].forEach(function (q) {
      var m = Math.hypot(q[0][0], q[0][1]) || 1, o = rad + 11;
      svg.appendChild(text(px[q[1]] + q[0][0] / m * o, py[q[1]] + q[0][1] / m * o + 4, q[2], {anchor: "middle", fill: INK, weight: 500}));
    });
    if (mi.mfe !== undefined) {
      svg.appendChild(text(W - 6, L.d2.top + 12, T("mi.d2.mfe", {mfe: mi.mfe}), {anchor: "end"}));
    }
  }
  function miDraw(svg, W, mi, L) {
    var n = mi.L, xEnd = hpX(n), i, k;
    if (L.d2) { miDrawD2(svg, W, mi, L); }
    var cov = mi.cov || [], top = L.reads ? maxOf(cov) : 0;
    if (L.reads) {
      var tk = ticks(top, 3);
      svg.appendChild(text(HP.x0, L.covTitle, T("mi.axis.cov"), {}));
      yGrid(svg, HP.x0, xEnd, L.covBase, HP.cov, -1, tk);
      var d = "";
      for (i = 0; i < n; i++) {
        if (!(cov[i] > 0)) { continue; }
        var hgt = cov[i] / tk.top * HP.cov;
        d += "M" + f2(hpX(i)) + "," + f2(L.covBase) + "h" + (HP.cell - 1) + "v" + f2(-hgt) + "h" + (1 - HP.cell) + "Z";
      }
      if (d) { svg.appendChild(sv("path", {d: d, fill: ACCENT, "class": "mk-area up"})); }
      zeroLine(svg, HP.x0, xEnd, L.covBase);
    }
    // segment bands, their names and the tinted letters
    var segs = [["m", mi.m], ["lp", mi.lp], ["s", mi.s]];
    segs.forEach(function (sg) {
      var g = sg[1];
      if (!g) { return; }
      svg.appendChild(sv("rect", {x: f2(hpX(g[0])), y: L.band, width: (g[1] - g[0]) * HP.cell, height: 4,
        fill: MIR_FILL[sg[0]], "class": "mk-band"}));
      svg.appendChild(text(hpX(g[0]), L.lab, T(MIR_NAME[sg[0]]), {fill: INK}));
    });
    // runs of one segment kind, each with its tint
    var a = 0;
    while (a < n) {
      k = miSeg(mi, a);
      var b = a + 1;
      while (b < n && miSeg(mi, b) === k) { b++; }
      if (k) {
        svg.appendChild(sv("rect", {x: f2(hpX(a)), y: L.seq, width: (b - a) * HP.cell, height: 18,
          fill: MIR_TINT[k], "class": "mk-tint"}));
      }
      svg.appendChild(mono(a, L.seq + 13, mi.seq.slice(a, b), INK));
      a = b;
    }
    svg.appendChild(text(HP.x0 - 6, L.seq + 13, "5′", {anchor: "end"}));
    svg.appendChild(text(xEnd + 6, L.seq + 13, "3′", {}));
    if (mi.db) { svg.appendChild(mono(0, L.db, mi.db, INK2)); }
    for (i = 9; i < n; i += 10) {
      svg.appendChild(text(hpX(i) + HP.cell / 2, L.tick, String(i + 1), {anchor: "middle", cls: "tk-x"}));
    }
    if (!L.reads) { return; }
    var mx = 0;
    mi.rd.forEach(function (r) { if (r[3] > mx) { mx = r[3]; } });
    mi.rd.forEach(function (r, j) {
      var y = L.rd + j * HP.row, s = r[2], ra = r[0], rb = r[1];
      svg.appendChild(sv("rect", {x: f2(hpX(ra)), y: y + 1, width: (rb - ra) * HP.cell, height: HP.row - 3,
        fill: r[4] ? MIR_TINT[r[4]] : "#f1f5f9", "class": "mk-read"}));
      var c = 0;
      while (c < s.length) {
        var ch = s.charAt(c), mis = ch !== ch.toLowerCase();
        if (mis) {
          svg.appendChild(sv("rect", {x: f2(hpX(ra + c)), y: y + 1, width: HP.cell, height: HP.row - 3,
            fill: INK, "class": "mk-mis"}));
          svg.appendChild(sv("text", {x: f2(hpX(ra + c) + HP.cell / 2), y: y + 13, "text-anchor": "middle",
            "font-family": MONO, fill: SURF, "data-ink": "1", "font-weight": 500, "pointer-events": "none"}, ch));
          c++;
          continue;
        }
        var e = c + 1;
        while (e < s.length && s.charAt(e) === s.charAt(e).toLowerCase()) { e++; }
        svg.appendChild(mono(ra + c, y + 13, s.slice(c, e), INK2));
        c = e;
      }
      var bw = Math.max(2, Math.round(r[3] / mx * 54));
      svg.appendChild(sv("rect", {x: f2(xEnd + 10), y: y + 5, width: bw, height: 8, fill: INK2, "class": "mk-bar"}));
      svg.appendChild(text(xEnd + 72, y + 13, exact(r[3]), {}));
      if (r[4]) { svg.appendChild(text(xEnd + 136, y + 13, T(MIR_NAME[r[4]]), {fill: INK})); }
    });
    if (L.more) {
      svg.appendChild(text(HP.x0, L.more, T("mi.more", {n: num(mi.more), reads: num(mi.mr)}), {}));
    }
  }
  function miTable(mi) {
    var rows = [], n = mi.L, i, pr = mi.db ? miPairs(mi.db) : null;
    for (i = 0; i < n; i++) {
      var k = miSeg(mi, i);
      rows.push([String(i + 1), mi.seq.charAt(i), mi.db ? mi.db.charAt(i) : T("ui.na"),
        pr ? (pr[i] >= 0 ? String(pr[i] + 1) : "") : T("ui.na"),
        k ? T(MIR_NAME[k]) : "", mi.cov ? exact(mi.cov[i]) : T("ui.na")]);
    }
    var more = [];
    if (mi.rd && mi.rd.length) {
      var rr = mi.rd.map(function (r) {
        return [String(r[0] + 1), String(r[1]), r[2], exact(r[3]), r[4] ? T(MIR_NAME[r[4]]) : T("mi.kind.other")];
      });
      var head = ["mi.th.start", "mi.th.end", "mi.th.seq", "mi.th.reads", "mi.th.type"].map(function (t) { return T(t); });
      more.push(dataTable(head, rr));
    }
    return {exact: true, head: ["mi.th.pos", "mi.th.base", "mi.th.fold", "mi.th.pair", "mi.th.seg", "mi.th.depth"].map(function (t) { return T(t); }),
      rows: rows, more: more};
  }
  function miFigure(mi, r) {
    var L = miLayout(mi), tot = 0;
    (mi.rd || []).forEach(function (q) { tot += q[3]; });
    tot += mi.mr || 0;
    var legend = ["m", "lp", "s"].filter(function (k) {
      return k === "lp" ? mi.lp : mi[k];
    }).map(function (k) { return {kind: "sw", color: MIR_FILL[k], label: T(MIR_NAME[k])}; });
    if (L.reads) { legend.push({kind: "sw", color: INK, label: T("mi.lg.mismatch")}); }
    var el = figure({id: "hairpin", title: T("mi.title"), tip: T("mi.tip"), n: L.reads ? T("ui.fig.n", {n: num(tot)}) : null,
      legend: legend, height: L.H, width: hpX(mi.L) + HP.right, wide: true, file: fileName(r, "hairpin"),
      draw: function (svg, W) { miDraw(svg, W, mi, L); }, table: function () { return miTable(mi); }});
    var strand = T("mi.strand", {st: mi.st === "-" ? "−" : "+", chrom: r.h});
    var notes = [note(strand)];
    if (mi.db) {
      if (mi.mfe !== undefined) { notes.push(note(T("mi.struct", {mfe: mi.mfe}))); }
      if (L.d2) { notes.push(note(T("mi.d2.note"))); }
    } else {
      notes.push(note(T("mi.nostruct." + (mi.ns || "m"))));
    }
    if (mi.nowin) {
      notes.push(note(T("mi.nowin")));
    } else if (mi.nostrand) {
      notes.push(note(T("mi.nostrand")));
    } else if (mi.err) {
      notes.push(note(T("ui.fig.err")));
    } else if (!L.reads) {
      notes.push(note(T("mi.noreads")));
    } else {
      notes.push(note(T("mi.reads.note")));
      if (mi.anti || mi.ovh) { notes.push(note(T("mi.anti", {anti: num(mi.anti || 0), ovh: num(mi.ovh || 0)}))); }
    }
    add(el, notes);
    return el;
  }
  function miValue(k, v) {
    if (v === null || v === undefined) { return T("ui.na"); }
    if (typeof v === "string") { return v; }
    if (k === "precision") { return T("mi.obs.pct", {p: exact(v)}); }
    return k === "mature_reads" || k === "star_reads" ? exact(v) : String(v);
  }
  function miNumbers(mi) {
    var tb = h("tbody");
    mi.t.forEach(function (p) {
      var name = h("th", "", T("mi.num." + p[0]));
      name.setAttribute("scope", "row");
      tb.appendChild(h("tr", "", [name, h("td", "rpt-srna-mono", miValue(p[0], p[1]))]));
    });
    var head = h("tr", "", [h("th", "", T("ck.th.check")), h("th", "", T("ck.th.observed"))]);
    var crit = h("tbody");
    mi.pc.forEach(function (p) {
      var name = h("th", "", T("mi.crit." + p[0]));
      name.setAttribute("scope", "row");
      var st = p[1];
      crit.appendChild(h("tr", "", [name, h("td", "rpt-srna-mono", p[2]), h("td", "", p[3]),
        h("td", "rpt-srna-ck__st s" + st, [icon(CK_PATH[st] || CK_PATH[2], 14, "0 0 16 16"),
          h("span", "", T("ck.st." + st))])]));
    });
    var chead = h("tr", "", ["ck.th.check", "ck.th.observed", "ck.th.rule", "ck.th.status"].map(function (k) {
      return h("th", "", T(k));
    }));
    return h("div", "rpt-srna-fig--wide", [
      h("section", "rpt-srna-ck", [h("h3", "rpt-srna-h", T("mi.num.title")), note(T("mi.num.note")),
        h("div", "rpt-srna-scroll", h("table", "rpt-srna-table", [h("thead", "", head), tb]))]),
      h("section", "rpt-srna-ck", [h("h3", "rpt-srna-h", T("mi.crit.title")), note(T("mi.crit.note")),
        h("div", "rpt-srna-scroll", h("table", "rpt-srna-table", [h("thead", "", chead), crit]))])]);
  }
  function mirnaBlock(mi, r) {
    return h("div", "rpt-srna-figs", [mi.seq ? miFigure(mi, r) : note(T("mi.noseq")), miNumbers(mi)]);
  }

  function openCluster(i) {
    var r = rows[i], d = D.details[String(i)];
    var box = clear(document.getElementById("sec-detail"));
    live = [];
    var back = button(T("ui.detail.back"));
    back.addEventListener("click", function () {
      show("clusters");
      if (listH2) { listH2.setAttribute("tabindex", "-1"); listH2.focus({preventScroll: true}); }
    });
    box.appendChild(back);
    var conf = h("span", "rpt-srna-conf", T("ui.detail.conflabel", {v: confText(r)}));
    conf.setAttribute("title", r.c === "ambiguous" ? T("ui.detail.conf.amb") :
      r.c === "miRNA" ? T("ui.detail.conf.mirna") : META.mode === "ml" ?
      T("ui.detail.conf.ml") : T(META.kingdom === "plant" ? "ui.detail.conf.plant" : "ui.detail.conf",
        {min: META.conf_min, k: META.conf_k}));
    var flag = null;
    if (r.rv) {
      flag = h("span", "rpt-srna-flag", T("ui.review"));
      flag.setAttribute("title", T("ui.review.tip"));
    }
    var h2 = h("h2", "", r.n);
    h2.setAttribute("tabindex", "-1");
    var anEl = null;
    if (d && d.an) {
      var chip = AN_OK.indexOf(r.an) >= 0 ? h("span", "rpt-srna-an " + r.an, r.an) : null;
      anEl = h("p", "rpt-srna-loc", [h("code", "", d.an), chip ? " " : "", chip]);
    }
    add(box, [
      h("div", "rpt-srna-head", [h2, tag(r.c), conf, flag]),
      h("p", "rpt-srna-loc", T("ui.detail.sub", {where: where(r), reads: num(r.r)}) +
        (typeof r.z === "number" ? T("ui.detail.size", {z: r.z}) : "")),
      anEl, META.an && d ? note(T("ui.an.legend")) : null]);
    if (!d) {
      box.appendChild(emptyMsg(T("ui.row.nodetail")));
      show("detail");
      h2.focus({preventScroll: true});
      return;
    }
    add(box, [h("p", "rpt-srna-say", T(d.v, {reason: r.u || ""})), strip(d, r), checklist(d, r)]);
    show("detail");
    h2.focus({preventScroll: true});
    if (d.mi) { box.appendChild(mirnaBlock(d.mi, r)); }
    box.appendChild(figures(d, r));
    live.forEach(function (f) { f.render(); });
    watch();
  }

  buildTabs();
  renderAbout();
  renderOverview();
  renderList();
  show("overview");
  ovLive.forEach(function (f) { f.render(); });
  watch();
  // the indicator slides only after the first placement (and once the fonts settle the widths)
  if (document.fonts && document.fonts.ready) { document.fonts.ready.then(placeIndicator); }
  window.setTimeout(function () { indicator.classList.remove("still"); }, 400);
})();
