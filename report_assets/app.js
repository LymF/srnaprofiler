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
  }

  // ---- cluster list ----
  var view = {cls: "all", sort: "conf", page: 0, floor: META.floor};
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
    add(ctl, [sel, h("label", "rpt-srna-lbl", [T("ui.floor.label"), inp])]);
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
        h("td", "num", num(r.r)),
        h("td", "rpt-srna-sub", last)]);
      tb.appendChild(tr);
    });
    var head = h("tr", "", [T("ui.th.cluster"), T("ui.th.class"), T("ui.th.conf"), T("ui.th.evidence"), T("ui.th.reads"), ""].map(function (t, j) {
      return h("th", j === 4 ? "num" : "", t);
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
    var plot = h("div", "rpt-srna-fig__plot");
    plot.style.height = o.height + "px";
    var tbl = h("div", "rpt-srna-fig__table");
    tbl.setAttribute("hidden", "");
    tbl.setAttribute("tabindex", "0");
    tbl.setAttribute("role", "region");
    tbl.setAttribute("aria-label", o.letter + " " + o.title + " " + T("ui.fig.table"));
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
      svg = svgRoot(cw, o.height, o.letter + " " + o.title);
      drawCtx = {marks: {}, focus: {}, sets: 0};
      o.draw(svg, cw, o.height);
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
          add(tbl, [note(t.exact ? T("ui.fig.exact") : T("ui.fig.precision")), dataTable(t.head, t.rows)]);
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
    var title = h("h3", "rpt-srna-fig__title", [h("span", "rpt-srna-fig__letter", o.letter),
      h("span", "", o.title), o.n ? h("span", "rpt-srna-fig__n", o.n) : null]);
    if (o.tip) { title.setAttribute("title", o.tip); }
    var el = h("section", "rpt-srna-fig" + (o.wide ? " rpt-srna-fig--wide" : ""), [
      h("div", "rpt-srna-fig__head", [title, tools]), leg,
      o.body || plot, o.body ? null : tbl]);
    el.setAttribute("data-fig", o.letter);
    if (!o.body) { live.push({plot: plot, render: render}); }
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
      ul.appendChild(h("li", "", [key, it.label]));
    });
    return ul;
  }
  function dataTable(head, body) {
    var tr = h("tr", "", head.map(function (t, j) { return h("th", j ? "num" : "", t); }));
    var tb = h("tbody");
    body.forEach(function (r) {
      tb.appendChild(h("tr", "", r.map(function (c, j) { return h("td", j ? "num" : "", c); })));
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
    var rlo = lr.length ? Math.floor(Math.min.apply(null, lr) * 4) / 4 : 0;
    var rhi = lr.length ? Math.ceil(Math.max.apply(null, lr) * 4 + 1e-9) / 4 : 1;
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
  function fileName(r, letter) {
    var base = String(r.n).replace(/[^\w.-]+/g, "_").replace(/^_+|_+$/g, "").slice(0, 80);
    return (base || T("ui.brand")) + "_" + letter;
  }

  function figures(d, r) {
    var sig = d.sig;
    if (!sig) { return emptyMsg(T("ui.fig.na")); }
    if (sig.err) { return emptyMsg(T("ui.fig.err")); }
    if (sig.none) { return emptyMsg(T("ui.fig.noreads")); }
    var dd = decode(sig);
    return h("div", "rpt-srna-figs", [
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
      T("ui.detail.conf.ml") : T("ui.detail.conf", {min: META.conf_min, k: META.conf_k}));
    var flag = null;
    if (r.rv) {
      flag = h("span", "rpt-srna-flag", T("ui.review"));
      flag.setAttribute("title", T("ui.review.tip"));
    }
    var h2 = h("h2", "", r.n);
    h2.setAttribute("tabindex", "-1");
    add(box, [
      h("div", "rpt-srna-head", [h2, tag(r.c), conf, flag]),
      h("p", "rpt-srna-loc", T("ui.detail.sub", {where: where(r), reads: num(r.r)}) +
        (typeof r.z === "number" ? T("ui.detail.size", {z: r.z}) : ""))]);
    if (!d) {
      box.appendChild(emptyMsg(T("ui.row.nodetail")));
      show("detail");
      h2.focus({preventScroll: true});
      return;
    }
    add(box, [h("p", "rpt-srna-say", T(d.v, {reason: r.u || ""})), strip(d, r), checklist(d, r)]);
    show("detail");
    h2.focus({preventScroll: true});
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
