// Runs report_assets/app.js under JavaScriptCore (osascript -l JavaScript) against a stub DOM.
// The Python side prepends `var CFG = {app, data, probe}`; the last expression is a JSON string.
// It implements only the DOM surface app.js touches: a missing method throws, on purpose.
(function () {
  function El(tag, ns) {
    this.tagName = tag; this.ns = ns || null; this.children = []; this.parent = null;
    this.attrs = {}; this.listeners = {}; this.className = ""; this.value = ""; this._text = null;
  }
  El.prototype = {
    get firstChild() { return this.children[0] || null; },
    get textContent() {
      if (this._text !== null) { return this._text; }
      return this.children.map(function (c) { return c.textContent; }).join("");
    },
    set textContent(v) { this.children = []; this._text = String(v); },
    appendChild: function (c) {
      if (!c || !(c instanceof El || c instanceof Txt)) { throw new TypeError("appendChild: not a node"); }
      c.parent = this; this.children.push(c); return c;
    },
    removeChild: function (c) {
      var i = this.children.indexOf(c);
      if (i < 0) { throw new Error("removeChild: not a child"); }
      this.children.splice(i, 1); c.parent = null; return c;
    },
    setAttribute: function (k, v) {
      this.attrs[k] = String(v);
      if (k === "value") { this.value = String(v); }
    },
    getAttribute: function (k) { return k in this.attrs ? this.attrs[k] : null; },
    addEventListener: function (t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); }
  };
  function Txt(s) { this.data = s; this.parent = null; }
  Txt.prototype = { get textContent() { return this.data; } };

  var ids = {};
  function root(id, tag) { var e = new El(tag || "div"); ids[id] = e; return e; }
  ["tabs", "navmeta", "sec-overview", "sec-clusters", "sec-detail", "sec-glossary",
    "sec-references"].forEach(function (id) { root(id); });
  root("report-data", "script")._text = CFG.data;

  var scrolls = 0;
  var document = {
    createElement: function (t) { return new El(t); },
    createElementNS: function (ns, t) { return new El(t, ns); },
    createTextNode: function (s) { return new Txt(String(s)); },
    getElementById: function (id) { return ids[id] || null; }
  };
  var window = { scrollTo: function () { scrolls += 1; } };

  // ---- helpers for the probe ----
  function walk(e, f) {
    f(e);
    (e.children || []).forEach(function (c) { if (c instanceof El) { walk(c, f); } });
  }
  function cls(e) { return e.className || e.attrs["class"] || ""; }
  var H = {
    id: function (id) { return ids[id]; },
    all: function (scope, tag) {
      var out = [];
      walk(scope, function (e) { if (!tag || e.tagName === tag) { out.push(e); } });
      return out;
    },
    count: function (scope, tag) { return H.all(scope, tag).length; },
    counts: function (scope) {
      var o = {};
      ["svg", "rect", "path", "circle", "line", "polyline", "tr", "button", "table"].forEach(
        function (t) { o[t] = H.count(scope, t); });
      return o;
    },
    byClass: function (scope, c, tag) {
      return H.all(scope, tag).filter(function (e) { return (" " + cls(e) + " ").indexOf(" " + c + " ") >= 0; });
    },
    text: function (e) { return e.textContent; },
    buttons: function (scope, label) {
      return H.all(scope, "button").filter(function (b) { return b.textContent === label; });
    },
    fire: function (e, type, ev) {
      var fs = e.listeners[type] || [];
      if (!fs.length) { throw new Error("no " + type + " handler on <" + e.tagName + ">"); }
      fs.slice().forEach(function (f) { f.call(e, ev || {type: type}); });
    },
    click: function (e) { H.fire(e, "click"); },
    setValue: function (e, v, type) { e.value = String(v); H.fire(e, type || "change"); },
    active: function () {
      return ["overview", "clusters", "detail", "glossary", "references"].filter(function (s) {
        return / active$/.test(ids["sec-" + s].className);
      });
    },
    scrolls: function () { return scrolls; }
  };

  var out;
  try {
    new Function("document", "window", CFG.app)(document, window);
    var probe = new Function("H", "document", CFG.probe);
    out = {ok: true, result: probe(H, document)};
  } catch (err) {
    out = {ok: false, error: String(err && err.name ? err.name + ": " + err.message : err)};
  }
  return JSON.stringify(out);
})()
