// Runs report_assets/app.js under JavaScriptCore (osascript -l JavaScript) against a stub DOM.
// The Python side prepends `var CFG = {app, data, probe, width}`; the last expression is a JSON string.
// It implements only the DOM surface app.js touches: a missing method throws, on purpose.
(function () {
  var WIDTH = CFG.width || 720;
  var stats = {downloads: [], focus: 0, timers: [], raf: []};

  function ClassList(e) { this.e = e; }
  ClassList.prototype = {
    _get: function () { return (this.e.ns ? this.e.attrs["class"] || "" : this.e.className).split(/\s+/).filter(Boolean); },
    _set: function (a) { if (this.e.ns) { this.e.attrs["class"] = a.join(" "); } else { this.e.className = a.join(" "); } },
    contains: function (c) { return this._get().indexOf(c) >= 0; },
    add: function (c) { var a = this._get(); if (a.indexOf(c) < 0) { a.push(c); } this._set(a); },
    remove: function (c) { this._set(this._get().filter(function (x) { return x !== c; })); },
    toggle: function (c, force) {
      var on = force === undefined ? !this.contains(c) : !!force;
      if (on) { this.add(c); } else { this.remove(c); }
      return on;
    }
  };

  function El(tag, ns) {
    this.tagName = tag; this.ns = ns || null; this.children = []; this.parent = null;
    this.attrs = {}; this.listeners = {}; this.className = ""; this.value = ""; this._text = null;
    this.style = {}; this.classList = new ClassList(this);
  }
  El.prototype = {
    get firstChild() { return this.children[0] || null; },
    get parentNode() { return this.parent; },
    get textContent() {
      if (this._text !== null) { return this._text; }
      return this.children.map(function (c) { return c.textContent; }).join("");
    },
    set textContent(v) { this.children = []; this._text = String(v); },
    // project rule: no markup from strings
    get innerHTML() { throw new Error("innerHTML is not allowed"); },
    set innerHTML(v) { throw new Error("innerHTML is not allowed"); },
    get clientWidth() { return WIDTH; },
    get offsetWidth() { return this.tagName === "button" ? 96 : 200; },
    get offsetHeight() { return 60; },
    get offsetLeft() { return this.parent ? this.parent.children.indexOf(this) * 96 : 0; },
    getBoundingClientRect: function () { return {left: 10, top: 20, right: 34, bottom: 44, width: 24, height: 24}; },
    appendChild: function (c) {
      if (!c || !(c instanceof El || c instanceof Txt)) { throw new TypeError("appendChild: not a node"); }
      if (c.parent) { c.parent.removeChild(c); }
      this._text = null;
      c.parent = this; this.children.push(c); return c;
    },
    insertBefore: function (c, ref) {
      if (!ref) { return this.appendChild(c); }
      var i = this.children.indexOf(ref);
      if (i < 0) { throw new Error("insertBefore: not a child"); }
      if (c.parent) { c.parent.removeChild(c); i = this.children.indexOf(ref); }
      c.parent = this; this.children.splice(i, 0, c); return c;
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
    removeAttribute: function (k) { delete this.attrs[k]; },
    cloneNode: function (deep) {
      var c = new El(this.tagName, this.ns), self = this;
      Object.keys(this.attrs).forEach(function (k) { c.attrs[k] = self.attrs[k]; });
      c.className = this.className; c._text = this._text;
      if (deep) {
        this.children.forEach(function (k) {
          c.appendChild(k instanceof Txt ? new Txt(k.data) : k.cloneNode(true));
        });
      }
      return c;
    },
    getAttribute: function (k) { return k in this.attrs ? this.attrs[k] : null; },
    hasAttribute: function (k) { return k in this.attrs; },
    addEventListener: function (t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); },
    removeEventListener: function (t, f) {
      this.listeners[t] = (this.listeners[t] || []).filter(function (x) { return x !== f; });
    },
    focus: function () {
      var prev = document.activeElement;
      if (prev === this) { return; }
      if (prev && prev.listeners) { (prev.listeners.blur || []).forEach(function (f) { f.call(prev, {type: "blur"}); }); }
      document.activeElement = this;
      stats.focus += 1;
      var self = this;
      (this.listeners.focus || []).forEach(function (f) { f.call(self, {type: "focus", target: self}); });
    },
    click: function () {
      var self = this;
      (this.listeners.click || []).forEach(function (f) { f.call(self, {type: "click", target: self}); });
      if (this.tagName === "a" && "download" in this.attrs) {
        stats.downloads.push({name: this.attrs.download, href: this.attrs.href || this.href});
      }
    }
  };
  function Txt(s) { this.data = s; this.parent = null; }
  Txt.prototype = { get textContent() { return this.data; } };

  var ids = {};
  function root(id, tag) { var e = new El(tag || "div"); ids[id] = e; return e; }
  ["tabs", "navmeta", "sec-about", "sec-overview", "sec-clusters", "sec-detail", "tip"].forEach(
    function (id) { root(id); });
  root("report-data", "script")._text = CFG.data;
  var body = new El("body");

  var document = {
    body: body,
    activeElement: null,
    listeners: {},
    createElement: function (t) { return new El(t); },
    createElementNS: function (ns, t) { return new El(t, ns); },
    createTextNode: function (s) { return new Txt(String(s)); },
    getElementById: function (id) { return ids[id] || null; },
    addEventListener: function (t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); }
  };

  var observers = [];
  function ResizeObserver(cb) { this.cb = cb; this.targets = []; observers.push(this); }
  ResizeObserver.prototype = {
    observe: function (e) { this.targets.push(e); },
    unobserve: function (e) { this.targets = this.targets.filter(function (x) { return x !== e; }); },
    disconnect: function () { this.targets = []; }
  };

  function esc(s, attr) {
    s = String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    return attr ? s.replace(/"/g, "&quot;") : s;
  }
  function ser(e, top) {
    if (e instanceof Txt) { return esc(e.data); }
    var a = "";
    if (top && e.ns && !("xmlns" in e.attrs)) { a += ' xmlns="' + e.ns + '"'; }
    Object.keys(e.attrs).forEach(function (k) { a += " " + k + '="' + esc(e.attrs[k], true) + '"'; });
    if (!e.ns && e.className) { a += ' class="' + esc(e.className, true) + '"'; }
    var inner = e._text !== null ? esc(e._text) : e.children.map(function (c) { return ser(c); }).join("");
    return "<" + e.tagName + a + ">" + inner + "</" + e.tagName + ">";
  }
  function XMLSerializer() {}
  XMLSerializer.prototype.serializeToString = function (e) { return ser(e, true); };
  function Blob(parts, opts) { this.parts = parts; this.type = opts && opts.type; }

  var scrolls = 0;
  var window = {
    innerWidth: 1280, innerHeight: 800,
    listeners: {},
    scrollTo: function () { scrolls += 1; },
    addEventListener: function (t, f) { (this.listeners[t] = this.listeners[t] || []).push(f); },
    setTimeout: function (f) { stats.timers.push(f); return stats.timers.length; },
    clearTimeout: function () {},
    // frames run only when the probe says so (H.frames), so a throttled redraw is observable
    requestAnimationFrame: function (f) { stats.rafId = (stats.rafId || 0) + 1; stats.raf.push({id: stats.rafId, f: f}); return stats.rafId; },
    cancelAnimationFrame: function (id) { stats.raf = stats.raf.filter(function (q) { return q.id !== id; }); },
    ResizeObserver: ResizeObserver,
    XMLSerializer: XMLSerializer,
    Blob: Blob,
    URL: {
      createObjectURL: function (b) { stats.lastBlob = b.parts.join(""); return "blob:stub"; },
      revokeObjectURL: function () {}
    }
  };

  // ---- helpers for the probe ----
  function walk(e, f) {
    f(e);
    (e.children || []).forEach(function (c) { if (c instanceof El) { walk(c, f); } });
  }
  function cls(e) { return e.ns ? e.attrs["class"] || "" : e.className || ""; }
  var H = {
    id: function (id) { return ids[id]; },
    all: function (scope, tag) {
      var out = [];
      walk(scope, function (e) { if (!tag || e.tagName === tag) { out.push(e); } });
      return out;
    },
    // every element the app can reach from the page roots (the data block is the page's own)
    doc: function (tag) {
      var out = [];
      Object.keys(ids).filter(function (k) { return k !== "report-data"; }).concat(["body"]).forEach(function (k) {
        out = out.concat(H.all(k === "body" ? body : ids[k], tag));
      });
      return out;
    },
    count: function (scope, tag) { return H.all(scope, tag).length; },
    counts: function (scope) {
      var o = {};
      ["svg", "rect", "path", "circle", "line", "polyline", "tr", "button", "table", "text"].forEach(
        function (t) { o[t] = H.count(scope, t); });
      return o;
    },
    byClass: function (scope, c, tag) {
      return H.all(scope, tag).filter(function (e) { return (" " + cls(e) + " ").indexOf(" " + c + " ") >= 0; });
    },
    cls: cls,
    text: function (e) { return e.textContent; },
    buttons: function (scope, label) {
      return H.all(scope, "button").filter(function (b) { return b.textContent === label; });
    },
    fire: function (e, type, ev) {
      var fs = e.listeners[type] || [];
      if (!fs.length) { throw new Error("no " + type + " handler on <" + e.tagName + ">"); }
      ev = ev || {};
      ev.type = type; ev.target = ev.target || e;
      ev.preventDefault = ev.preventDefault || function () { ev.prevented = true; };
      fs.slice().forEach(function (f) { f.call(e, ev); });
      return ev;
    },
    has: function (e, type) { return (e.listeners[type] || []).length > 0; },
    click: function (e) { H.fire(e, "click"); },
    key: function (e, key) { return H.fire(e, "keydown", {key: key}); },
    setValue: function (e, v, type) { e.value = String(v); H.fire(e, type || "change"); },
    active: function () {
      return ["about", "overview", "clusters", "detail"].filter(function (s) {
        return / active$/.test(ids["sec-" + s].className);
      });
    },
    resize: function (w) {
      WIDTH = w;
      observers.forEach(function (o) {
        if (o.targets.length) { o.cb(o.targets.map(function (t) { return {target: t}; })); }
      });
    },
    observed: function () {
      var n = 0;
      observers.forEach(function (o) { n += o.targets.length; });
      return n;
    },
    serialize: function (e) { return ser(e, true); },
    stats: function () { return stats; },
    // run the queued animation frames; returns how many callbacks ran
    frames: function () {
      var n = 0, q;
      while (stats.raf.length && n < 100) {
        q = stats.raf; stats.raf = [];
        q.forEach(function (it) { it.f(0); n += 1; });
      }
      return n;
    },
    pendingFrames: function () { return stats.raf.length; },
    scrolls: function () { return scrolls; },
    activeEl: function () { return document.activeElement; }
  };

  var out;
  try {
    new Function("document", "window", CFG.app)(document, window);
    var probe = new Function("H", "document", "window", CFG.probe);
    out = {ok: true, result: probe(H, document, window)};
  } catch (err) {
    out = {ok: false, error: String(err && err.name ? err.name + ": " + err.message : err)};
  }
  return JSON.stringify(out);
})()
