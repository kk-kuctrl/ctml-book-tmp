"use strict";

// Minimal canvas-based plotting engine covering what the ported figures
// need: line/scatter/fill-between series, linear or log axes, legends and
// axis labels. Math text (labels) is rendered with KaTeX into HTML overlays
// -- same division of labor as loopshaper (canvas/SVG for data, KaTeX for
// text) -- instead of drawing math glyphs onto the canvas ourselves.

const COLORS = {
  blue: "#2563eb",
  orange: "#f59e0b",
  green: "#16a34a",
  red: "#dc2626",
  purple: "#9333ea",
  black: "#111827",
  gray: "#9ca3af",
};

function resolveColor(c) {
  if (!c) return COLORS.blue;
  if (Array.isArray(c)) {
    const [r, g, b] = c.map((v) => Math.round(v * 255));
    return `rgb(${r},${g},${b})`;
  }
  return COLORS[c] || c;
}

function renderLatex(el, tex) {
  const cleaned = String(tex).replace(/^\$|\$$/g, "");
  try {
    window.katex.render(cleaned, el, { throwOnError: false, displayMode: false });
  } catch (e) {
    el.textContent = cleaned;
  }
}

const MARGIN = { left: 62, right: 18, top: 16, bottom: 46 };

function makeCanvas(container, width, height) {
  container.innerHTML = "";
  container.style.position = "relative";
  const wrap = document.createElement("div");
  wrap.className = "chart-wrap";
  const canvas = document.createElement("canvas");
  const dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr;
  canvas.height = height * dpr;
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  wrap.appendChild(canvas);
  container.appendChild(wrap);
  return { wrap, canvas, ctx, width, height };
}

// Charts are drawn at a fixed logical pixel size (for crisp, predictable
// canvas math), but the grid layout can hand a card less width than that
// (e.g. a 3-per-row layout on a narrower window) -- without this, the
// card's `overflow:auto` just clips the right edge instead of shrinking
// the chart to fit. Scales the whole rendered chart (canvas + legend/axis
// label overlays) down as one unit via a CSS transform, so everything
// stays aligned, and sizes a wrapper div to the resulting on-screen box so
// the card lays out correctly (a transform alone doesn't shrink the
// element's contribution to layout/overflow, only its paint).
//
// The actual measuring/scaling is deferred (see finalizeResponsiveLayout)
// rather than done here immediately: a figure with several cards inserts
// them into the grid one at a time, and CSS Grid's auto-fit recomputes
// column widths as items arrive -- measuring container.clientWidth right
// after inserting just the first card sees a stale, too-wide layout (that
// one card temporarily filling the whole row) instead of the final
// multi-column layout, so its scale decision ends up wrong.
const pendingResponsive = [];
function makeResponsive(container, wrap) {
  pendingResponsive.push({ container, wrap });
}

// Call once after a figure's render() function has finished inserting all
// of its cards, so every container's clientWidth reflects the final grid
// layout (all sibling cards already present) before any scaling decision
// is made.
function finalizeResponsiveLayout() {
  for (const { container, wrap } of pendingResponsive) {
    const naturalW = wrap.offsetWidth;
    const naturalH = wrap.offsetHeight;
    const outer = document.createElement("div");
    outer.style.overflow = "hidden";
    outer.style.margin = "0 auto"; // center the chart when it's narrower than its card
    container.insertBefore(outer, wrap);
    outer.appendChild(wrap);
    const available = container.clientWidth;
    const scale = available > 0 && naturalW > available ? available / naturalW : 1;
    if (scale < 1) {
      wrap.style.transformOrigin = "top left";
      wrap.style.transform = `scale(${scale})`;
    }
    outer.style.width = Math.round(naturalW * scale) + "px";
    outer.style.height = Math.round(naturalH * scale) + "px";
  }
  pendingResponsive.length = 0;
}

function niceTicks(min, max, count = 5) {
  if (min === max) return [min];
  const span = max - min;
  const step0 = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  let step;
  if (norm < 1.5) step = 1 * mag;
  else if (norm < 3) step = 2 * mag;
  else if (norm < 7) step = 5 * mag;
  else step = 10 * mag;
  const start = Math.ceil(min / step) * step;
  const ticks = [];
  for (let v = start; v <= max + step * 1e-9; v += step) ticks.push(Math.round(v / step) * step);
  return ticks;
}

function logTicks(min, max) {
  const lo = Math.floor(Math.log10(min));
  const hi = Math.ceil(Math.log10(max));
  const ticks = [];
  for (let e = lo; e <= hi; e++) ticks.push(Math.pow(10, e));
  return ticks;
}

function fmtNum(v) {
  if (Math.abs(v) < 1e-9) return "0";
  const r = Math.round(v * 1000) / 1000;
  return String(r);
}

function fmtExp(v) {
  const e = Math.round(Math.log10(v));
  return `10^{${e}}`;
}

class Chart2D {
  constructor(container, opts) {
    this.opts = Object.assign(
      {
        width: 460,
        height: 380,
        xlim: [0, 1],
        ylim: [0, 1],
        xscale: "linear",
        yscale: "linear",
        grid: true,
        xlabel: "",
        ylabel: "",
        legendLoc: "upper right",
      },
      opts
    );
    const { wrap, canvas, ctx, width, height } = makeCanvas(container, this.opts.width, this.opts.height);
    this.container = container;
    this.wrap = wrap;
    this.canvas = canvas;
    this.ctx = ctx;
    this.width = width;
    this.height = height;
    this.plotW = width - MARGIN.left - MARGIN.right;
    this.plotH = height - MARGIN.top - MARGIN.bottom;
    this.legendEntries = [];
    this._drawFrame();
  }

  _fx(x) {
    return this.opts.xscale === "log" ? Math.log10(x) : x;
  }
  _fy(y) {
    return this.opts.yscale === "log" ? Math.log10(y) : y;
  }

  toPx(x, y) {
    const [x0, x1] = this.opts.xlim;
    const [y0, y1] = this.opts.ylim;
    const fx0 = this._fx(x0),
      fx1 = this._fx(x1),
      fy0 = this._fy(y0),
      fy1 = this._fy(y1);
    const px = MARGIN.left + ((this._fx(x) - fx0) / (fx1 - fx0)) * this.plotW;
    const py = MARGIN.top + (1 - (this._fy(y) - fy0) / (fy1 - fy0)) * this.plotH;
    return [px, py];
  }

  _drawFrame() {
    const { ctx, width, height, opts } = this;
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    ctx.save();
    ctx.strokeStyle = "#e2e5ea";
    ctx.lineWidth = 1;
    ctx.font = "11px -apple-system, sans-serif";
    ctx.fillStyle = "#6b7280";

    // logTicks rounds outward to whole powers of 10, which can land past
    // xlim/ylim's actual max (e.g. xlim=[0.001, pi] still gets a tick at 10)
    // -- filter those back out so no gridline is drawn outside the axis box.
    const inRange = (t, lim) => t >= lim[0] && t <= lim[1];
    const xTicks = (opts.xscale === "log" ? logTicks(opts.xlim[0], opts.xlim[1]) : niceTicks(opts.xlim[0], opts.xlim[1])).filter((t) =>
      inRange(t, opts.xlim)
    );
    const yTicks = (opts.yscale === "log" ? logTicks(opts.ylim[0], opts.ylim[1]) : niceTicks(opts.ylim[0], opts.ylim[1])).filter((t) =>
      inRange(t, opts.ylim)
    );

    if (opts.grid) {
      for (const t of xTicks) {
        const [px] = this.toPx(t, opts.ylim[0]);
        ctx.beginPath();
        ctx.moveTo(px, MARGIN.top);
        ctx.lineTo(px, MARGIN.top + this.plotH);
        ctx.stroke();
      }
      for (const t of yTicks) {
        const [, py] = this.toPx(opts.xlim[0], t);
        ctx.beginPath();
        ctx.moveTo(MARGIN.left, py);
        ctx.lineTo(MARGIN.left + this.plotW, py);
        ctx.stroke();
      }
    }

    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const t of xTicks) {
      const [px] = this.toPx(t, opts.ylim[0]);
      ctx.fillText(fmtNum(t), px, MARGIN.top + this.plotH + 6);
    }
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const t of yTicks) {
      const [, py] = this.toPx(opts.xlim[0], t);
      ctx.fillText(opts.yscale === "log" ? `1e${Math.round(Math.log10(t))}` : fmtNum(t), MARGIN.left - 8, py);
    }

    ctx.strokeStyle = "#c7cbd1";
    ctx.strokeRect(MARGIN.left, MARGIN.top, this.plotW, this.plotH);
    ctx.restore();
  }

  _clip(fn) {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(MARGIN.left, MARGIN.top, this.plotW, this.plotH);
    ctx.clip();
    fn();
    ctx.restore();
  }

  line(xs, ys, o = {}) {
    const { ctx } = this;
    this._clip(() => {
      ctx.strokeStyle = resolveColor(o.color);
      ctx.globalAlpha = o.alpha !== undefined ? o.alpha : 1;
      ctx.lineWidth = o.lineWidth !== undefined ? o.lineWidth : 1.5;
      if (o.dash) ctx.setLineDash(o.dash);
      ctx.beginPath();
      for (let i = 0; i < xs.length; i++) {
        const [px, py] = this.toPx(xs[i], ys[i]);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
    });
    if (o.label) this.legendEntries.push({ type: "line", color: resolveColor(o.color), label: o.label, dash: o.dash });
    return this;
  }

  scatter(xs, ys, o = {}) {
    const { ctx } = this;
    const size = o.size || 4;
    this._clip(() => {
      ctx.strokeStyle = resolveColor(o.color);
      ctx.fillStyle = resolveColor(o.color);
      for (let i = 0; i < xs.length; i++) {
        const [px, py] = this.toPx(xs[i], ys[i]);
        ctx.beginPath();
        if (o.marker === "star") {
          this._star(px, py, size);
        } else {
          ctx.arc(px, py, size, 0, 2 * Math.PI);
        }
        if (o.filled === false) ctx.stroke();
        else ctx.fill();
      }
    });
    if (o.label)
      this.legendEntries.push({ type: "marker", color: resolveColor(o.color), label: o.label, marker: o.marker, filled: o.filled });
    return this;
  }

  _star(cx, cy, r) {
    const { ctx } = this;
    const spikes = 5;
    const step = Math.PI / spikes;
    let rot = -Math.PI / 2;
    ctx.moveTo(cx + Math.cos(rot) * r, cy + Math.sin(rot) * r);
    for (let i = 0; i < spikes; i++) {
      ctx.lineTo(cx + Math.cos(rot) * r, cy + Math.sin(rot) * r);
      rot += step;
      ctx.lineTo(cx + Math.cos(rot) * r * 0.45, cy + Math.sin(rot) * r * 0.45);
      rot += step;
    }
    ctx.closePath();
  }

  fillBetween(xs, y1, y2, o = {}) {
    const { ctx } = this;
    this._clip(() => {
      ctx.fillStyle = resolveColor(o.color);
      ctx.globalAlpha = o.alpha !== undefined ? o.alpha : 0.25;
      ctx.beginPath();
      for (let i = 0; i < xs.length; i++) {
        const [px, py] = this.toPx(xs[i], y1[i]);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      for (let i = xs.length - 1; i >= 0; i--) {
        const [px, py] = this.toPx(xs[i], y2[i]);
        ctx.lineTo(px, py);
      }
      ctx.closePath();
      ctx.fill();
      ctx.globalAlpha = 1;
    });
    return this;
  }

  finish() {
    const { opts } = this;
    if (opts.xlabel) {
      const el = document.createElement("div");
      el.className = "chart-xlabel";
      renderLatex(el, opts.xlabel);
      this.wrap.appendChild(el);
    }
    if (opts.ylabel) {
      const el = document.createElement("div");
      el.className = "chart-ylabel";
      renderLatex(el, opts.ylabel);
      this.wrap.appendChild(el);
    }
    if (this.legendEntries.length) {
      const box = document.createElement("div");
      box.className = "chart-legend " + (opts.legendLoc === "upper left" ? "loc-ul" : "loc-ur");
      for (const entry of this.legendEntries) {
        const row = document.createElement("div");
        row.className = "chart-legend-row";
        const swatch = document.createElement("span");
        swatch.className = "chart-legend-swatch";
        if (entry.type === "marker") {
          // The base .chart-legend-swatch box is sized for a line swatch
          // (16x3, wide and flat) -- border-radius:50% on that squashes a
          // circle into an ellipse, so give marker swatches their own
          // square footprint instead.
          swatch.style.width = "10px";
          swatch.style.height = "10px";
          swatch.style.background = entry.filled === false ? "transparent" : entry.color;
          swatch.style.border = "2px solid " + entry.color;
          swatch.style.borderRadius = "50%";
        } else {
          swatch.style.background = entry.color;
          if (entry.dash) swatch.style.backgroundImage = `repeating-linear-gradient(90deg, ${entry.color} 0 4px, transparent 4px 7px)`;
        }
        const label = document.createElement("span");
        label.className = "chart-legend-label";
        // Plain-word labels ("True", "white", "初期位置", ...) have no math
        // markup in them at all; running them through KaTeX math mode
        // anyway squashes the spacing and italicizes every letter. Simple
        // math like "u_k" (subscript, no backslash) should still go through
        // KaTeX -- only skip it when there's no LaTeX-ish syntax whatsoever.
        if (/[\\_^$]/.test(entry.label)) renderLatex(label, entry.label);
        else label.textContent = entry.label;
        row.appendChild(swatch);
        row.appendChild(label);
        box.appendChild(row);
      }
      this.wrap.appendChild(box);
    }
    makeResponsive(this.container, this.wrap);
    return this;
  }
}

function createChart(container, opts) {
  return new Chart2D(container, opts);
}

// Very small orthographic 3D line-plot renderer (matplotlib-style
// elev/azim view), sufficient for wireframe-ish 3D curves -- not a general
// 3D engine, just enough for Figure 3.1(a).
// xlim/ylim/zlim are optional -- pass them to fix the axis box (like
// matplotlib's ax.set_xlim/ylim/zlim) instead of auto-fitting to whatever
// this run's data happens to span, which matters for figures whose data is
// randomly generated and would otherwise rescale every run.
function createChart3D(container, opts) {
  const o = Object.assign(
    { width: 480, height: 420, elev: 26, azim: -107, xlabel: "", ylabel: "", zlabel: "", xlim: null, ylim: null, zlim: null },
    opts
  );
  const { wrap, canvas, ctx, width, height } = makeCanvas(container, o.width, o.height);
  // .chart-wrap's CSS padding (58px left / 34px bottom) exists for the 2D
  // Chart2D's external xlabel/ylabel divs. This 3D chart instead draws tick
  // labels on the canvas itself and positions its axis-name divs with
  // absolute left/top computed in the CANVAS's own coordinate space -- but
  // an absolutely-positioned child's (0,0) is the *padding* box's corner,
  // while the canvas (normal flow) starts padding-left/top past that, so
  // that inherited padding both wastes space around the canvas AND shifts
  // every axis label off of where its own math says it should be. Zero it
  // out here rather than in the shared CSS, so Chart2D is unaffected.
  wrap.style.padding = "0";

  const elev = (o.elev * Math.PI) / 180;
  const azim = (o.azim * Math.PI) / 180;

  const series = [];
  function addLine(xs, ys, zs, color, opts) {
    const lineWidth = opts && opts.lineWidth != null ? opts.lineWidth : 1.2;
    series.push({ xs, ys, zs, color: resolveColor(color), lineWidth });
  }

  function finish() {
    // Find each axis's own data range first (raw values, no projection yet)
    // -- matplotlib's default 3D "auto" box aspect scales x/y/z independently
    // to fill the box instead of preserving their relative numeric scale, so
    // e.g. Figure 8.8's k in [-100,100] alongside a state in [-3,3] still
    // both fill the plot rather than the state axis collapsing to a sliver.
    let xMin = Infinity,
      xMax = -Infinity,
      yMin = Infinity,
      yMax = -Infinity,
      zMin = Infinity,
      zMax = -Infinity;
    for (const s of series) {
      for (let i = 0; i < s.xs.length; i++) {
        xMin = Math.min(xMin, s.xs[i]);
        xMax = Math.max(xMax, s.xs[i]);
        yMin = Math.min(yMin, s.ys[i]);
        yMax = Math.max(yMax, s.ys[i]);
        zMin = Math.min(zMin, s.zs[i]);
        zMax = Math.max(zMax, s.zs[i]);
      }
    }
    if (o.xlim) [xMin, xMax] = o.xlim;
    if (o.ylim) [yMin, yMax] = o.ylim;
    if (o.zlim) [zMin, zMax] = o.zlim;
    const xMid = (xMin + xMax) / 2,
      yMid = (yMin + yMax) / 2,
      zMid = (zMin + zMax) / 2;
    const xHalf = (xMax - xMin) / 2 || 1,
      yHalf = (yMax - yMin) / 2 || 1,
      zHalf = (zMax - zMin) / 2 || 1;

    // Normalize each axis to a common +/-1 range before rotating, so the
    // rotation (and the final uniform screen scale) treats all three axes
    // as equally "wide" regardless of their actual numeric units. `xStretch`
    // (etc.) then optionally exaggerates one axis's normalized range beyond
    // +/-1 -- unlike changing the canvas width/height (which only adds or
    // removes margin around an unchanged, aspect-preserved plot), this
    // actually widens that axis's contribution to the rendered geometry.
    const xStretch = o.xStretch || 1,
      yStretch = o.yStretch || 1,
      zStretch = o.zStretch || 1;
    function project(x, y, z) {
      const nx = ((x - xMid) / xHalf) * xStretch,
        ny = ((y - yMid) / yHalf) * yStretch,
        nz = ((z - zMid) / zHalf) * zStretch;
      const x1 = nx * Math.cos(azim) - ny * Math.sin(azim);
      const y1 = nx * Math.sin(azim) + ny * Math.cos(azim);
      const z1 = nz;
      const y2 = y1 * Math.cos(elev) - z1 * Math.sin(elev);
      const z2 = y1 * Math.sin(elev) + z1 * Math.cos(elev);
      return { sx: x1, sy: z2, depth: y2 };
    }

    // Screen-space bounding box, from the box's own corners (not wherever
    // the data happens to land) so the whole axis range is always visible.
    let minSx = Infinity,
      maxSx = -Infinity,
      minSy = Infinity,
      maxSy = -Infinity;
    for (const xv of [xMin, xMax])
      for (const yv of [yMin, yMax])
        for (const zv of [zMin, zMax]) {
          const p = project(xv, yv, zv);
          minSx = Math.min(minSx, p.sx);
          maxSx = Math.max(maxSx, p.sx);
          minSy = Math.min(minSy, p.sy);
          maxSy = Math.max(maxSy, p.sy);
        }

    const projected = series.map((s) => s.xs.map((x, i) => project(x, s.ys[i], s.zs[i])));

    const margin = o.margin != null ? o.margin : 50;
    const plotW = width - 2 * margin;
    const plotH = height - 2 * margin;
    // Default: one shared scale for both screen axes, so shapes that must
    // stay geometrically correct (e.g. Figure 11.1's covariance ellipses)
    // never get stretched. Opt into `stretchToFill` for figures like 8.8
    // where the box is deliberately wide/short and the point is to fill it.
    // `contentScale` shrinks the rendered content within its plot area
    // (independent of canvas size / margin), leaving extra blank space
    // around it -- e.g. a canvas at 90% size with content at 80% fill.
    const contentScale = o.contentScale != null ? o.contentScale : 1;
    const baseScaleX = o.stretchToFill ? plotW / (maxSx - minSx || 1) : Math.min(plotW / (maxSx - minSx || 1), plotH / (maxSy - minSy || 1));
    const baseScaleY = o.stretchToFill ? plotH / (maxSy - minSy || 1) : baseScaleX;
    const scaleX = baseScaleX * contentScale;
    const scaleY = baseScaleY * contentScale;
    const cx = width / 2,
      cy = height / 2;
    const midSx = (minSx + maxSx) / 2,
      midSy = (minSy + maxSy) / 2;
    const toPx = (p) => [cx + (p.sx - midSx) * scaleX, cy - (p.sy - midSy) * scaleY];

    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    // Grid on the three "floor/wall" panes (x=xMin, y=yMin, z=zMin), the
    // way matplotlib's 3D axes shade grid lines onto the background panes.
    const xTicks = niceTicks(xMin, xMax);
    const yTicks = niceTicks(yMin, yMax);
    const zTicks = niceTicks(zMin, zMax);
    ctx.save();
    ctx.strokeStyle = "#d6dae0";
    ctx.lineWidth = 1;
    function gridLine(p0, p1) {
      const a = toPx(project(p0[0], p0[1], p0[2]));
      const b = toPx(project(p1[0], p1[1], p1[2]));
      ctx.beginPath();
      ctx.moveTo(a[0], a[1]);
      ctx.lineTo(b[0], b[1]);
      ctx.stroke();
    }
    // Floor (z = zMin): lines along y for each x-tick, along x for each y-tick.
    for (const xt of xTicks) gridLine([xt, yMin, zMin], [xt, yMax, zMin]);
    for (const yt of yTicks) gridLine([xMin, yt, zMin], [xMax, yt, zMin]);
    // Wall (x = xMin): lines along z for each y-tick, along y for each z-tick.
    for (const yt of yTicks) gridLine([xMin, yt, zMin], [xMin, yt, zMax]);
    for (const zt of zTicks) gridLine([xMin, yMin, zt], [xMin, yMax, zt]);
    // Wall (y = yMin): lines along z for each x-tick, along x for each z-tick.
    for (const xt of xTicks) gridLine([xt, yMin, zMin], [xt, yMin, zMax]);
    for (const zt of zTicks) gridLine([xMin, yMin, zt], [xMax, yMin, zt]);
    // Box edges, slightly darker, for a clean frame.
    ctx.strokeStyle = "#b7bcc4";
    const corners = {
      a: [xMin, yMin, zMin],
      b: [xMax, yMin, zMin],
      c: [xMin, yMax, zMin],
      d: [xMin, yMin, zMax],
    };
    gridLine(corners.a, corners.b);
    gridLine(corners.a, corners.c);
    gridLine(corners.a, corners.d);
    ctx.restore();

    // Tick number labels along the three box edges meeting at (xMin,yMin,zMin).
    ctx.save();
    ctx.font = "10.5px -apple-system, sans-serif";
    ctx.fillStyle = "#6b7280";
    function tickLabel(p, text, align, baseline, dx, dy) {
      const [px, py] = toPx(project(p[0], p[1], p[2]));
      ctx.textAlign = align;
      ctx.textBaseline = baseline;
      ctx.fillText(text, px + dx, py + dy);
    }
    for (const xt of xTicks) tickLabel([xt, yMin, zMin], fmtNum(xt), "center", "top", 0, 4);
    for (const yt of yTicks) tickLabel([xMax, yt, zMin], fmtNum(yt), "left", "middle", 4, 0);
    for (const zt of zTicks) tickLabel([xMin, yMin, zt], fmtNum(zt), "right", "middle", -4, 0);
    ctx.restore();

    for (let s = 0; s < series.length; s++) {
      const pts = projected[s];
      ctx.strokeStyle = series[s].color;
      ctx.lineWidth = series[s].lineWidth;
      ctx.beginPath();
      pts.forEach((p, i) => {
        const [px, py] = toPx(p);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      });
      ctx.stroke();
    }

    // Axis-name labels follow the actual projected geometry (same as the
    // tick numbers above) rather than sitting at fixed canvas corners --
    // otherwise they drift away from their own axis whenever elev/azim
    // changes which edge ends up where on screen.
    function axisLabel(midPoint, text, dx, dy) {
      const [px, py] = toPx(project(midPoint[0], midPoint[1], midPoint[2]));
      const el = document.createElement("div");
      el.className = "chart-3d-axis-label";
      el.style.left = px + dx + "px";
      el.style.top = py + dy + "px";
      renderLatex(el, text);
      wrap.appendChild(el);
    }
    if (o.xlabel) axisLabel([(xMin + xMax) / 2, yMin, zMin], o.xlabel, 0, 26);
    if (o.ylabel) axisLabel([xMax, (yMin + yMax) / 2, zMin], o.ylabel, 30, 0);
    if (o.zlabel) axisLabel([xMin, yMin, (zMin + zMax) / 2], o.zlabel, -30, 0);
    makeResponsive(container, wrap);
  }

  return { addLine, finish };
}

// value: 2D array [row][col] in [0,1]-normalizable range -> heatmap + colorbar,
// standing in for matplotlib's contourf (a filled heatmap reads the same way
// for this kind of smooth scalar field, without needing a marching-squares
// contour implementation).
function turbo(t) {
  t = Math.min(1, Math.max(0, t));
  const r = Math.round(255 * Math.min(1, Math.max(0, 1.5 - Math.abs(4 * t - 3))));
  const g = Math.round(255 * Math.min(1, Math.max(0, 1.5 - Math.abs(4 * t - 2))));
  const b = Math.round(255 * Math.min(1, Math.max(0, 1.5 - Math.abs(4 * t - 1))));
  return [r, g, b];
}

function createHeatmap(container, grid, opts) {
  const o = Object.assign({ width: 420, height: 380 }, opts);
  container.innerHTML = "";
  container.style.position = "relative";
  const wrap = document.createElement("div");
  wrap.className = "chart-wrap heatmap-wrap";

  const rows = grid.length,
    cols = grid[0].length;
  let vmin = Infinity,
    vmax = -Infinity;
  for (const row of grid) for (const v of row) {
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
  }

  const plotSize = o.width - MARGIN.left - MARGIN.right;
  // grid[i][j] is (x=i, y=j) -- so the offscreen buffer's WIDTH spans i
  // (rows.length) and its HEIGHT spans j (cols.length), with pixel
  // (col=i, row=j). Getting i/j swapped here previously put j on the x-axis
  // instead of i, mismatching the simulation's own state[:,0]=x, state[:,1]=y
  // convention used elsewhere in the same figure.
  const off = document.createElement("canvas");
  off.width = rows;
  off.height = cols;
  const octx = off.getContext("2d");
  const img = octx.createImageData(rows, cols);
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      const t = (grid[i][j] - vmin) / (vmax - vmin || 1);
      const [r, g, b] = turbo(t);
      const idx = (j * rows + i) * 4;
      img.data[idx] = r;
      img.data[idx + 1] = g;
      img.data[idx + 2] = b;
      img.data[idx + 3] = 255;
    }
  }
  octx.putImageData(img, 0, 0);

  const canvas = document.createElement("canvas");
  const dpr = window.devicePixelRatio || 1;
  canvas.width = o.width * dpr;
  canvas.height = o.height * dpr;
  canvas.style.width = o.width + "px";
  canvas.style.height = o.height + "px";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.imageSmoothingEnabled = true;
  // Flip vertically so j (y) increases upward like matplotlib, matching the
  // i=x, j=y orientation the offscreen buffer was built with above.
  ctx.save();
  ctx.translate(MARGIN.left, MARGIN.top + plotSize);
  ctx.scale(1, -1);
  ctx.drawImage(off, 0, 0, rows, cols, 0, 0, plotSize, plotSize);
  ctx.restore();
  ctx.strokeStyle = "#c7cbd1";
  ctx.strokeRect(MARGIN.left, MARGIN.top, plotSize, plotSize);

  // Colorbar
  const cbW = 14;
  const cbX = MARGIN.left + plotSize + 14;
  const cbGrad = ctx.createLinearGradient(0, MARGIN.top, 0, MARGIN.top + plotSize);
  for (let s = 0; s <= 1; s += 0.05) {
    const [r, g, b] = turbo(1 - s);
    cbGrad.addColorStop(s, `rgb(${r},${g},${b})`);
  }
  ctx.fillStyle = cbGrad;
  ctx.fillRect(cbX, MARGIN.top, cbW, plotSize);
  ctx.strokeRect(cbX, MARGIN.top, cbW, plotSize);
  ctx.fillStyle = "#6b7280";
  ctx.font = "11px -apple-system, sans-serif";
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  for (let s = 0; s <= 1.001; s += 0.25) {
    const val = vmin + s * (vmax - vmin);
    const py = MARGIN.top + plotSize - s * plotSize;
    ctx.fillText(fmtNum(val), cbX + cbW + 4, py);
  }

  wrap.appendChild(canvas);
  container.appendChild(wrap);
  makeResponsive(container, wrap);
  return { canvas, ctx, plotSize, margin: MARGIN };
}

function makeCard(outputGrid, subtitle) {
  const card = document.createElement("div");
  card.className = "output-card";
  if (subtitle) {
    const h = document.createElement("div");
    h.className = "output-card-title";
    h.textContent = subtitle;
    card.appendChild(h);
  }
  const body = document.createElement("div");
  body.className = "output-card-body";
  card.appendChild(body);
  outputGrid.appendChild(card);
  return body;
}

window.plotlib = {
  createChart,
  createChart3D,
  createHeatmap,
  COLORS,
  resolveColor,
  renderLatex,
  makeCard,
  finalizeResponsiveLayout,
};
