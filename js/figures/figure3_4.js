"use strict";

// Port of code_distribute/python/Figure3_4.py
// F_c = 1/(s+a) has the trivial minimal realization Ac=[-a], Bc=[1],
// Cc=[1], Dc=[0] (dy/dt = -a*y + u), so no general TF->SS conversion is
// needed -- just the Tustin discretization + discrete-Lyapunov
// normalization + frequency response, all via lib/linalg.js.
// (The book's original fixes a=0.3; here it's an adjustable parameter.)
window.figureLib = window.figureLib || {};
window.figureLib.figure3_4 = function (outputGrid, params) {
  const Ts = params.ts;
  const kBar = Math.max(2, Math.round(params.k_bar));
  const a = params.a;

  const Ac = [[-a]],
    Bc = [[1]],
    Cc = [[1]],
    Dc = [[0]];
  const { Ad, Bd, Cd } = window.linalg.c2dTustin(Ac, Bc, Cc, Dc, Ts);

  // Normalize so the stationary output variance is 1: solve the discrete
  // Lyapunov equation P = Ad P Ad^T + Bd Bd^T, then scale Bd by 1/sqrt(Cd P Cd^T).
  const Q = window.linalg.matMul(Bd, window.linalg.transpose(Bd));
  const P = window.linalg.discreteLyapunov(Ad, Q);
  const varAmp = window.linalg.matMul(Cd, window.linalg.matMul(P, window.linalg.transpose(Cd)))[0][0];
  const Bnorm = window.linalg.scale(Bd, 1 / Math.sqrt(varAmp));

  // n=1 system -- work with plain scalars from here on.
  const Aw = Ad[0][0];
  const Bw = Bnorm[0][0];
  const Cw = Cd[0][0];

  // Turns an array of N step-heights into the doubled x/y arrays that draw
  // as a proper staircase (matplotlib's plt.stairs) when passed to chart.line.
  function stairsXY(values) {
    const xs = [0];
    const ys = [values[0]];
    for (let i = 0; i < values.length; i++) {
      xs.push(i + 1);
      ys.push(values[i]);
      if (i + 1 < values.length) {
        xs.push(i + 1);
        ys.push(values[i + 1]);
      }
    }
    return { xs, ys };
  }

  // ---------------------------------------------------------------
  // Figure 3.4(a): frequency response magnitude of the (normalized)
  // discrete noise-shaping filter, H(e^{jw}) = Cw / (e^{jw} - Aw).
  // ---------------------------------------------------------------
  {
    const omegas = [];
    for (let w = 0.001; w < Math.PI; w += 0.001) omegas.push(w);
    // H(e^{jw}) = Cw * (e^{jw} - Aw)^{-1} * Bw  (D=0), so the magnitude
    // needs the Bw factor too -- dropping it here previously left the
    // curve scaled up by 1/|Bw| (Bw < 1), pushing it off the top of the
    // [0,9] y-range for most of the low-frequency end.
    const mags = omegas.map((w) => {
      const re = Math.cos(w) - Aw;
      const im = Math.sin(w);
      return (Math.abs(Cw) * Math.abs(Bw)) / Math.sqrt(re * re + im * im);
    });

    // DC gain |H(1)| = |Cw*Bw| / |1-Aw| -- the actual peak of the fitted
    // filter, which does move with a (a narrower-band filter concentrates
    // the same unit stationary variance into a taller peak). The red
    // "prior information" line is a fixed design target from the book,
    // independent of a, so only the axis range adapts (to whichever is
    // taller: the fixed reference or the current curve) so neither one
    // gets clipped.
    const dcGain = (Math.abs(Cw) * Math.abs(Bw)) / Math.abs(1 - Aw);
    const priorHeight = 8;

    const body = window.plotlib.makeCard(outputGrid, "Figure 3.4(a) — 周波数重みの大きさ");
    const chart = window.plotlib.createChart(body, {
      xscale: "log",
      xlim: [0.001, Math.PI],
      ylim: [0, Math.max(dcGain, priorHeight) * 1.15],
      xlabel: "\\varpi",
    });
    chart.line(omegas, mags, { color: "blue", lineWidth: 1, label: "frequency weight" });
    chart.line([0.001, 0.01, 0.6, Math.PI], [priorHeight, priorHeight, 0, 0], {
      color: "red",
      lineWidth: 2,
      dash: [6, 4],
      label: "prior information",
    });
    chart.finish();
  }

  // ---------------------------------------------------------------
  // Figure 3.4(b): white noise v_k filtered through (Aw,Bw,Cw) to
  // produce colored noise y_k.
  // ---------------------------------------------------------------
  {
    const v = window.rnd.randnVec(kBar + 1);
    const xBar = new Array(kBar + 1).fill(0);
    const y = new Array(kBar + 1).fill(0);
    for (let k = 0; k < kBar; k++) {
      xBar[k + 1] = Aw * xBar[k] + Bw * v[k];
      y[k + 1] = Cw * xBar[k + 1];
    }

    const nShow = Math.min(201, kBar + 1);
    const vShow = v.slice(0, nShow);
    const yShow = y.slice(0, nShow);

    const body = window.plotlib.makeCard(outputGrid, "Figure 3.4(b) — 白色雑音と色付き雑音");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, Math.min(200, kBar)],
      ylim: [-2, 2],
      xlabel: "$k$",
      legendLoc: "upper right",
    });
    const whiteStairs = stairsXY(vShow);
    const coloredStairs = stairsXY(yShow);
    chart.line(whiteStairs.xs, whiteStairs.ys, { color: "blue", lineWidth: 1, label: "white" });
    chart.line(coloredStairs.xs, coloredStairs.ys, { color: "orange", lineWidth: 1.5, label: "colored" });
    chart.finish();
  }
};
