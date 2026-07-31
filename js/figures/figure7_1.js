"use strict";

// Port of code_distribute/python/Figure7_1.py
window.figureLib = window.figureLib || {};
window.figureLib.figure7_1 = function (outputGrid, params) {
  const nSampleA = Math.max(1, Math.round(params.n_sample_a));
  const nSampleB = Math.max(1, Math.round(params.n_sample_b));
  const sBar = Math.max(1, Math.round(params.s_bar));

  const nX = 100; // n_x - number of x-grid for plot
  const xP = Array.from({ length: nX }, (_, i) => i / (nX - 1)); // x grid points for plot
  const nF = 10; // number of features

  // Feature mapping phi(x) = [1, x, x^2, ..., x^9]
  function phi(x) {
    const out = new Array(nF);
    let p = 1;
    for (let k = 0; k < nF; k++) {
      out[k] = p;
      p *= x;
    }
    return out;
  }

  // dot(paramVec, phi(x)) evaluated over an array of x values
  function evalCurve(paramVec, xs) {
    return xs.map((x) => {
      const ph = phi(x);
      let s = 0;
      for (let k = 0; k < nF; k++) s += paramVec[k] * ph[k];
      return s;
    });
  }

  function fTrue(x) {
    return 2 * Math.sin(5 * x);
  }

  // Pad the observed min/max of a set of curves (and possibly extra points)
  // by ~10% to get a sensible ylim, since the original script leaves ylim
  // to matplotlib's auto-ranging.
  function computeYlim(curves, extraPoints, padFrac) {
    let lo = Infinity,
      hi = -Infinity;
    for (const c of curves) {
      for (const v of c) {
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    }
    if (extraPoints) {
      for (const v of extraPoints) {
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    }
    if (!isFinite(lo) || !isFinite(hi)) {
      lo = -1;
      hi = 1;
    }
    const pad = (hi - lo) * (padFrac !== undefined ? padFrac : 0.1) || 1;
    return [lo - pad, hi + pad];
  }

  // ---------------------------------------------------------------
  // Figure 7.1(a): prior sample functions f(x) = para . phi(x),
  // para ~ N(0, I_10)
  // ---------------------------------------------------------------
  {
    const zeroCurve = new Array(nX).fill(0);
    const sampleCurves = [];
    for (let i = 0; i < nSampleA; i++) {
      const para = window.rnd.randnVec(nF);
      sampleCurves.push(evalCurve(para, xP));
    }
    const ylim = computeYlim(sampleCurves.concat([zeroCurve]));

    const body = window.plotlib.makeCard(outputGrid, "Figure 7.1(a) — Prior sample functions");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, 1],
      ylim: ylim,
      xlabel: "\\rm x",
      ylabel: "f({\\rm x})",
    });

    chart.line(xP, zeroCurve, { lineWidth: 3 });
    for (const fx of sampleCurves) {
      chart.line(xP, fx, { lineWidth: 0.5, color: [0.7, 0.7, 0.7] });
    }
    chart.finish();
  }

  // ---------------------------------------------------------------
  // Figure 7.1(b): Bayesian linear regression posterior, three noise
  // variance settings sigma_sq = [0.5, 100, 1e-6] (fixed, not parametrized)
  // ---------------------------------------------------------------
  {
    const sigmaSq = [0.5, 100, 1e-6];

    // x_data: s_bar evenly spaced points over [0,1] (matches Python's
    // `x = np.linspace(0, 1, s_bar)` used as the data locations).
    const xData = Array.from({ length: sBar }, (_, i) => (sBar === 1 ? 0 : i / (sBar - 1)));

    // Phi: n_f x s_bar (features as ROWS, data points as COLUMNS) --
    // matches numpy's np.column_stack([phi(x_s) for x_s in x]) shape exactly.
    const Phi = Array.from({ length: nF }, (_, k) => xData.map((x) => phi(x)[k]));
    const PhiT = window.linalg.transpose(Phi); // s_bar x n_f

    const y = xData.map((x) => fTrue(x) + window.rnd.randn());

    const means = [];
    const covs = [];
    for (let l = 0; l < sigmaSq.length; l++) {
      const A = window.linalg.matMul(PhiT, Phi); // s_bar x s_bar  (Phi.T @ Phi)
      const Areg = window.linalg.addDiag(A, sigmaSq[l]); // + sigma_sq[l] * I
      const Ainv = window.linalg.inv(Areg);
      const PhiAinv = window.linalg.matMul(Phi, Ainv); // n_f x s_bar
      const tmp0 = window.linalg.matMul(PhiAinv, PhiT); // n_f x n_f = Phi @ Ainv @ Phi.T
      const eyeN = window.linalg.eye(nF);
      const tmp = eyeN.map((row, i) => row.map((v, j) => v - tmp0[i][j]));
      const tmpT = window.linalg.transpose(tmp);
      const cov = tmp.map((row, i) => row.map((v, j) => (v + tmpT[i][j]) / 2));
      const mean = window.linalg.matVec(Phi, window.linalg.matVec(Ainv, y)); // Phi @ Ainv @ y

      means.push(mean);
      covs.push(cov);
    }

    const meanCurves = [
      evalCurve(means[0], xP),
      evalCurve(means[1], xP),
      evalCurve(means[2], xP),
    ];

    const sampleCurves = [];
    for (let i = 0; i < nSampleB; i++) {
      const para = window.rnd.mvnSample(means[0], covs[0]);
      sampleCurves.push(evalCurve(para, xP));
    }

    const ylim = computeYlim(meanCurves.concat(sampleCurves), y);

    const body = window.plotlib.makeCard(outputGrid, "Figure 7.1(b) — Posterior sample functions");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, 1],
      ylim: ylim,
      xlabel: "\\rm x",
      ylabel: "f({\\rm x})",
      legendLoc: "upper right",
    });

    chart.line(xP, meanCurves[0], { lineWidth: 2, color: "blue", label: "\\sigma=0.5" });
    chart.line(xP, meanCurves[1], { lineWidth: 2, color: "orange", label: "\\sigma=10" });
    chart.line(xP, meanCurves[2], { lineWidth: 2, color: "green", label: "\\sigma=10^{-3}" });

    for (const fx of sampleCurves) {
      chart.line(xP, fx, { lineWidth: 0.3, color: [0.7, 0.7, 0.7] });
    }

    chart.scatter(xData, y, { marker: "circle", color: "blue", filled: false, label: "{\\rm y}_s" });

    chart.finish();
  }
};
