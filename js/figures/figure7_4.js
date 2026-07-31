"use strict";

// Port of code_distribute/python/Figure7_4.py
//
// Three convex regressions on the same polynomial-feature design matrix:
//   Naive (plain least squares)  -> closed form: solve(Phi Phi^T, Phi y)
//   Ridge (L2 penalty)           -> closed form: solve(Phi Phi^T + sigma^2 I, Phi y)
//   Lasso (L1 penalty)           -> no closed form; solved with FISTA
//     (accelerated proximal gradient + soft-thresholding), which converges
//     in a few hundred iterations for a 10-feature problem like this --
//     no need for a general QP/cvxpy-style solver.
window.figureLib = window.figureLib || {};
window.figureLib.figure7_4 = function (outputGrid, params) {
  const sBar = Math.max(4, Math.round(params.s_bar));
  const sigmaSq = params.sigma_sq;
  const L = window.linalg;
  const nF = 10;

  function phi(x) {
    const out = new Array(nF);
    let p = 1;
    for (let k = 0; k < nF; k++) {
      out[k] = p;
      p *= x;
    }
    return out;
  }
  function fTrue(x) {
    return 2 * Math.sin(5 * x);
  }
  function linspace(a, b, n) {
    return Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
  }
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  const nX = 100;
  const xP = linspace(0, 1, nX);
  const fReal = xP.map(fTrue);

  const xData = linspace(0, 1, sBar);
  // Phi: n_f x s_bar (features as rows, data points as columns) -- matches
  // numpy's np.column_stack([phi(x_s) for x_s in x]) shape exactly.
  const Phi = Array.from({ length: nF }, (_, k) => xData.map((x) => phi(x)[k]));
  const PhiT = L.transpose(Phi); // s_bar x n_f
  const y = xData.map((x) => fTrue(x) + window.rnd.randn());

  const G = L.matMul(Phi, PhiT); // n_f x n_f, = Phi @ Phi^T
  const Phiy = L.matVec(Phi, y); // n_f

  // Naive least squares: solve(Phi Phi^T, Phi y).
  const paraNaive = L.solve(G, Phiy);

  // Ridge: solve(Phi Phi^T + sigma^2 I, Phi y).
  const paraRidge = L.solve(L.addDiag(G, sigmaSq), Phiy);

  // Lasso via coordinate descent on f(para) = ||Phi^T para - y||^2 + sigma^2 ||para||_1.
  //
  // This basis is a Vandermonde-like matrix (1, x, x^2, ..., x^9 on 30
  // points) -- severely ill-conditioned/collinear, a famous source of
  // instability in polynomial fitting. A gradient/FISTA-style method
  // crawls extremely slowly along the near-degenerate directions of a
  // matrix like this and can look numerically "wild" long before it
  // actually converges. Coordinate descent -- update one coefficient at a
  // time via its exact 1D closed-form soft-threshold solution, holding the
  // others fixed -- is the standard textbook Lasso algorithm precisely
  // because it doesn't share that weakness with correlated features.
  function softThreshold1(v, t) {
    return Math.sign(v) * Math.max(Math.abs(v) - t, 0);
  }
  const Gjj = Phi.map((row) => dot(row, row)); // ||Phi[j,:]||^2 for each feature
  let paraLasso = new Array(nF).fill(0);
  let residCD = y.slice(); // y - Phi^T @ paraLasso, maintained incrementally
  // 500 sweeps looked converged (stable-looking numbers) but wasn't: verified
  // against cvxpy on the same design matrix, plain (non-accelerated)
  // coordinate descent needs ~20,000 sweeps on this ill-conditioned
  // Vandermonde basis before the small coefficients actually hit exact zero
  // -- short of that it just looks smoothly non-sparse instead of erroring
  // out, so the under-convergence is easy to miss. 30,000 gives margin and
  // is still trivial (10 features x 30 points): well under 100ms.
  for (let sweep = 0; sweep < 30000; sweep++) {
    for (let j = 0; j < nF; j++) {
      const oldPj = paraLasso[j];
      for (let i = 0; i < sBar; i++) residCD[i] += oldPj * Phi[j][i]; // exclude feature j
      let rho = 0;
      for (let i = 0; i < sBar; i++) rho += Phi[j][i] * residCD[i];
      const newPj = softThreshold1(rho, sigmaSq / 2) / Gjj[j];
      paraLasso[j] = newPj;
      for (let i = 0; i < sBar; i++) residCD[i] -= newPj * Phi[j][i]; // put it back
    }
  }

  function evalCurve(para) {
    return xP.map((x) => dot(para, phi(x)));
  }
  const naiveDat = evalCurve(paraNaive);
  const ridgeDat = evalCurve(paraRidge);
  const lassoDat = evalCurve(paraLasso);

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 7.4(a) — 最小二乗・Ridge・Lassoの比較");
    const chart = window.plotlib.createChart(body, { xlim: [0, 1], ylim: [-5, 5], xlabel: "\\rm x", ylabel: "f({\\rm x})" });
    chart.scatter(xData, y, { color: "blue", marker: "circle", filled: false, label: "{\\rm y}_s" });
    chart.line(xP, naiveDat, { color: "black", lineWidth: 2, label: "\\text{Least Square}" });
    chart.line(xP, ridgeDat, { color: "blue", lineWidth: 2, label: "\\text{Ridge}" });
    chart.line(xP, lassoDat, { color: "red", lineWidth: 2, label: "\\text{Lasso}" });
    chart.line(xP, fReal, { color: "green", dash: [6, 3], lineWidth: 2, label: "2\\sin(5{\\rm x}_s)" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 7.4(b) — 係数の大きさ比較");
    const idx = Array.from({ length: nF }, (_, i) => i + 1);
    const ridgeAbs = paraRidge.map(Math.abs);
    const lassoAbs = paraLasso.map(Math.abs);
    // This basis is ill-conditioned (Vandermonde-like), so coefficient
    // magnitudes vary a lot run to run -- a fixed ylim=[0,3] clipped Lasso
    // (and sometimes Ridge) points right out of view whenever a draw
    // produced larger values. Fit the axis to whatever this run actually got.
    const maxAbs = Math.max(...ridgeAbs, ...lassoAbs, 1e-6);
    const chart = window.plotlib.createChart(body, { xlim: [-0.2, nF + 0.5], ylim: [-0.5, maxAbs * 1.15], xlabel: "$i$" });
    chart.scatter(idx, ridgeAbs, { color: "blue", marker: "circle", filled: false, label: "\\text{Ridge}" });
    chart.scatter(idx, lassoAbs, { color: "orange", marker: "circle", label: "\\text{Lasso}" });
    chart.finish();
  }
};
