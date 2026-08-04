"use strict";

// Port of code_distribute/python/Figure5_4.py
//
// The textbook's original is a QP with an auxiliary variable t_i >=
// |residual_i| per time step, minimizing sum(t_i) (the standard epigraph
// trick for an L1 penalty, since w is assumed Laplace-distributed -> L1 is
// the MLE-optimal residual loss). At the optimum t_i always equals
// |residual_i| exactly, so t can be eliminated analytically and the default
// problem reduces to an unconstrained-except-for-a-box-on-u minimization:
//
//   min_{x0, u}  ||x0-mu0||^2/2 + ||u||^2/2 + sum_i |Aux1_i.x0 + Aux2_i.u - y_i|
//   s.t. -1 <= u_k <= 1
//
// Since x0/u/w's true distributions are now user-selectable (gaussian /
// laplace / uniform each), each term's penalty is generalized to that
// distribution's own negative-log-likelihood (gaussian -> L2, laplace ->
// L1, uniform -> a hard box constraint with no smooth penalty at all) --
// see nllGrad below -- solved throughout with projected Adam on a smoothed
// absolute value where needed (sqrt(r^2+eps), whose gradient
// r/sqrt(r^2+eps) -> sign(r) as eps -> 0), instead of a general QP/SOCP
// solver.
window.figureLib = window.figureLib || {};
window.figureLib.figure5_4 = function (outputGrid, params) {
  const kBar = Math.max(10, Math.round(params.k_bar));

  const A = [
    [-0.39, -0.67, -0.34],
    [0.71, -0.51, 0.11],
    [-0.46, -0.35, -0.12],
  ];
  const Bv = [0, 1, 0]; // flattened (x_dim x 1)
  const C = [1, 0, 0]; // flattened (1 x x_dim)
  const mu0 = [params.init_mean, params.init_mean, params.init_mean];
  const xDim = 3;

  // Draws one scalar sample from the chosen distribution, centered at
  // `mean`. `scale` plays the role of the support's half-width for
  // "uniform" and of the VARIANCE for "gaussian"/"laplace" (Laplace(b) has
  // variance 2b^2, so b is back-derived from the requested variance) --
  // this only changes how the TRUE x_0 / w_k are generated for the
  // simulation, not the (fixed, L2/L1) estimator objective below.
  function sampleDist(kind, scale, mean) {
    if (kind === "uniform") return window.rnd.uniformSample(mean - scale, mean + scale);
    if (kind === "laplace") return mean + window.rnd.laplaceSample(Math.sqrt(scale / 2));
    return mean + Math.sqrt(scale) * window.rnd.randn();
  }

  // Process noise v_k: unlike w (above), its support is a param INDEPENDENT
  // of the chosen shape -- "uniform" just uses it directly as the range
  // (variance/mean meaningless there), while "gaussian"/"laplace" draw from
  // the full (mean-0) distribution and then either truncate to
  // [-support, support] via rejection sampling (matching the original
  // truncNormSample(-1,1), which is exactly this with kind="gaussian",
  // scale=1, support=1) or, if proc_unbounded is set, skip truncation
  // entirely and return the raw draw.
  function sampleProcessNoise(kind, scale, support, unbounded) {
    if (kind === "uniform") return window.rnd.uniformSample(-support, support);
    const draw = () => (kind === "laplace" ? window.rnd.laplaceSample(Math.sqrt(scale / 2)) : Math.sqrt(scale) * window.rnd.randn());
    if (unbounded) return draw();
    for (let i = 0; i < 1000; i++) {
      const v = draw();
      if (v >= -support && v <= support) return v;
    }
    return 0; // fallback, practically unreachable for reasonable scale/support combos
  }

  function matVec3(M, v) {
    return [M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2], M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2], M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]];
  }
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  // ---------------------------------------------------------------
  // Generate the true trajectory + noisy measurements.
  // ---------------------------------------------------------------
  const x0True = mu0.map((m) => sampleDist(params.init_dist, params.init_scale, m));
  const xData = [x0True.slice()];
  const yData = new Array(kBar + 1).fill(0);
  const vData = new Array(kBar).fill(0);
  const wData = new Array(kBar + 1).fill(0);

  let x = x0True.slice();
  for (let k = 0; k < kBar; k++) {
    const vk = sampleProcessNoise(params.proc_dist, params.proc_scale, params.proc_support, params.proc_unbounded);
    const wk = sampleDist(params.noise_dist, params.noise_scale, 0);
    yData[k] = dot(C, x) + wk;
    vData[k] = vk;
    wData[k] = wk;
    x = matVec3(A, x).map((val, i) => val + Bv[i] * vk);
    xData.push(x.slice());
  }
  wData[kBar] = sampleDist(params.noise_dist, params.noise_scale, 0);
  yData[kBar] = dot(C, x) + wData[kBar];

  // ---------------------------------------------------------------
  // Aux1[i] = C @ A^i  (row, length x_dim)
  // Aux2[i][j] = C @ A^(i-1-j) @ Bv  for j < i, else 0
  // ---------------------------------------------------------------
  function rowTimesMat(rowVec, M) {
    // rowVec (1xn) @ M (nxn) -> row (1xn)
    const n = M.length;
    const out = new Array(n).fill(0);
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += rowVec[k] * M[k][j];
      out[j] = s;
    }
    return out;
  }

  const Aux1 = [];
  const Aux2 = [];
  const powers = [window.linalg.eye(xDim)]; // powers[i] = A^i
  for (let i = 1; i <= kBar; i++) powers.push(window.linalg.matMul(powers[i - 1], A));

  for (let i = 0; i <= kBar; i++) Aux1.push(rowTimesMat(C, powers[i]));

  for (let i = 0; i <= kBar; i++) {
    const row = new Array(kBar).fill(0);
    for (let j = 0; j < i; j++) {
      const CAB = rowTimesMat(C, powers[i - 1 - j]); // C @ A^(i-1-j), length x_dim
      row[j] = dot(CAB, Bv);
    }
    Aux2.push(row);
  }

  // ---------------------------------------------------------------
  // Projected Adam on the reduced (t-eliminated) objective -- generalized
  // to a NEGATIVE LOG-LIKELIHOOD gradient per term, matching whichever
  // distribution was selected for x0 / w / v, instead of always assuming
  // the textbook's original (Gaussian prior on x0, Laplace/L1 likelihood
  // on w, Gaussian/L2 prior on u):
  //   gaussian: -log p(r) ~ r^2/(2*variance)      -> grad = r/variance
  //   laplace:  -log p(r) ~ |r-mean|/b, b=sqrt(variance/2) -> grad = sign(r)/b
  //   uniform:  -log p(r) is CONSTANT inside the support -> grad = 0 (no
  //             smooth signal at all; the support is instead enforced as a
  //             hard projection after each step, same idea as the box
  //             constraint on u the original problem already had)
  //
  // A single small fixed eps isn't enough for the laplace/L1 case: this
  // problem has more free parameters (x0: 3, u: 60) than residual
  // constraints (61), so the true L1 solution is SPARSE -- most residuals
  // land exactly on zero (the fit silently absorbs typical noise into x0/u)
  // while a few genuine outliers stay large. A smoothed |r| ~ sqrt(r^2+eps)
  // has no exact zero, so with one fixed eps the optimizer spreads small
  // corrections over every residual instead of concentrating on the
  // outliers. Annealing eps from loose to sharp (graduated non-convexity)
  // converges to a much better approximation of that sparse structure than
  // a single small eps does; gaussian terms just ignore eps2 entirely
  // (their gradient is already exactly smooth), so sharing one schedule
  // across all three distribution choices is harmless.
  function nllGrad(r, dist, variance, eps2) {
    if (dist === "uniform") return 0;
    if (dist === "laplace") return r / Math.sqrt(r * r + eps2) / Math.sqrt(variance / 2);
    return r / variance;
  }

  function gradients(x0v, u, eps2) {
    const gx0 = x0v.map((v, i) => nllGrad(v - mu0[i], params.init_dist, params.init_scale, eps2));
    const gu = u.map((uk) => nllGrad(uk, params.proc_dist, params.proc_scale, eps2));
    for (let i = 0; i <= kBar; i++) {
      const r = dot(Aux1[i], x0v) + dot(Aux2[i], u) - yData[i];
      const s = nllGrad(r, params.noise_dist, params.noise_scale, eps2);
      for (let j = 0; j < xDim; j++) gx0[j] += s * Aux1[i][j];
      for (let j = 0; j < kBar; j++) gu[j] += s * Aux2[i][j];
    }
    return { gx0, gu };
  }

  function makeAdam(n) {
    return { m: new Array(n).fill(0), v: new Array(n).fill(0), t: 0 };
  }
  function adamStep(x, g, state, lr) {
    state.t++;
    const b1 = 0.9,
      b2 = 0.999,
      eps = 1e-8;
    for (let i = 0; i < x.length; i++) {
      state.m[i] = b1 * state.m[i] + (1 - b1) * g[i];
      state.v[i] = b2 * state.v[i] + (1 - b2) * g[i] * g[i];
      const mHat = state.m[i] / (1 - Math.pow(b1, state.t));
      const vHat = state.v[i] / (1 - Math.pow(b2, state.t));
      x[i] -= (lr * mHat) / (Math.sqrt(vHat) + eps);
    }
  }

  // Whether to hard-project onto a bounded support after each step: always
  // for "uniform" (its NLL gradient is always 0, so nothing else would ever
  // pull the estimate back inside the support), and additionally for
  // gaussian/laplace on u unless proc_unbounded is checked (matching the
  // same "unbounded" semantics used for the true data generation above).
  const projectU = params.proc_dist === "uniform" || !params.proc_unbounded;
  const uSupport = params.proc_support;
  const projectX0 = params.init_dist === "uniform";

  let x0Est = mu0.slice();
  let uEst = new Array(kBar).fill(0);
  const adamX0 = makeAdam(xDim);
  const adamU = makeAdam(kBar);
  const nIter = 6000;
  const epsStart = 1e-1;
  const epsEnd = 1e-10;
  for (let it = 0; it < nIter; it++) {
    const frac = it / (nIter - 1);
    const eps2 = epsStart * Math.pow(epsEnd / epsStart, frac); // geometric anneal
    const lr = 0.05 * (1 - 0.9 * frac); // smaller steps as the objective sharpens
    const { gx0, gu } = gradients(x0Est, uEst, eps2);
    adamStep(x0Est, gx0, adamX0, lr);
    adamStep(uEst, gu, adamU, lr);
    if (projectU) for (let k = 0; k < kBar; k++) uEst[k] = Math.max(-uSupport, Math.min(uSupport, uEst[k]));
    if (projectX0) for (let j = 0; j < xDim; j++) x0Est[j] = Math.max(mu0[j] - params.init_scale, Math.min(mu0[j] + params.init_scale, x0Est[j]));
  }

  // Reconstruct the estimated state trajectory from x0Est, uEst.
  const xhat = [x0Est.slice()];
  let xk = x0Est.slice();
  for (let k = 0; k < kBar; k++) {
    xk = matVec3(A, xk).map((val, i) => val + Bv[i] * uEst[k]);
    xhat.push(xk.slice());
  }

  // ---------------------------------------------------------------
  // Plotting (matches the Python's figure5_4a-d)
  // ---------------------------------------------------------------
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

  const xMax = 6;
  const ks = Array.from({ length: kBar + 1 }, (_, k) => k);
  const col = (arr, j) => arr.map((row) => row[j]);

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.4(a) — 状態1の真値・観測・推定値");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar], ylim: [-xMax, xMax], xlabel: "$k$" });
    chart.line(ks, col(xData, 0), { color: "blue", lineWidth: 1.5, label: "\\text{True state } (x_k)_1" });
    chart.line(ks, yData, { color: "green", lineWidth: 1.5, dash: [6, 4], label: "\\text{Measurements } y" });
    chart.line(ks, col(xhat, 0), { color: "red", lineWidth: 1.5, dash: [2, 2], label: "\\text{Estimated } (\\hat{x}_k)_1" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.4(b) — 状態3の真値・推定値");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar], ylim: [-xMax, xMax], xlabel: "$k$" });
    chart.line(ks, col(xData, 2), { color: "blue", lineWidth: 1.5, label: "\\text{True state } (x_k)_3" });
    chart.line(ks, col(xhat, 2), { color: "red", lineWidth: 1.5, dash: [2, 2], label: "\\text{Estimated } (\\hat{x}_k)_3" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.4(c) — 外乱の真値・推定値");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar - 1], ylim: [-xMax, xMax], xlabel: "$k$" });
    const vS = stairsXY(vData);
    const uS = stairsXY(uEst);
    chart.line(vS.xs, vS.ys, { color: "blue", lineWidth: 1.5, label: "\\text{Disturbance } v_k" });
    chart.line(uS.xs, uS.ys, { color: "red", lineWidth: 1.5, dash: [2, 2], label: "\\text{Estimated } \\hat{v}_k:=u(k)" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.4(d) — 観測雑音の真値・推定値");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar], ylim: [-xMax, xMax], xlabel: "$k$" });
    const wS = stairsXY(wData);
    const westEst = yData.map((y, i) => y - xhat[i][0]);
    const wEstS = stairsXY(westEst);
    chart.line(wS.xs, wS.ys, { color: "blue", lineWidth: 1.5, label: "\\text{Noise } w_k" });
    chart.line(wEstS.xs, wEstS.ys, { color: "red", lineWidth: 1.5, dash: [2, 2], label: "\\text{Estimated } \\hat{w}_k" });
    chart.finish();
  }
};
