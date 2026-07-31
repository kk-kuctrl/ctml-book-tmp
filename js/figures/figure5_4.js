"use strict";

// Port of code_distribute/python/Figure5_4.py
//
// The original solves a QP with an auxiliary variable t_i >= |residual_i|
// per time step and minimizes sum(t_i) (the standard epigraph trick for an
// L1 penalty). At the optimum t_i always equals |residual_i| exactly (there
// is never a reason to make it larger), so t can be eliminated analytically
// and the problem reduces to an unconstrained-except-for-a-box-on-u convex
// minimization:
//
//   min_{x0, u}  ||x0-mu0||^2/2 + ||u||^2/2 + sum_i |Aux1_i.x0 + Aux2_i.u - y_i|
//   s.t. -1 <= u_k <= 1
//
// which we solve with projected Adam on a smoothed absolute value
// (sqrt(r^2+eps), whose gradient r/sqrt(r^2+eps) -> sign(r) as eps -> 0),
// instead of needing a general QP/SOCP solver.
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
  const mu0 = [1, 1, 1];
  const xDim = 3;

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
  const x0True = window.rnd.mvnSample(mu0, window.linalg.eye(xDim));
  const xData = [x0True.slice()];
  const yData = new Array(kBar + 1).fill(0);
  const vData = new Array(kBar).fill(0);
  const wData = new Array(kBar + 1).fill(0);

  let x = x0True.slice();
  for (let k = 0; k < kBar; k++) {
    const vk = window.rnd.truncNormSample(-1, 1);
    const wk = window.rnd.laplaceSample(1);
    yData[k] = dot(C, x) + wk;
    vData[k] = vk;
    wData[k] = wk;
    x = matVec3(A, x).map((val, i) => val + Bv[i] * vk);
    xData.push(x.slice());
  }
  wData[kBar] = window.rnd.laplaceSample(1);
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
  // Projected Adam on the reduced (t-eliminated) smooth-L1 objective.
  //
  // A single small fixed eps isn't enough: this problem has more free
  // parameters (x0: 3, u: 60) than residual constraints (61), so the true
  // L1 solution is SPARSE -- most residuals land exactly on zero (the fit
  // silently absorbs typical noise into x0/u) while a few genuine outliers
  // stay large. A smoothed |r| ~ sqrt(r^2+eps) has no exact zero, so with
  // one fixed eps the optimizer spreads small corrections over every
  // residual instead of concentrating on the outliers. Annealing eps from
  // loose to sharp (graduated non-convexity: solve the easy, nearly-quadratic
  // version first, then sharpen toward true L1) converges to a much better
  // approximation of that sparse structure than a single small eps does.
  function gradients(x0v, u, eps2) {
    const gx0 = x0v.map((v, i) => v - mu0[i]);
    const gu = u.slice();
    for (let i = 0; i <= kBar; i++) {
      const r = dot(Aux1[i], x0v) + dot(Aux2[i], u) - yData[i];
      const s = r / Math.sqrt(r * r + eps2);
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
    for (let k = 0; k < kBar; k++) uEst[k] = Math.max(-1, Math.min(1, uEst[k]));
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
