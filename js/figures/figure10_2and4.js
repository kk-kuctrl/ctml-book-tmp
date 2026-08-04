"use strict";

// Port of code_distribute/python/Figure10_2to4.py (figure10_2b + figure10_4).
//
// scipy.optimize.minimize -> both L_KL and L_IRL are convex in v (each is a
// sum/combination of an affine term and a log-sum-exp of affine terms, which
// is convex), so instead of a general-purpose QP/NLP solver they're solved
// here with plain gradient descent (Adam) on a hand-derived analytic
// gradient -- no numerical-diff, no external optimizer needed.
window.figureLib = window.figureLib || {};
window.figureLib.figure10_2and4 = function (outputGrid, params) {
  const L = window.linalg;
  const beta = 0.8;

  // P is free-text (like Figure 6.1/8.8's A) and any square size -- its own
  // row count IS the number of states. cost must then match that size,
  // falling back to a generated 1..n (rather than the fixed 4-dim default)
  // if it doesn't parse at n's current size.
  const DEFAULT_P_TEXT = "0.3333,0.3333,0,0\n0,0.3333,0.3333,0\n0,0.3333,0.3333,0.3333\n0.6667,0,0.3333,0.6667";
  function parseSquareMatrixAny(text) {
    const rows = String(text)
      .trim()
      .split(/[\n;]+/)
      .map((r) => r.trim())
      .filter((r) => r.length > 0);
    const m = rows.length;
    if (m === 0) return null;
    const M = rows.map((r) =>
      r
        .split(/[,\s]+/)
        .filter((s) => s.length > 0)
        .map(Number)
    );
    if (M.some((row) => row.length !== m || row.some((v) => !Number.isFinite(v)))) return null;
    return M;
  }
  function parseVecFixed(text, m) {
    const v = String(text)
      .trim()
      .split(/[,\s]+/)
      .filter((s) => s.length > 0)
      .map(Number);
    if (v.length !== m || v.some((x) => !Number.isFinite(x))) return null;
    return v;
  }
  function defaultCost(m) {
    return Array.from({ length: m }, (_, i) => i + 1);
  }

  let Praw = parseSquareMatrixAny(params.P_text);
  if (!Praw) Praw = parseSquareMatrixAny(DEFAULT_P_TEXT);
  const n = Praw.length;

  let cost = parseVecFixed(params.cost_text, n);
  if (!cost) cost = defaultCost(n);

  // Nominal transition matrix P(i,j) = Prob(current=j -> next=i), built from
  // the typed-in raw weight matrix: each from-state j's column of weights is
  // renormalized to sum to 1, so any typed matrix gives a valid stochastic
  // matrix -- a 0 entry removes that edge entirely, so the topology itself
  // (not just the probabilities) is fully adjustable. If every weight in a
  // column is 0 (all outgoing edges removed), fall back to a uniform column
  // so P never degenerates into an invalid (all-zero) column.
  function buildNominalP(raw) {
    const P = L.zeros(n, n);
    for (let j = 0; j < n; j++) {
      let total = 0;
      for (let i = 0; i < n; i++) total += raw[i][j];
      for (let i = 0; i < n; i++) P[i][j] = total <= 1e-9 ? 1 / n : raw[i][j] / total;
    }
    return P;
  }
  const P = buildNominalP(Praw);

  // w[a] = sum_b P[b][a] * z[b] = (P^T @ z)[a], shared by both gradients below.
  function wOf(z) {
    const w = new Array(n).fill(0);
    for (let a = 0; a < n; a++) {
      let s = 0;
      for (let b = 0; b < n; b++) s += P[b][a] * z[b];
      w[a] = s;
    }
    return w;
  }

  // grad of L_KL(v) = sum_a (v_a - cost_a + log(w_a(v)))^2, derived by hand:
  //   r_a = v_a - cost_a + log(w_a),  dr_a/dv_k = delta(a,k) - beta*z_k*P[k,a]/w_a
  //   dL/dv_k = 2*r_k - 2*beta*z_k * sum_a (P[k,a]/w_a) * r_a
  function gradLKL(v) {
    const z = v.map((vi) => Math.exp(-beta * vi));
    const w = wOf(z);
    const r = v.map((vi, a) => vi - cost[a] + Math.log(w[a]));
    const grad = new Array(n).fill(0);
    for (let k = 0; k < n; k++) {
      let s = 0;
      for (let a = 0; a < n; a++) s += (P[k][a] / w[a]) * r[a];
      grad[k] = 2 * r[k] - 2 * beta * z[k] * s;
    }
    return grad;
  }

  // grad of L_IRL(v) = beta*(cntA . v) + cntB . log(w(v)), same w_a(v) as above:
  //   dL/dv_k = beta*cntA_k - beta*z_k * sum_a cntB_a*P[k,a]/w_a
  function gradLIRL(v, cntA, cntB) {
    const z = v.map((vi) => Math.exp(-beta * vi));
    const w = wOf(z);
    const grad = new Array(n).fill(0);
    for (let k = 0; k < n; k++) {
      let s = 0;
      for (let a = 0; a < n; a++) s += (cntB[a] * P[k][a]) / w[a];
      grad[k] = beta * cntA[k] - beta * z[k] * s;
    }
    return grad;
  }

  // Adam on an analytic gradient. `mask` freezes coordinates (used to pin
  // v[0] = 0, matching the Python's equality constraint x[0] = offset = 0).
  function adamMinimize(gradFn, x0, iters, lr, mask) {
    const dim = x0.length;
    const m = mask || new Array(dim).fill(true);
    const x = x0.slice();
    const mm = new Array(dim).fill(0);
    const vv = new Array(dim).fill(0);
    const b1 = 0.9,
      b2 = 0.999,
      eps = 1e-8;
    for (let t = 1; t <= iters; t++) {
      const g = gradFn(x);
      for (let i = 0; i < dim; i++) {
        if (!m[i]) continue;
        mm[i] = b1 * mm[i] + (1 - b1) * g[i];
        vv[i] = b2 * vv[i] + (1 - b2) * g[i] * g[i];
        const mHat = mm[i] / (1 - Math.pow(b1, t));
        const vHat = vv[i] / (1 - Math.pow(b2, t));
        x[i] -= (lr * mHat) / (Math.sqrt(vHat) + eps);
      }
    }
    return x;
  }

  function cumsumCols(M) {
    const out = L.zeros(n, n);
    for (let j = 0; j < n; j++) {
      let acc = 0;
      for (let i = 0; i < n; i++) {
        acc += M[i][j];
        out[i][j] = acc;
      }
    }
    return out;
  }
  function simulateChain(Pchain, kBar, startState) {
    const Paccum = cumsumCols(Pchain);
    const states = new Array(kBar + 1).fill(0);
    states[0] = startState;
    for (let k = 1; k <= kBar; k++) {
      const u = Math.random();
      const col = states[k - 1] - 1;
      let next = n;
      for (let i = 0; i < n; i++) {
        if (u <= Paccum[i][col]) {
          next = i + 1;
          break;
        }
      }
      states[k] = next;
    }
    return states;
  }
  function stairsXY(values) {
    const xs = [],
      ys = [];
    for (let k = 0; k < values.length; k++) {
      xs.push(k, k + 1);
      ys.push(values[k], values[k]);
    }
    return { xs, ys };
  }

  // Stationary distribution of a column-stochastic P, via power iteration
  // (repeatedly apply P to a uniform start) instead of solving (I-P)x=0
  // directly -- with the topology now fully user-adjustable, P can easily
  // become reducible or periodic (e.g. a 2-cycle), where (I-P) is singular
  // or the plain iterate never settles. Averaging the trailing iterates
  // washes out periodic oscillation while still converging to the true
  // stationary distribution in the ordinary (irreducible, aperiodic) case.
  function stationaryDistribution(Pin) {
    let p = new Array(n).fill(1 / n);
    const tailLen = 30;
    const tail = [];
    for (let k = 0; k < 400; k++) {
      p = L.matVec(Pin, p);
      if (k >= 400 - tailLen) tail.push(p.slice());
    }
    const avg = new Array(n).fill(0);
    for (const t of tail) for (let i = 0; i < n; i++) avg[i] += t[i] / tail.length;
    const s = avg.reduce((a, b) => a + b, 0) || 1;
    return avg.map((v) => v / s);
  }

  // ---- Figure 10.2(b): sample path of the autonomous chain under P ----
  const kBar2 = Math.max(10, Math.round(params.k_bar_2));
  {
    const states = simulateChain(P, kBar2, n);
    const { xs, ys } = stairsXY(states);
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.2(b) — ノミナル遷移行列によるサンプル軌道");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar2], ylim: [0.8, n + 0.2], xlabel: "$k$" });
    chart.line(xs, ys, { color: "blue", lineWidth: 2 });
    chart.finish();
  }

  // ---- Figure 10.4: KL-control optimal policy + online IRL cost recovery ----
  // Step 1: solve L_KL(v) = 0 once for the optimal value function, then get
  // the KL-optimal transition matrix P_opt in closed form (same simplification
  // as the Python's P_opt[i,j] = P[i,j]*z_opt[i]/w(j), just precomputed).
  const vStar = adamMinimize(gradLKL, new Array(n).fill(5), 5000, 0.05, null);
  const zOpt = vStar.map((vi) => Math.exp(-beta * vi));
  const wOpt = wOf(zOpt);
  const Popt = L.zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Popt[i][j] = (P[i][j] * zOpt[i]) / wOpt[j];

  {
    const pStationary = stationaryDistribution(P);
    const pOptStationary = stationaryDistribution(Popt);
    // Long-run average cost actually incurred under each policy -- the
    // expectation of cost(state) against that policy's own stationary
    // distribution (nominal P just drifts on its own; P_opt is what KL
    // control steers it to).
    const nominalCost = pStationary.reduce((s, p, i) => s + p * cost[i], 0);
    const optimalCost = pOptStationary.reduce((s, p, i) => s + p * cost[i], 0);
    const fmtVec = (v) => "(" + v.map((x) => x.toFixed(3)).join(", ") + ")";
    const rows = [
      ["ノミナル P の定常分布 π", fmtVec(pStationary)],
      ["最適 P* の定常分布 π*", fmtVec(pOptStationary)],
      ["ノミナルコスト（π の下でのコストの定常平均）", nominalCost.toFixed(4)],
      ["最適コスト（π* の下でのコストの定常平均）", optimalCost.toFixed(4)],
    ];
    const body = window.plotlib.makeCard(outputGrid, "定常分布と最適コスト");
    for (const [label, value] of rows) {
      const row = document.createElement("p");
      row.className = "hint";
      row.style.margin = "4px 0";
      row.textContent = `${label}: ${value}`;
      body.appendChild(row);
    }
  }

  // Step 2: simulate the optimally-controlled chain while an IRL observer
  // (who only knows the nominal model P, not the true cost) re-estimates the
  // cost vector every step from state-visitation counts alone -- each step
  // re-solves L_IRL from scratch (not warm-started), matching the Python.
  const kBar4 = Math.max(50, Math.round(params.k_bar_4));
  const Paccum4 = cumsumCols(Popt);
  const states4 = new Array(kBar4 + 1).fill(0);
  states4[0] = n;
  const cntA = new Array(n).fill(0);
  cntA[n - 1] = 1; // visits into each next-state (the initial state counts as one visit)
  const cntB = new Array(n).fill(0); // visits out of each current-state
  const invLHist = Array.from({ length: n }, () => new Array(kBar4));
  const irlMask = Array.from({ length: n }, (_, i) => i !== 0); // v[0] pinned to offset = 0
  const irlV0 = Array.from({ length: n }, (_, i) => (i === 0 ? 0 : 5));
  for (let k = 0; k < kBar4; k++) {
    const v = adamMinimize((vv) => gradLIRL(vv, cntA, cntB), irlV0, 400, 0.1, irlMask);
    const invZ = v.map((vi) => Math.exp(-beta * vi));
    const w = wOf(invZ);
    for (let a = 0; a < n; a++) invLHist[a][k] = v[a] + Math.log(w[a]);

    const u = Math.random();
    const col = states4[k] - 1;
    let next = n;
    for (let i = 0; i < n; i++) {
      if (u <= Paccum4[i][col]) {
        next = i + 1;
        break;
      }
    }
    states4[k + 1] = next;
    cntA[next - 1] += 1;
    cntB[states4[k] - 1] += 1;
  }
  // Cost is only identifiable up to an additive shift here, so every column
  // is re-anchored to ell_1 = 1 (matches inv_l_hist - inv_l_hist[0,:] + 1).
  const normHist = invLHist.map((row, a) => row.map((v, k) => v - invLHist[0][k] + 1));

  {
    const kShow = Math.min(50, kBar4);
    const states4head = states4.slice(0, kShow + 1);
    const { xs, ys } = stairsXY(states4head);
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.4(a) — KL制御最適方策によるサンプル軌道");
    const chart = window.plotlib.createChart(body, { xlim: [0, kShow], ylim: [0.8, n + 0.2], xlabel: "$k$" });
    chart.line(xs, ys, { color: "blue", lineWidth: 2 });
    chart.finish();
  }
  {
    // The curves should settle near cost_a - cost_1 + 1 (identifiable only up
    // to the additive shift that pins ell_1 = 1), but the first few steps --
    // solved from near-empty visitation counts -- can swing wildly before
    // enough data accumulates. Sizing the axis off the actual min/max would
    // let that early transient dominate and flatten the converged region
    // (which is the part worth looking at), so size it off where the curves
    // are expected to land instead -- same idea as the Python's fixed
    // [0, 6.5], generalized to track wherever the cost sliders put it.
    const targets = cost.map((c) => c - cost[0] + 1);
    const tMin = Math.min(...targets, 1);
    const tMax = Math.max(...targets, 1);
    const ks = Array.from({ length: kBar4 }, (_, i) => i);
    const palette = ["blue", "orange", "green", "red"];
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.4(b) — IRLによる状態コストの推定値の推移");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar4], ylim: [tMin - 1, tMax + 2.5], xlabel: "$k$" });
    for (let a = 0; a < n; a++) chart.line(ks, normHist[a], { color: palette[a % palette.length], lineWidth: 1.5, label: `\\ell_${a + 1}` });
    chart.finish();
  }
};
