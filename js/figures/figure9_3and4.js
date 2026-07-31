"use strict";

// Port of code_distribute/python/Figure9_3and4.py
//
// control.dlqr -> same Riccati-to-convergence trick as figure6_1.js.
// scipy's solve_discrete_lyapunov -> linalg.discreteLyapunov. The
// eigenvalue-based instability check -> linalg.spectralRadius. The TD/RLS
// update for the Q-function parameters is the same recursive-least-squares
// pattern as figure8_2.js's sysid_module, just with a different regressor.
window.figureLib = window.figureLib || {};
window.figureLib.figure9_3and4 = function (outputGrid, params) {
  const kUpdate = Math.max(5, Math.round(params.k_update));
  const nPath = Math.max(2, Math.round(params.n_path));
  const L = window.linalg;

  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function outer(a, b) {
    return a.map((ai) => b.map((bj) => ai * bj));
  }
  function vecNorm(v) {
    return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
  }

  const beta = 0.95;
  const A = [
    [0.8, 0.9, 0.86],
    [0.3, 0.25, 1.0],
    [0.1, 0.55, 0.5],
  ];
  const Bvec = [1, 0, 0];
  const xDim = 3,
    uDim = 1,
    nAug = xDim + uDim;
  const pDim = (nAug * (nAug + 1)) / 2 + 1;

  // K_opt: converge the discounted Riccati recursion for (sqrt(beta)A,
  // sqrt(beta)B) -- once pre-scaled by sqrt(beta), this is an ordinary
  // (undiscounted) DARE, same trick as figure6_1.js's P_opt.
  function lqrConverge(Ain, Bin) {
    let PI = L.zeros(xDim, xDim);
    let K = new Array(xDim).fill(0);
    for (let it = 0; it < 3000; it++) {
      const PIB = L.matVec(PI, Bin);
      const Rt = 1 + dot(Bin, PIB); // R = I(1x1) = 1
      const St = L.matVec(L.transpose(Ain), PIB);
      const Qt = L.add(L.eye(xDim), L.matMul(L.transpose(Ain), L.matMul(PI, Ain))); // Q = I
      K = St.map((v) => v / Rt);
      PI = L.sub(Qt, outer(St, K));
    }
    return K;
  }
  const sqrtBeta = Math.sqrt(beta);
  const Kopt = lqrConverge(
    L.scale(A, sqrtBeta),
    Bvec.map((v) => v * sqrtBeta)
  );

  // phi(x, u): upper-triangular (incl. diagonal) entries of [x;u][x;u]^T,
  // row-major, plus a constant 1 bias term -- eq (9.35).
  function phi(x, u) {
    const xu = x.concat([u]); // length nAug
    const out = [];
    for (let i = 0; i < nAug; i++) for (let j = i; j < nAug; j++) out.push(xu[i] * xu[j]);
    out.push(1.0);
    return out;
  }
  // Inverse packing: p (length pDim) -> symmetric nAug x nAug matrix.
  function pToUps(p) {
    const H = L.zeros(nAug, nAug);
    let idx = 0;
    for (let i = 0; i < nAug; i++) {
      for (let j = i; j < nAug; j++) {
        H[i][j] = p[idx++];
      }
    }
    const Ht = L.transpose(H);
    return H.map((row, i) => row.map((v, j) => (v + Ht[i][j]) / 2));
  }

  const KIni = [4.11, 11.7519, 19.2184];

  function rl4lqr(sigma, iterGain) {
    const xNormHistP = [];
    const KErrHistP = [];
    const UpsErrHistP = [];

    for (let j = 0; j < nPath; j++) {
      let K = KIni.slice();
      let x = new Array(xDim).fill(0);
      const xHist = [x.slice()];
      const KErrHist = [];
      const UpsErrHist = [];
      let unstable = false;

      for (let g = 0; g < iterGain; g++) {
        const Acl = L.sub(A, outer(Bvec, K));
        const MtM = L.add(L.eye(xDim), outer(K, K)); // [I,-K^T] @ [I;-K] = I + K^T K
        const Pi = L.discreteLyapunov(L.scale(L.transpose(Acl), sqrtBeta), MtM);

        // AB = [A, B] (xDim x nAug); Ups_true = beta*AB^T@Pi@AB + blockdiag(Q,R)
        const AB = A.map((row, i) => row.concat([Bvec[i]]));
        const ABt = L.transpose(AB);
        const UpsTrue = L.scale(L.matMul(ABt, L.matMul(Pi, AB)), beta);
        for (let i = 0; i < xDim; i++) UpsTrue[i][i] += 1; // + Q (=I3)
        UpsTrue[xDim][xDim] += 1; // + R (=1)

        let p = new Array(pDim).fill(0);
        let Sigma = L.scale(L.eye(pDim), 10);

        for (let k = 0; k < kUpdate; k++) {
          KErrHist.push(vecNorm(K.map((v, i) => v - Kopt[i])) / vecNorm(Kopt));
          const upsDiff = [];
          const Ups = pToUps(p);
          for (let i = 0; i < nAug; i++) for (let jj = 0; jj < nAug; jj++) upsDiff.push(Ups[i][jj] - UpsTrue[i][jj]);
          const upsTrueFlat = [];
          for (let i = 0; i < nAug; i++) for (let jj = 0; jj < nAug; jj++) upsTrueFlat.push(UpsTrue[i][jj]);
          UpsErrHist.push(vecNorm(upsDiff) / vecNorm(upsTrueFlat));

          const u = -dot(K, x) + window.rnd.randn() * sigma;
          const cost = dot(x, x) + u * u; // Q=I, R=1

          const phiPre = phi(x, u);
          const xNext = L.matVec(A, x).map((v, i) => v + Bvec[i] * u);
          const uNext = -dot(K, xNext);
          const phiNext = phi(xNext, uNext);
          const q = phiPre.map((v, i) => v - beta * phiNext[i]);

          const Sq = L.matVec(Sigma, q);
          const denom = 1 + dot(q, Sq);
          const innov = cost - dot(q, p);
          p = p.map((v, i) => v + (Sq[i] * innov) / denom);
          const SqOuter = outer(Sq, Sq);
          Sigma = Sigma.map((row, i) => row.map((v, jj) => v - SqOuter[i][jj] / denom));

          x = xNext;
          xHist.push(x.slice());
        }

        const Ups = pToUps(p);
        const Uuu = Ups[xDim][xDim];
        K = Ups[xDim].slice(0, xDim).map((v) => v / Uuu);

        if (L.spectralRadius(L.sub(A, outer(Bvec, K))) > 1) {
          unstable = true;
          break;
        }
      }

      if (!unstable) {
        xNormHistP.push(xHist.map(vecNorm));
        KErrHistP.push(KErrHist);
        UpsErrHistP.push(UpsErrHist);
      }
    }
    return { xNormHistP, KErrHistP, UpsErrHistP, iterGain };
  }

  function meanCurve(rows) {
    const n = rows[0].length;
    return Array.from({ length: n }, (_, i) => rows.reduce((s, r) => s + r[i], 0) / rows.length);
  }

  function plotScenario(title, sigma, iterGain) {
    const { xNormHistP, KErrHistP, UpsErrHistP } = rl4lqr(sigma, iterGain);
    if (xNormHistP.length === 0) {
      // Every one of n_path runs hit an unstable gain after some update and
      // got dropped (matches the Python's `if not IsUnstable: append(...)`)
      // -- with nothing left to plot, say so instead of silently showing
      // no cards at all for this scenario.
      const body = window.plotlib.makeCard(outputGrid, `${title} — 表示なし`);
      const msg = document.createElement("p");
      msg.className = "hint";
      msg.style.margin = "0";
      msg.textContent = `全 ${nPath} パスが更新後に不安定なゲインとなったため、表示できるデータがありません（探索雑音レベルを上げるか、ステップ数を増やしてみてください）。`;
      body.appendChild(msg);
      return;
    }

    const nT = xNormHistP[0].length;
    const ks = Array.from({ length: nT }, (_, i) => i);

    {
      // ylim fit to whatever this run's exploration noise actually produced,
      // instead of a fixed [-1,8] that could clip (or look loosely padded
      // relative to) large sigma / iter_gain combinations.
      let normMax = 0;
      for (const row of xNormHistP) for (const v of row) if (v > normMax) normMax = v;
      const body = window.plotlib.makeCard(outputGrid, `${title} — 状態ノルムの推移`);
      const chart = window.plotlib.createChart(body, { xlim: [0, nT - 1], ylim: [-normMax * 0.05, normMax * 1.15], xlabel: "$k$" });
      for (const row of xNormHistP) chart.line(ks, row, { color: [0.7, 0.7, 0.7], lineWidth: 1 });
      chart.line(ks, meanCurve(xNormHistP), { color: "black", lineWidth: 2, label: "\\text{mean}" });
      chart.finish();
    }
    {
      const nU = UpsErrHistP[0].length;
      const ksU = Array.from({ length: nU }, (_, i) => i);
      const body = window.plotlib.makeCard(outputGrid, `${title} — Q関数（Upsilon）誤差`);
      const chart = window.plotlib.createChart(body, { xlim: [0, nU - 1], ylim: [-0.2, 1.2], xlabel: "$k$" });
      for (const row of UpsErrHistP) chart.line(ksU, row, { color: [0.7, 0.7, 0.7], lineWidth: 1 });
      chart.line(ksU, meanCurve(UpsErrHistP), { color: "black", lineWidth: 1.5, label: "\\text{mean}" });
      chart.finish();
    }
    if (iterGain > 1) {
      const nK = KErrHistP[0].length;
      const ksK = Array.from({ length: nK }, (_, i) => i);
      // Fixed on purpose (not data-fit like the panels above), but 12 was
      // much taller than this curve can ever get: K starts from the same
      // fixed K_ini every run, compared against a fixed K_opt, so the value
      // at k=0 -- and thus the curve's ceiling -- is always the same ~7.1
      // regardless of any parameter. A [0,12] axis left ~40% of the chart
      // permanently empty, visually flattening the part that *does* change
      // (the decay shape/rate) whenever a parameter changes. 8 fits the
      // actual ceiling with a little headroom instead.
      const body = window.plotlib.makeCard(outputGrid, `${title} — ゲイン誤差`);
      const chart = window.plotlib.createChart(body, { xlim: [0, nK - 1], ylim: [0, 8], xlabel: "$k$" });
      chart.line(ksK, meanCurve(KErrHistP), { color: "red", lineWidth: 2, label: "\\text{Gain error}" });
      chart.finish();
    }
  }

  // 9.3(a)/(b) are single-gain-update baselines contrasting exploration
  // noise level; 9.4 reuses 9.3(b)'s noise level to show the effect of
  // repeated gain updates -- iter_gain is only adjustable for 9.4 so that
  // distinction (single vs. repeated updates) stays intact.
  plotScenario("Figure 9.3(a)", params.sigma_a, 1);
  plotScenario("Figure 9.3(b)", params.sigma_b, 1);
  plotScenario("Figure 9.4", params.sigma_b, Math.max(1, Math.round(params.iter_gain)));
};
