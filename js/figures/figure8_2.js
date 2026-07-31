"use strict";

// Port of code_distribute/python/Figure8_2.py
window.figureLib = window.figureLib || {};
window.figureLib.figure8_2 = function (outputGrid, params) {
  const kBarA = Math.max(3, Math.round(params.k_bar_a));
  const kBarB = Math.max(3, Math.round(params.k_bar_b));
  const kBarCdef = Math.max(3, Math.round(params.k_bar_cdef));

  // ---------------------------------------------------------------
  // Shared helpers (ports of initialization() / sysid_module() in the
  // Python source).
  // ---------------------------------------------------------------

  // ARX system parameters, shared by all sub-figures.
  function initialization() {
    const a = [1.2, -0.47, 0.06]; // coeffs of (z-0.5)(z-0.4)(z-0.3)
    const b = [1.0, 2.0]; // b2, b1
    const q0 = [0, 0, 0, 1, 1]; // initial states/inputs, note u0=u1=1
    const p0 = [0, 0, 0, 0, 0]; // initial parameter estimate
    const Sigma0 = window.linalg.eye(5).map((row) => row.map((v) => v * 1e4));
    return { a, b, q0, p0, Sigma0 };
  }

  function dot(x, y) {
    let s = 0;
    for (let i = 0; i < x.length; i++) s += x[i] * y[i];
    return s;
  }

  function sumSq(x, y) {
    let s = 0;
    for (let i = 0; i < x.length; i++) s += (x[i] - y[i]) ** 2;
    return s;
  }

  function scaleMat(A, s) {
    return A.map((row) => row.map((v) => v / s));
  }

  // Recursive least-squares system identification.
  // p_star: (k_dat x p_dim) array-of-arrays, the TRUE parameter sequence.
  // a_dim: number of "a" parameters (p_dim - a_dim = number of "b" parameters).
  // q0, u, v, p0, Sigma0, alpha: as in the Python sysid_module().
  // Returns { aErr, bErr, trSigma, pEst, y } (plain arrays / arrays-of-arrays).
  function sysidModule(pStar, aDim, q0, u, v, p0, Sigma0, alpha) {
    const kDat = pStar.length;
    const pDim = pStar[0].length;

    const y = new Array(kDat - 1);
    const q = new Array(kDat);
    q[0] = q0.slice();

    // Open-loop simulation to generate y data.
    // q_k = [y_{k-1}, y_{k-2}, y_{k-3}, u_{k-1}, u_{k-2}]  (for a_dim=3, p_dim=5)
    // Matches numpy: q[k+1,:] = hstack([y[k], q[k,:a_dim-1], u[k], q[k,a_dim:p_dim-1]])
    for (let k = 0; k < kDat - 1; k++) {
      y[k] = dot(pStar[k], q[k]) + v[k];
      q[k + 1] = [y[k]]
        .concat(q[k].slice(0, aDim - 1))
        .concat([u[k]])
        .concat(q[k].slice(aDim, pDim - 1));
    }

    const pEst = new Array(kDat);
    pEst[0] = p0.slice();
    const aErr = new Array(kDat);
    const bErr = new Array(kDat);
    const trSigma = new Array(kDat);

    let Sigma = Sigma0.map((row) => row.slice());

    aErr[0] = sumSq(pStar[0].slice(0, aDim), p0.slice(0, aDim));
    bErr[0] = sumSq(pStar[0].slice(aDim), p0.slice(aDim));
    trSigma[0] = window.linalg.trace(Sigma);

    for (let k = 0; k < kDat - 1; k++) {
      const qk = q[k];
      const SigmaQk = window.linalg.matVec(Sigma, qk);
      const denom = alpha + dot(qk, SigmaQk);
      const H = SigmaQk.map((v2) => v2 / denom);
      const innovation = y[k] - dot(pEst[k], qk);
      pEst[k + 1] = pEst[k].map((v2, i) => v2 + H[i] * innovation);

      const outerHq = H.map((hi) => qk.map((qj) => hi * qj));
      const prod = window.linalg.matMul(outerHq, Sigma);
      const newSigma = Sigma.map((row, i) => row.map((val, j) => (val - prod[i][j]) / alpha));
      Sigma = newSigma.map((row, i) => row.map((val, j) => 0.5 * (val + newSigma[j][i])));

      aErr[k + 1] = sumSq(pStar[k + 1].slice(0, aDim), pEst[k + 1].slice(0, aDim));
      bErr[k + 1] = sumSq(pStar[k + 1].slice(aDim), pEst[k + 1].slice(aDim));
      trSigma[k + 1] = window.linalg.trace(Sigma);
    }

    return { aErr, bErr, trSigma, pEst, y };
  }

  // ---------------------------------------------------------------
  // Figure 8.2(a): open-loop identification, two noise/input settings.
  // ---------------------------------------------------------------
  {
    const { a, b, q0, p0, Sigma0 } = initialization();
    const aDim = a.length;
    const pStarRow = a.concat(b);
    const pStar = Array.from({ length: kBarA }, () => pStarRow.slice());

    const alpha = 1;

    // u_k = 1, v_k ~ N(0,1)
    let sigmaV = 1;
    let v = window.rnd.randnVec(kBarA).map((z) => z * sigmaV);
    let u = new Array(kBarA).fill(1);
    const res1 = sysidModule(pStar, aDim, q0, u, v, p0, scaleMat(Sigma0, sigmaV), alpha);

    // u_k ~ N(1,1), v_k ~ N(0,0.01)
    sigmaV = 0.1;
    v = window.rnd.randnVec(kBarA).map((z) => z * sigmaV);
    u = window.rnd.randnVec(kBarA).map((z) => 1 + z);
    const res2 = sysidModule(pStar, aDim, q0, u, v, p0, scaleMat(Sigma0, sigmaV), alpha);

    const body = window.plotlib.makeCard(outputGrid, "Figure 8.2(a) — RLS output trajectories");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBarA],
      ylim: (function () {
        let lo = Infinity,
          hi = -Infinity;
        for (const arr of [res1.y, res2.y]) {
          for (const val of arr) {
            if (val < lo) lo = val;
            if (val > hi) hi = val;
          }
        }
        const pad = (hi - lo) * 0.1 || 1;
        return [lo - pad, hi + pad];
      })(),
      xlabel: "$k$",
      ylabel: "$y_k$",
    });
    const ks1 = Array.from({ length: res1.y.length }, (_, i) => i);
    const ks2 = Array.from({ length: res2.y.length }, (_, i) => i);
    chart.line(ks1, res1.y, { lineWidth: 0.5, color: "blue", label: "u_k=1,v_k\\sim {\\mathcal N}(0,1)" });
    chart.line(ks2, res2.y, {
      lineWidth: 0.5,
      color: "orange",
      label: "u_k \\sim {\\mathcal N}(1,1), v_k\\sim {\\mathcal N}(0,0.01)",
    });
    chart.finish();
  }

  // ---------------------------------------------------------------
  // Figure 8.2(b): time-varying true parameter, three forgetting factors.
  // ---------------------------------------------------------------
  {
    const { a, b, q0, p0, Sigma0 } = initialization();
    const aDim = a.length;
    const aChanged = [1.0, -0.47, 0.06];

    const pStar = new Array(kBarB);
    for (let r = 0; r < kBarB; r++) {
      pStar[r] = (r % 2000 > 1000 ? aChanged : a).concat(b);
    }

    const sigmaV = 0.1;
    const v = window.rnd.randnVec(kBarB).map((z) => z * sigmaV);
    const u = window.rnd.randnVec(kBarB);
    const Sigma0Div = scaleMat(Sigma0, sigmaV);

    const resAlpha1 = sysidModule(pStar, aDim, q0, u, v, p0, Sigma0Div, 1.0);
    const resAlpha995 = sysidModule(pStar, aDim, q0, u, v, p0, Sigma0Div, 0.995);
    const resAlpha8 = sysidModule(pStar, aDim, q0, u, v, p0, Sigma0Div, 0.8);

    const body = window.plotlib.makeCard(outputGrid, "Figure 8.2(b) — Tracking a time-varying parameter");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBarB],
      ylim: [0.8, 1.4],
      xlabel: "$k$",
      ylabel: "{\\rm a}_1",
    });
    const ks = Array.from({ length: kBarB }, (_, i) => i);
    chart.line(ks, resAlpha1.pEst.map((row) => row[0]), { lineWidth: 0.5, color: "blue", label: "\\alpha=1" });
    chart.line(ks, resAlpha995.pEst.map((row) => row[0]), { lineWidth: 0.5, color: "orange", label: "\\alpha=0.995" });
    chart.line(ks, resAlpha8.pEst.map((row) => row[0]), { lineWidth: 0.5, color: "green", label: "\\alpha=0.8" });
    chart.line(ks, pStar.map((row) => row[0]), { lineWidth: 1.5, color: "black", label: "True" });
    chart.finish();
  }

  // ---------------------------------------------------------------
  // Figure 8.2(c)-(f): convergence of estimation error / covariance trace
  // on log-log axes, under four (u, v) noise combinations.
  // ---------------------------------------------------------------
  const cdefSpecs = [
    { label: "c", title: "$u_k=1$, $v_k\\sim {\\mathcal N}(0,1)$" },
    { label: "d", title: "$u_k\\sim {\\mathcal N}(1,1)$, $v_k\\sim {\\mathcal N}(0,1)$" },
    { label: "e", title: "$u_k=1$, $v_k\\sim {\\mathcal N}(0,0.01)$" },
    { label: "f", title: "$u_k\\sim {\\mathcal N}(1,1)$, $v_k\\sim {\\mathcal N}(0,0.01)$" },
  ];

  for (const spec of cdefSpecs) {
    const label = spec.label;
    const { a, b, q0, p0, Sigma0 } = initialization();
    const aDim = a.length;
    const pStarRow = a.concat(b);
    const pStar = Array.from({ length: kBarCdef }, () => pStarRow.slice());

    const sigmaV = label === "c" || label === "d" ? 1 : 0.1;
    let u, v;
    if (label === "c" || label === "e") {
      u = new Array(kBarCdef).fill(1);
      v = window.rnd.randnVec(kBarCdef).map((z) => z * sigmaV);
    } else {
      u = window.rnd.randnVec(kBarCdef).map((z) => 1 + z);
      v = window.rnd.randnVec(kBarCdef).map((z) => z * sigmaV);
    }

    const alpha = 1.0;
    const result = sysidModule(pStar, aDim, q0, u, v, p0, scaleMat(Sigma0, sigmaV), alpha);

    const body = window.plotlib.makeCard(outputGrid, `Figure 8.2(${label}) — Convergence (${spec.title})`);
    const chart = window.plotlib.createChart(body, {
      xscale: "log",
      yscale: "log",
      xlim: [1, kBarCdef],
      ylim: [1e-10, 1e5],
      xlabel: "$k$",
    });
    const ks = Array.from({ length: kBarCdef }, (_, i) => i + 1);
    chart.line(ks, result.aErr, { lineWidth: 1, color: "blue", label: "\\|\\check{\\rm p}^a - {\\rm a}^*\\|^2" });
    chart.line(ks, result.bErr, { lineWidth: 1, color: "orange", label: "\\|\\check{\\rm p}^b - {\\rm b}^*\\|^2" });
    chart.line(
      ks,
      result.trSigma.map((t) => t * sigmaV),
      { lineWidth: 1, color: "green", label: "{\\rm Trace}(\\check\\Sigma)" }
    );
    chart.finish();
  }
};
