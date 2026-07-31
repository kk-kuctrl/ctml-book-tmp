"use strict";

// Port of code_distribute/python/Figure5_3.py
//
// Everything here that looked like it needed the `control` package reduces
// to code already built for Figure 3.4 / Figure 4.2:
//   - create_double_integrator() is a ZOH discretization -> linalg.c2dZOH
//   - create_noise_model() is IDENTICAL to Figure 3.4's version (Tustin +
//     discrete-Lyapunov normalization) -> linalg.c2dTustin + discreteLyapunov
// The LQR Riccati recursion and the Kalman filter/LQG simulation never
// touched `control` at all in the original -- just numpy matrix algebra,
// simplified here using the fact that this system always has a single
// control input and a single measured output (so several matrix products
// collapse to plain vectors/scalars).
window.figureLib = window.figureLib || {};
window.figureLib.figure5_3 = function (outputGrid, params) {
  const kBar = Math.max(20, Math.round(params.k_bar));
  const Ts = 0.1;
  const aNoise = params.a;

  const L = window.linalg;

  function outer(a, b) {
    return a.map((ai) => b.map((bj) => ai * bj));
  }

  // ---------------------------------------------------------------
  // Noise model: F_c = 1/(s+a) -> Tustin discretize -> normalize
  // (identical construction to figure3_4.js's create_noise_model; the
  // book's original fixes a=0.3, here it's an adjustable parameter).
  // ---------------------------------------------------------------
  function createNoiseModel(Ts, a) {
    const Ac = [[-a]],
      Bc = [[1]],
      Cc = [[1]],
      Dc = [[0]];
    const { Ad, Bd, Cd } = L.c2dTustin(Ac, Bc, Cc, Dc, Ts);
    const Q = L.matMul(Bd, L.transpose(Bd));
    const P = L.discreteLyapunov(Ad, Q);
    const varAmp = L.matMul(Cd, L.matMul(P, L.transpose(Cd)))[0][0];
    const Bnorm = L.scale(Bd, 1 / Math.sqrt(varAmp));
    return { Aw: Ad[0][0], Bw: Bnorm[0][0], Cw: Cd[0][0] };
  }

  // Double integrator P = 1/s^2, ZOH discretized.
  function createDoubleIntegrator(Ts) {
    const Ac = [
      [0, 1],
      [0, 0],
    ];
    const Bc = [[0], [1]];
    const Cc = [1, 0]; // unchanged by ZOH
    const { Ad, Bd } = L.c2dZOH(Ac, Bc, Ts);
    return { Ar: Ad, Br: Bd.map((r) => r[0]), Cr: Cc };
  }

  const { Aw, Bw, Cw } = createNoiseModel(Ts, aNoise);
  const { Ar, Br, Cr } = createDoubleIntegrator(Ts);
  const nr = 2,
    nw = 1,
    xDim = 3;

  // Augmented system: state = [position, velocity, noise-filter state].
  const A = [
    [Ar[0][0], Ar[0][1], Br[0] * Cw],
    [Ar[1][0], Ar[1][1], Br[1] * Cw],
    [0, 0, Aw],
  ];
  const Bu = [Br[0], Br[1], 0]; // control input matrix, flattened (single input)
  const Bv = [0, 0, Bw]; // disturbance input matrix, flattened
  const Cvec = [Cr[0], Cr[1], 0]; // output matrix, flattened (single output)

  const Q = [
    [1, 0, 0],
    [0, 10, 0],
    [0, 0, 0],
  ];
  const R = 1e-4;
  const Qf = Q;
  const Rv = 1;
  const Rw = 1e-4;

  // ---------------------------------------------------------------
  // Finite-horizon discrete-time LQR (single-input Riccati recursion).
  // A: n x n, Bu: n-vector, Q/Qf: n x n, R: scalar.
  // Returns K[k] (n-vector, the row gain) and Sigma[k] (n x n), k=0..kBar.
  // ---------------------------------------------------------------
  function lqrControl(A, Bu, Q, R, Qf, kBar) {
    const n = A.length;
    const Sigma = new Array(kBar + 1);
    const K = new Array(kBar);
    Sigma[kBar] = Qf;
    const At = L.transpose(A);
    for (let k = kBar - 1; k >= 0; k--) {
      const SBu = L.matVec(Sigma[k + 1], Bu); // Sigma[k+1] @ Bu
      const Rtilde = R + dot(Bu, SBu);
      const Stilde = L.matVec(At, SBu); // A^T Sigma[k+1] Bu  (n-vector)
      const Qtilde = L.add(Q, L.matMul(At, L.matMul(Sigma[k + 1], A)));
      const Kk = Stilde.map((v) => v / Rtilde);
      K[k] = Kk;
      Sigma[k] = L.sub(Qtilde, outer(Stilde, Kk));
    }
    return { K, Sigma };
  }

  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  // ---------------------------------------------------------------
  // LQG / LQR simulation (single-input, single-output Kalman filter).
  // mode: 'lqr' | 'lqg_pred' | 'lqg_kalman'
  // ---------------------------------------------------------------
  function simulateLQ(A, Bu, Bv, C, mu, Sigma0, K, kBar, x0, mode, Rw, Rv, v, w) {
    const n = A.length;
    const xTrue = [x0.slice()];
    const xHat = [mu.slice()];
    const sigmasDiag = [L.diag(Sigma0)];
    const xCheck = [];
    const yArr = [];

    let Sigma = Sigma0.map((row) => row.slice());

    for (let k = 0; k < kBar; k++) {
      let uk;
      if (mode === "lqr") {
        uk = -dot(K[k], xTrue[k]);
      } else {
        if (mode === "lqg_pred") uk = -dot(K[k], xHat[k]);

        const yk = dot(C, xTrue[k]) + w[k];
        yArr.push(yk);

        const Mtilde = dot(C, L.matVec(Sigma, C)) + Rw;
        const Lcheck = L.matVec(Sigma, C);
        const Hcheck = Lcheck.map((v2) => v2 / Mtilde);

        const innov = yk - dot(C, xHat[k]);
        const xCheckK = xHat[k].map((v2, i) => v2 + Hcheck[i] * innov);
        xCheck.push(xCheckK);
        const SigmaCheck = L.sub(Sigma, outer(Hcheck, Lcheck));

        if (mode === "lqg_kalman") uk = -dot(K[k], xCheckK);

        const xHatNext = A.map((row, i) => dot(row, xCheckK) + Bu[i] * uk);
        xHat.push(xHatNext);
        Sigma = L.add(L.matMul(A, L.matMul(SigmaCheck, L.transpose(A))), L.scale(outer(Bv, Bv), Rv));
      }

      const xTrueNext = A.map((row, i) => dot(row, xTrue[k]) + Bu[i] * uk + Bv[i] * v[k]);
      xTrue.push(xTrueNext);
      sigmasDiag.push(L.diag(Sigma));
    }

    return { xTrue, xHat, y: yArr, sigmasDiag };
  }

  // ---------------------------------------------------------------
  // Set up the scenario and run all four simulations.
  // ---------------------------------------------------------------
  const mu = new Array(xDim).fill(0);
  const Sigma0 = L.eye(xDim);
  const x0 = window.rnd.mvnSample(mu, Sigma0).map((v) => v * 0.1);

  const processNoise = window.rnd.randnVec(kBar).map((v) => Math.sqrt(Rv) * v);
  const obsNoise = window.rnd.randnVec(kBar).map((v) => Math.sqrt(Rw) * v);

  const { K } = lqrControl(A, Bu, Q, R, Qf, kBar);
  const Q2 = [
    [1, 0],
    [0, 10],
  ];
  const { K: Klqr } = lqrControl(Ar, Br, Q2, R, Q2, kBar);

  const lqg = simulateLQ(A, Bu, Bv, Cvec, mu, Sigma0, K, kBar, x0, "lqg_pred", Rw, Rv, processNoise, obsNoise);
  const lqgK = simulateLQ(A, Bu, Bv, Cvec, mu, Sigma0, K, kBar, x0, "lqg_kalman", Rw, Rv, processNoise, obsNoise);

  const mu2 = mu.slice(0, 2);
  const Sigma02 = [
    [1, 0],
    [0, 1],
  ];
  const x02 = x0.slice(0, 2);
  const lqrWhite = simulateLQ(Ar, Br, Br, Cr, mu2, Sigma02, Klqr, kBar, x02, "lqr", Rw, Rv, processNoise, obsNoise);

  const KlqrAug = Klqr.map((row) => row.concat(new Array(xDim - 2).fill(0)));
  const lqrColored = simulateLQ(A, Bu, Bv, Cvec, mu, Sigma0, KlqrAug, kBar, x0, "lqr", Rw, Rv, processNoise, obsNoise);

  // ---------------------------------------------------------------
  // Plotting
  // ---------------------------------------------------------------
  const col = (arr, j) => arr.map((row) => row[j]);
  const ks = Array.from({ length: kBar + 1 }, (_, k) => k);

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.3(a) — LQR/LQGの位置比較");
    const chart = window.plotlib.createChart(body, { xlim: [0, kBar], ylim: [-0.5, 0.5], xlabel: "$k$" });
    chart.line(ks, col(lqrWhite.xTrue, 0), { color: [0.6, 0.6, 0.6], dash: [6, 4], lineWidth: 1.5, label: "\\text{LQR (white noise)}" });
    chart.line(ks, col(lqrColored.xTrue, 0), { color: "black", lineWidth: 1.5, label: "\\text{LQR}" });
    chart.line(ks, col(lqg.xTrue, 0), { color: "blue", lineWidth: 1.5, label: "\\text{LQG}" });
    chart.line(ks, col(lqgK.xTrue, 0), { color: "red", lineWidth: 1.5, label: "\\text{LQG Kalman}" });
    chart.finish();
  }

  const xlimBC = [0, Math.min(100, kBar)];
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.3(b) — 位置の推定");
    const chart = window.plotlib.createChart(body, { xlim: xlimBC, ylim: [-0.2, 0.2], xlabel: "$k$" });
    const sd = lqg.sigmasDiag.map((d) => Math.sqrt(Math.max(0, d[0])));
    const hat0 = col(lqg.xHat, 0);
    chart.fillBetween(
      ks,
      hat0.map((v, i) => v - sd[i]),
      hat0.map((v, i) => v + sd[i]),
      { color: "blue", alpha: 0.2 }
    );
    chart.line(ks, col(lqg.xTrue, 0), { color: "blue", lineWidth: 1.5, label: "\\text{Position}" });
    chart.line(ks.slice(0, kBar), lqg.y, { color: "red", lineWidth: 1, label: "\\text{Measurements}" });
    chart.line(ks, hat0, { color: "green", dash: [5, 4], lineWidth: 1.5, label: "\\text{Estimate}" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.3(c) — 速度の推定");
    const chart = window.plotlib.createChart(body, { xlim: xlimBC, ylim: [-0.5, 0.5], xlabel: "$k$" });
    const sd = lqg.sigmasDiag.map((d) => Math.sqrt(Math.max(0, d[1])));
    const hat1 = col(lqg.xHat, 1);
    chart.fillBetween(
      ks,
      hat1.map((v, i) => v - sd[i]),
      hat1.map((v, i) => v + sd[i]),
      { color: "blue", alpha: 0.2 }
    );
    chart.line(ks, col(lqg.xTrue, 1), { color: "blue", lineWidth: 1.5, label: "\\text{Velocity}" });
    chart.line(ks, hat1, { color: "red", dash: [5, 4], lineWidth: 1.5, label: "\\text{Estimate}" });
    chart.finish();
  }
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 5.3(d) — 色付き雑音成分の推定");
    const chart = window.plotlib.createChart(body, { xlim: xlimBC, ylim: [-2, 2], xlabel: "$k$" });
    const sd = lqg.sigmasDiag.map((d) => Math.sqrt(Math.max(0, d[2])));
    const hat2 = col(lqg.xHat, 2);
    chart.fillBetween(
      ks,
      hat2.map((v, i) => v - sd[i]),
      hat2.map((v, i) => v + sd[i]),
      { color: "blue", alpha: 0.2 }
    );
    chart.line(ks, col(lqg.xTrue, 2), { color: "blue", lineWidth: 1.5, label: "\\text{Colored noise}" });
    chart.line(ks, hat2, { color: "red", dash: [5, 4], lineWidth: 1.5, label: "\\text{Estimate}" });
    chart.finish();
  }
};
