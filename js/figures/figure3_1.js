"use strict";
// Port of code_distribute/python/Figure3_1.py
window.figureLib = window.figureLib || {};
window.figureLib.figure3_1 = function (outputGrid, params) {
  const kBarA = Math.max(1, Math.round(params.k_bar_a));
  const kBarB = Math.max(2, Math.round(params.k_bar_b));
  const nSampleB = Math.max(1, Math.round(params.n_sample_b));

  // ---- Figure 3.1(a): distribution propagation via a discretized transition matrix ----
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.1(a) — 分布の時間発展（3次元表示）");

    // NOTE: the original Python script uses n_tmp=1000 (a 2001x2001 dense
    // transition matrix). For an instant, in-browser demo we use a coarser
    // grid, n_tmp=200 (401x401), which keeps the shape of the propagated
    // distribution while being fast to build and multiply.
    const nTmp = 200;
    const nX = 2 * nTmp + 1;
    const x = Array.from({ length: nX }, (_, i) => -1 + (2 * i) / (nX - 1));
    const dx = x[1] - x[0];

    // Build transition probability matrix P (uniform noise over an interval).
    const P = window.linalg.zeros(nX, nX);
    for (let i = 0; i < nX; i++) {
      const xi = x[i];
      const fb = xi + 0.1 * (xi - xi ** 3);
      const c = 0.5 * (1 - Math.abs(xi));
      const upper = fb + c;
      const lower = fb - c;
      let idxU = Math.floor(upper * nTmp + 0.5) + nTmp;
      let idxL = Math.ceil(lower * nTmp - 0.5) + nTmp;
      // Defensive clamp against floating-point edge effects at the boundary.
      idxU = Math.min(nX - 1, Math.max(0, idxU));
      idxL = Math.min(nX - 1, Math.max(0, idxL));
      if (idxL > idxU) [idxL, idxU] = [idxU, idxL];
      const w = 1 / (idxU - idxL + 1);
      for (let j = idxL; j <= idxU; j++) P[j][i] += w;
    }

    // Initial distribution: uniform on [-0.5, 0.5], normalized so sum*dx = 1.
    let init = x.map((xi) => (xi >= -0.5 && xi <= 0.5 ? 1.0 : 0.0));
    const initSum = init.reduce((s, v) => s + v, 0) * dx;
    init = init.map((v) => v / initSum);

    // Propagate: phi[:, k+1] = P @ phi[:, k]. Implemented directly (rather
    // than via window.linalg.matVec) since this loop is on the hot path for
    // a 401x401 matrix repeated k_bar_a times.
    const phi = [init.slice()];
    for (let k = 0; k < kBarA; k++) {
      const prev = phi[k];
      const next = new Array(nX).fill(0);
      for (let j = 0; j < nX; j++) {
        const row = P[j];
        let s = 0;
        for (let i = 0; i < nX; i++) s += row[i] * prev[i];
        next[j] = s;
      }
      phi.push(next);
    }

    const chart3d = window.plotlib.createChart3D(body, {
      elev: 26,
      azim: -107,
      xlabel: "$k$",
      ylabel: "${\\rm x}$",
      zlabel: "$\\varphi_{x_k}$",
    });
    const palette = ["blue", "orange", "green", "red", "purple"];
    for (let k = 0; k <= kBarA; k++) {
      const ks = new Array(nX).fill(k);
      chart3d.addLine(ks, x, phi[k], palette[k % palette.length]);
    }
    chart3d.finish();
  }

  // ---- Figure 3.1(b): sample trajectories under uniform noise ----
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.1(b) — サンプル軌道");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBarB - 1],
      ylim: [-1, 1],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });
    const ks = Array.from({ length: kBarB }, (_, k) => k);
    for (let s = 0; s < nSampleB; s++) {
      const xArr = new Array(kBarB).fill(0);
      xArr[0] = Math.random() - 0.5;
      for (let k = 0; k < kBarB - 1; k++) {
        const noise = Math.random() - 0.5;
        xArr[k + 1] = xArr[k] + 0.1 * (xArr[k] - xArr[k] ** 3) + noise * (1 - Math.abs(xArr[k]));
      }
      chart.line(ks, xArr, { lineWidth: 2 });
    }
    chart.finish();
  }
};
