"use strict";

// Port of code_distribute/python/Figure8_8.py
// solve_discrete_lyapunov -> linalg.discreteLyapunov (already built for
// Figure 3.4/5.3); the only new piece is linalg.spectralRadius, used to
// replicate `A / max(|eigvals(A)|) * 0.95` without a general eigensolver.
window.figureLib = window.figureLib || {};
window.figureLib.figure8_8 = function (outputGrid, params) {
  const kBar = Math.max(20, Math.round(params.k_bar));
  const nSample = Math.max(1, Math.round(params.n_sample));
  const xDim = 5;
  const L = window.linalg;

  let A = Array.from({ length: xDim }, () => Array.from({ length: xDim }, () => Math.random()));
  const rho = L.spectralRadius(A);
  A = L.scale(A, 0.95 / rho);
  const B = Array.from({ length: xDim }, () => [Math.random()]);
  const BBt = L.matMul(B, L.transpose(B));
  const X = L.discreteLyapunov(A, BBt);
  const Bvec = B.map((row) => row[0]);

  const kHalf = Math.floor(kBar / 2);
  const body = window.plotlib.makeCard(outputGrid, "Figure 8.8 — 5次元系のサンプル軌道（最初の2成分を3D表示）");
  const chart3d = window.plotlib.createChart3D(body, {
    width: 960,
    stretchToFill: true,
    elev: 10,
    azim: 3,
    xlabel: "$k$",
    ylabel: "$({\\rm y})_1$",
    zlabel: "$({\\rm y})_2$",
    // Fixed axis box (matches the Python's ax.set_xlim/ylim/zlim) instead
    // of auto-fitting to this run's random trajectories, so the framing
    // doesn't rescale (and isn't comparable) from run to run.
    xlim: [-kHalf, kHalf],
    ylim: [-3, 3],
    zlim: [-3, 3],
  });
  const kVals = Array.from({ length: kBar + 1 }, (_, i) => i - kHalf);
  const palette = ["blue", "orange", "green", "red", "purple"];

  for (let s = 0; s < nSample; s++) {
    let x = window.rnd.mvnSample(new Array(xDim).fill(0), X);
    const y1 = [x[0]];
    const y2 = [x[1]];
    for (let t = 0; t < kBar; t++) {
      const noise = window.rnd.randn();
      x = L.matVec(A, x).map((v, i) => v + Bvec[i] * noise);
      y1.push(x[0]);
      y2.push(x[1]);
    }
    chart3d.addLine(kVals, y1, y2, palette[s % palette.length]);
  }
  chart3d.finish();
};
