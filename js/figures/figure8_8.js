"use strict";

// Port of code_distribute/python/Figure8_8.py
// solve_discrete_lyapunov -> linalg.discreteLyapunov (already built for
// Figure 3.4/5.3); the only new piece is linalg.spectralRadius, used to
// replicate `A / max(|eigvals(A)|) * 0.95` without a general eigensolver.
window.figureLib = window.figureLib || {};
window.figureLib.figure8_8 = function (outputGrid, params) {
  const kBar = Math.max(20, Math.round(params.k_bar));
  const nSample = Math.max(1, Math.round(params.n_sample));
  const yDim = 2; // output dimension is fixed: the 3D plot needs exactly (y)_1, (y)_2
  const L = window.linalg;

  // A is free-text (typed in, like Figure 6.1's A/B/Q) and any square size --
  // its own row count IS the state dimension, so typing a 3x3 or an 8x8
  // both just work. B stays randomized per Run (sized to match) so there's
  // still some run-to-run variation. C must still produce exactly 2 outputs
  // (the 3D plot needs (y)_1, (y)_2) but its column count follows A's size.
  const DEFAULT_A_TEXT =
    "0.42,0.72,0.00,0.30,0.15\n0.09,0.19,0.35,0.40,0.54\n0.42,0.69,0.20,0.88,0.03\n0.67,0.42,0.56,0.14,0.20\n0.80,0.97,0.31,0.69,0.88";
  function parseSquareMatrixAny(text) {
    const rows = String(text)
      .trim()
      .split(/[\n;]+/)
      .map((r) => r.trim())
      .filter((r) => r.length > 0);
    const n = rows.length;
    if (n === 0) return null;
    const M = rows.map((r) =>
      r
        .split(/[,\s]+/)
        .filter((s) => s.length > 0)
        .map(Number)
    );
    if (M.some((row) => row.length !== n || row.some((v) => !Number.isFinite(v)))) return null;
    return M;
  }
  function parseMatrixFixed(text, rows, cols) {
    const rowsArr = String(text)
      .trim()
      .split(/[\n;]+/)
      .map((r) => r.trim())
      .filter((r) => r.length > 0);
    if (rowsArr.length !== rows) return null;
    const M = rowsArr.map((r) =>
      r
        .split(/[,\s]+/)
        .filter((s) => s.length > 0)
        .map(Number)
    );
    if (M.some((row) => row.length !== cols || row.some((v) => !Number.isFinite(v)))) return null;
    return M;
  }
  // Falls back to "pick out the first two state components" (Ups the
  // Python's implicit behavior) at whatever size xDim turns out to be --
  // there's no fixed-text default for C since its width depends on A.
  function defaultC(n) {
    const C = L.zeros(yDim, n);
    for (let i = 0; i < yDim && i < n; i++) C[i][i] = 1;
    return C;
  }

  let Araw = parseSquareMatrixAny(params.A_text);
  if (!Araw) Araw = parseSquareMatrixAny(DEFAULT_A_TEXT);
  const xDim = Araw.length;

  let C = parseMatrixFixed(params.C_text, yDim, xDim);
  if (!C) C = defaultC(xDim);

  // Only rescaled (to spectral radius 0.95, matching the Python's
  // A/max(|eig|)*0.95) when the typed-in A is actually unstable -- an
  // already-stable A is used exactly as typed, instead of always forcing
  // every A to the same spectral radius regardless of what was entered.
  const rho = L.spectralRadius(Araw);
  const A = rho >= 1 ? L.scale(Araw, 0.95 / rho) : Araw;
  const B = Array.from({ length: xDim }, () => [Math.random()]);
  const BBt = L.matMul(B, L.transpose(B));
  const X = L.discreteLyapunov(A, BBt);
  const Bvec = B.map((row) => row[0]);

  const kHalf = Math.floor(kBar / 2);
  const body = window.plotlib.makeCard(outputGrid, `Figure 8.8 — ${xDim}次元系のサンプル軌道（最初の2成分を3D表示）`);
  const chart3d = window.plotlib.createChart3D(body, {
    width: 960 * 0.9,
    height: 420 * 0.9,
    contentScale: 0.8,
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
    let y = L.matVec(C, x);
    const y1 = [y[0]];
    const y2 = [y[1]];
    for (let t = 0; t < kBar; t++) {
      const noise = window.rnd.randn();
      x = L.matVec(A, x).map((v, i) => v + Bvec[i] * noise);
      y = L.matVec(C, x);
      y1.push(y[0]);
      y2.push(y[1]);
    }
    chart3d.addLine(kVals, y1, y2, palette[s % palette.length]);
  }
  chart3d.finish();
};
