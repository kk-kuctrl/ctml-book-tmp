"use strict";

// Port of code_distribute/python/Figure6_1.py
//
// control.dlqr(...) (solving the discrete algebraic Riccati equation for
// P_opt) is replaced by running the exact same Riccati fixed-point
// recursion the script already hand-implements for "value iteration" --
// just for many more iterations, until it converges. solve_discrete_lyapunov
// -> linalg.discreteLyapunov (already built). The one genuinely new piece is
// checking "is this closed loop stable" without a general eigenvalue
// solver, which linalg.spectralRadius (power iteration) handles.
window.figureLib = window.figureLib || {};
window.figureLib.figure6_1 = function (outputGrid, params) {
  const beta = params.beta;
  const nIter = Math.max(3, Math.round(params.n_iter));
  const L = window.linalg;

  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function outer(a, b) {
    return a.map((ai) => b.map((bj) => ai * bj));
  }

  // A, B, Q are free-text (see the params panel) so they can be typed in
  // directly instead of via sliders -- a 3x3 matrix + a 3-vector is a lot of
  // numbers to expose as 15 separate range inputs. Rows are newline/`;`
  // separated, entries within a row are comma/space separated. Anything
  // that doesn't parse to a well-formed 3x3 matrix (or 3-vector) falls back
  // to the textbook's original (A, B, Q) as a whole, rather than risking a
  // dimension mismatch from mixing a parsed and a default value.
  const DEFAULT_A_TEXT = "0.8,0.9,0.86\n0.3,0.25,1\n0.1,0.55,0.5";
  const DEFAULT_B_TEXT = "1,0,0";
  const DEFAULT_Q_TEXT = "1,0,0\n0,1,0\n0,0,1";
  function parseSquareMatrix(text, n) {
    const rows = String(text)
      .trim()
      .split(/[\n;]+/)
      .map((r) => r.trim())
      .filter((r) => r.length > 0);
    if (rows.length !== n) return null;
    const M = rows.map((r) =>
      r
        .split(/[,\s]+/)
        .filter((s) => s.length > 0)
        .map(Number)
    );
    if (M.some((row) => row.length !== n || row.some((v) => !Number.isFinite(v)))) return null;
    return M;
  }
  function parseVec(text, n) {
    const v = String(text)
      .trim()
      .split(/[,\s]+/)
      .filter((s) => s.length > 0)
      .map(Number);
    if (v.length !== n || v.some((x) => !Number.isFinite(x))) return null;
    return v;
  }

  const xDim = 3;
  let A = parseSquareMatrix(params.A_text, xDim);
  let Bvec = parseVec(params.B_text, xDim);
  let Q = parseSquareMatrix(params.Q_text, xDim);
  if (!A || !Bvec || !Q) {
    A = parseSquareMatrix(DEFAULT_A_TEXT, xDim);
    Bvec = parseVec(DEFAULT_B_TEXT, xDim);
    Q = parseSquareMatrix(DEFAULT_Q_TEXT, xDim);
  }
  const R = 1;
  const Svec = [0, 0, 0];

  // One step of the (discounted) Riccati recursion: given the current PI,
  // returns the optimal gain K (as a length-xDim row vector) and the
  // updated PI. Mirrors figure5_3.js's lqrControl, generalized with a
  // discount factor beta and a (here zero) state/control cross term S.
  function riccatiStep(PI) {
    const APIA = L.matMul(L.transpose(A), L.matMul(PI, A));
    const PIBu = L.matVec(PI, Bvec);
    const Rt = R + beta * dot(Bvec, PIBu);
    const St = Svec.map((s, i) => s + beta * L.matVec(L.transpose(A), PIBu)[i]);
    const Qt = L.add(Q, L.scale(APIA, beta));
    const K = St.map((v) => v / Rt);
    const PInext = L.sub(Qt, outer(St, K));
    return { K, PInext };
  }

  function frobeniusNorm(M) {
    let s = 0;
    for (const row of M) for (const v of row) s += v * v;
    return Math.sqrt(s);
  }

  // Ground truth P_opt: iterate the same recursion to convergence (stands
  // in for control.dlqr's direct DARE solve).
  let PIopt = L.zeros(xDim, xDim);
  for (let i = 0; i < 3000; i++) PIopt = riccatiStep(PIopt).PInext;
  const Kopt = riccatiStep(PIopt).K;

  // ---------------------------------------------------------------
  // Value iteration: find a stabilizing starting point, then track the
  // Frobenius-norm error to P_opt over nIter further iterations.
  // ---------------------------------------------------------------
  let PI = L.zeros(xDim, xDim);
  let K, PInext;
  for (let guard = 0; guard < 2000; guard++) {
    ({ K, PInext } = riccatiStep(PI));
    const Acl = L.sub(A, outer(Bvec, K));
    const rho = L.spectralRadius(L.scale(Acl, Math.sqrt(beta)), 300);
    if (rho < 1) break;
    PI = PInext;
  }
  const PIini = PI;
  const Kini = K;

  const errListVI = new Array(nIter);
  PI = PIini;
  for (let i = 0; i < nIter; i++) {
    errListVI[i] = frobeniusNorm(L.sub(PIopt, PI));
    PI = riccatiStep(PI).PInext;
  }

  // ---------------------------------------------------------------
  // Policy iteration: evaluate the current policy exactly (discrete
  // Lyapunov equation) instead of iterating the Bellman map.
  // ---------------------------------------------------------------
  function policyEvalQterm(Krow) {
    // M = [I; -K]  (4 x 3), QRSfull = [[Q,S],[S^T,R]] (4 x 4).
    // Returns M^T @ QRSfull @ M (3 x 3).
    const M = [...L.eye(xDim), Krow.map((v) => -v)]; // 4 rows x 3 cols
    const QRSfull = [
      [Q[0][0], Q[0][1], Q[0][2], Svec[0]],
      [Q[1][0], Q[1][1], Q[1][2], Svec[1]],
      [Q[2][0], Q[2][1], Q[2][2], Svec[2]],
      [Svec[0], Svec[1], Svec[2], R],
    ];
    return L.matMul(L.transpose(M), L.matMul(QRSfull, M));
  }

  const errListPI = new Array(nIter);
  let Kpi = Kini;
  let PIpi = PIini;
  for (let i = 0; i < nIter; i++) {
    errListPI[i] = frobeniusNorm(L.sub(PIopt, PIpi));
    const Acl = L.sub(A, outer(Bvec, Kpi));
    const Qterm = policyEvalQterm(Kpi);
    const PIQ = L.discreteLyapunov(L.scale(L.transpose(Acl), Math.sqrt(beta)), Qterm);

    const PIQBu = L.matVec(PIQ, Bvec);
    const Rt = R + beta * dot(Bvec, PIQBu);
    const St = Svec.map((s, i2) => s + beta * L.matVec(L.transpose(A), PIQBu)[i2]);
    const Qt = L.add(Q, L.scale(L.matMul(L.transpose(A), L.matMul(PIQ, A)), beta));
    Kpi = St.map((v) => v / Rt);
    PIpi = L.sub(Qt, outer(St, Kpi));
  }

  // ---------------------------------------------------------------
  // Plot: convergence of both methods toward P_opt, log-scale.
  // ---------------------------------------------------------------
  const body = window.plotlib.makeCard(outputGrid, "Figure 6.1 — 価値反復と方策反復の収束比較");
  const chart = window.plotlib.createChart(body, {
    xlim: [0, nIter - 1],
    ylim: [1e-9, 1000],
    yscale: "log",
    xlabel: "\\text{Iterations}",
  });
  const iters = Array.from({ length: nIter }, (_, i) => i);
  chart.line(iters, errListVI, { color: "blue", lineWidth: 1.5, label: "\\text{Value Iteration}" });
  chart.line(iters, errListPI, { color: "red", lineWidth: 1.5, label: "\\text{Policy Iteration}" });
  chart.finish();

  const gainLine = document.createElement("p");
  gainLine.className = "hint";
  gainLine.style.margin = "8px 0 0";
  gainLine.textContent = `最適ゲイン K_opt = (${Kopt.map((v) => v.toFixed(4)).join(", ")})`;
  body.appendChild(gainLine);
};
