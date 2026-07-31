"use strict";

// Minimal dense linear algebra helpers -- everything these figures need
// (small-to-medium symmetric systems, Cholesky for GP sampling) and nothing
// more. Matrices are plain arrays of arrays (row-major).

function zeros(n, m) {
  const rows = m === undefined ? 1 : n;
  const cols = m === undefined ? n : m;
  return Array.from({ length: rows }, () => new Array(cols).fill(0));
}

function eye(n) {
  const I = zeros(n, n);
  for (let i = 0; i < n; i++) I[i][i] = 1;
  return I;
}

function matMul(A, B) {
  const n = A.length,
    k = B.length,
    m = B[0].length;
  const C = zeros(n, m);
  for (let i = 0; i < n; i++) {
    for (let p = 0; p < k; p++) {
      const a = A[i][p];
      if (a === 0) continue;
      for (let j = 0; j < m; j++) C[i][j] += a * B[p][j];
    }
  }
  return C;
}

function matVec(A, v) {
  return A.map((row) => row.reduce((s, a, j) => s + a * v[j], 0));
}

function transpose(A) {
  const n = A.length,
    m = A[0].length;
  const T = zeros(m, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) T[j][i] = A[i][j];
  return T;
}

function addDiag(A, eps) {
  const B = A.map((row) => row.slice());
  for (let i = 0; i < B.length; i++) B[i][i] += eps;
  return B;
}

// Solve A x = b via Gaussian elimination with partial pivoting.
// A: n x n, b: n-vector or n x k matrix (solves for each column).
function solve(A, b) {
  const n = A.length;
  const bIsVector = !Array.isArray(b[0]);
  const B = bIsVector ? b.map((v) => [v]) : b.map((row) => row.slice());
  const k = B[0].length;
  const M = A.map((row) => row.slice());

  for (let col = 0; col < n; col++) {
    let pivot = col;
    let best = Math.abs(M[col][col]);
    for (let r = col + 1; r < n; r++) {
      if (Math.abs(M[r][col]) > best) {
        best = Math.abs(M[r][col]);
        pivot = r;
      }
    }
    if (pivot !== col) {
      [M[col], M[pivot]] = [M[pivot], M[col]];
      [B[col], B[pivot]] = [B[pivot], B[col]];
    }
    const p = M[col][col];
    for (let j = col; j < n; j++) M[col][j] /= p;
    for (let j = 0; j < k; j++) B[col][j] /= p;

    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const factor = M[r][col];
      if (factor === 0) continue;
      for (let j = col; j < n; j++) M[r][j] -= factor * M[col][j];
      for (let j = 0; j < k; j++) B[r][j] -= factor * B[col][j];
    }
  }
  return bIsVector ? B.map((row) => row[0]) : B;
}

function inv(A) {
  return solve(A, eye(A.length));
}

function trace(A) {
  let s = 0;
  for (let i = 0; i < A.length; i++) s += A[i][i];
  return s;
}

function diag(A) {
  return A.map((row, i) => row[i]);
}

// Lower-triangular Cholesky factor: L L^T = A. Adds a tiny diagonal jitter
// for numerical robustness (kernel Gram matrices can be near-singular),
// matching what numpy.random.multivariate_normal effectively tolerates.
function cholesky(A, jitter = 1e-9) {
  const n = A.length;
  const L = zeros(n, n);
  const M = addDiag(A, jitter);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = M[i][j];
      for (let p = 0; p < j; p++) sum -= L[i][p] * L[j][p];
      if (i === j) {
        L[i][j] = Math.sqrt(Math.max(sum, 0));
      } else {
        L[i][j] = L[j][j] !== 0 ? sum / L[j][j] : 0;
      }
    }
  }
  return L;
}

function add(A, B) {
  return A.map((row, i) => row.map((v, j) => v + B[i][j]));
}

function sub(A, B) {
  return A.map((row, i) => row.map((v, j) => v - B[i][j]));
}

function scale(A, s) {
  return A.map((row) => row.map((v) => v * s));
}

function kron(A, B) {
  const [ra, ca] = [A.length, A[0].length];
  const [rb, cb] = [B.length, B[0].length];
  const K = zeros(ra * rb, ca * cb);
  for (let i = 0; i < ra; i++)
    for (let j = 0; j < ca; j++)
      for (let p = 0; p < rb; p++)
        for (let q = 0; q < cb; q++) K[i * rb + p][j * cb + q] = A[i][j] * B[p][q];
  return K;
}

// Matrix exponential of (A*t) via scaling-and-squaring + truncated Taylor
// series. Good enough for the small (1-3 state), well-scaled systems these
// figures use -- not a general-purpose expm.
function expm(A, t = 1) {
  const n = A.length;
  const M0 = scale(A, t);
  let normEst = 0;
  for (const row of M0) for (const v of row) normEst = Math.max(normEst, Math.abs(v));
  const s = Math.max(0, Math.ceil(Math.log2((normEst || 1) * 4)));
  const M = scale(M0, 1 / Math.pow(2, s));

  let term = eye(n);
  let result = eye(n);
  for (let k = 1; k <= 25; k++) {
    term = scale(matMul(term, M), 1 / k);
    result = add(result, term);
  }
  for (let i = 0; i < s; i++) result = matMul(result, result);
  return result;
}

// Zero-order-hold discretization of continuous (A,B) at sample time Ts,
// via the standard augmented-matrix-exponential trick: builds
// M = [[A, B], [0, 0]], and Ad,Bd are the top blocks of expm(M*Ts).
// (C, D pass through unchanged for ZOH -- not returned here.)
function c2dZOH(A, B, Ts) {
  const n = A.length;
  const m = B[0].length;
  const M = zeros(n + m, n + m);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) M[i][j] = A[i][j];
    for (let j = 0; j < m; j++) M[i][n + j] = B[i][j];
  }
  const E = expm(M, Ts);
  const Ad = E.slice(0, n).map((row) => row.slice(0, n));
  const Bd = E.slice(0, n).map((row) => row.slice(n, n + m));
  return { Ad, Bd };
}

// Bilinear (Tustin) discretization of continuous (A,B,C,D) at sample time
// Ts, matching scipy.signal.cont2discrete(method="bilinear", alpha=0.5) /
// python-control's c2d(..., method="tustin").
function c2dTustin(A, B, C, D, Ts) {
  const n = A.length;
  const I = eye(n);
  const ima = sub(I, scale(A, 0.5 * Ts)); // I - 0.5*A*Ts
  const imaInv = inv(ima);
  const Ad = matMul(imaInv, add(I, scale(A, 0.5 * Ts)));
  const Bd = matMul(imaInv, scale(B, Ts));
  const Cd = matMul(C, imaInv);
  const Dd = add(D, scale(matMul(C, Bd), 0.5));
  return { Ad, Bd, Cd, Dd };
}

// Solve the discrete Lyapunov (Stein) equation P = A P A^T + Q for P,
// via vectorization: vec(P) = (I - kron(A,A))^{-1} vec(Q).
function discreteLyapunov(A, Q) {
  const n = A.length;
  const vecQ = [];
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) vecQ.push(Q[i][j]); // column-major, matches kron convention
  const K = kron(A, A);
  const lhs = sub(eye(n * n), K);
  const vecP = solve(lhs, vecQ);
  const P = zeros(n, n);
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) P[i][j] = vecP[j * n + i];
  return P;
}

function vecNorm(v) {
  return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
}

// Estimates the spectral radius (largest |eigenvalue|) of a real square
// matrix via power iteration on the growth rate of ||M^k v||. This still
// converges to the true spectral radius even when the dominant eigenvalue
// is a complex-conjugate pair (the vector direction keeps rotating inside
// that invariant subspace, but the norm's growth factor still tends to
// rho(M)) -- exactly what's needed here to test "is this matrix stable"
// without implementing a general eigenvalue solver.
function spectralRadius(M, iters = 200) {
  const n = M.length;
  let v = new Array(n).fill(0).map(() => Math.random() - 0.5);
  v = v.map((x) => x / vecNorm(v));
  let growth = 1;
  for (let i = 0; i < iters; i++) {
    const w = matVec(M, v);
    const norm = vecNorm(w);
    if (norm < 1e-300) return 0;
    if (i >= iters - 20) growth *= norm; // average the last few growth factors
    v = w.map((x) => x / norm);
  }
  return Math.pow(growth, 1 / 20);
}

window.linalg = {
  zeros,
  eye,
  matMul,
  matVec,
  transpose,
  addDiag,
  solve,
  inv,
  trace,
  diag,
  cholesky,
  add,
  sub,
  scale,
  kron,
  expm,
  c2dZOH,
  c2dTustin,
  discreteLyapunov,
  spectralRadius,
};
