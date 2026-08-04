"use strict";

// Port of code_distribute/python/Figure11_1.py -- a covariance-steering SDP
// (minimize sum of control-input variances subject to steering the state
// covariance from Sigma_0 to a target Sigma_10, via a PSD/LMI relaxation).
//
// This is the one figure in the book that's a genuine semidefinite program,
// not just "convex enough for gradient descent": cvxpy hands it to SCS.
// Naive gradient descent on a reduced (non-lifted) parametrization was tried
// first and reliably lands in the wrong local optimum (verified against the
// cvxpy/SCS reference: true optimum ~2.825, gradient descent from many
// random restarts converges to 3.0-20+). A hand-rolled ADMM on the lifted
// (Sigma,P,M) formulation was tried next and is *correctly formulated*
// (verified step-by-step against independent solves) but numerically
// unstable for this problem's fully-degenerate optimal face (every one of
// the 10 LMI blocks is rank-deficient at the optimum).
//
// What actually works, and is what's implemented here: a primal-dual
// (infeasible-start) Newton method with a log-det barrier for the PSD/LMI
// constraints -- the standard interior-point approach real SDP solvers use
// internally. Verified to match the cvxpy/SCS reference to 5+ significant
// figures in both the unconstrained (a) and mid-horizon-constrained (b)
// cases.
window.figureLib = window.figureLib || {};
window.figureLib.figure11_1 = function (outputGrid, params) {
  const L = window.linalg;
  const nSample = Math.max(1, Math.round(params.n_sample));

  const A = [
    [1, 0.1],
    [-0.3, 1],
  ];
  const Bvec = [0.7, 0.4]; // B is 2x1; kept as a plain vector
  const Bmat = [[0.7], [0.4]];
  const N = [
    [0.1, 0],
    [0, 0.1],
  ];
  // Sigma0/Sigma10/mu0/mu10 are free-text (like Figure 6.1/8.8's A) instead
  // of individual sliders -- a full 2x2 covariance (not just an isotropic
  // sigma0*I or a diagonal Sigma10) and an arbitrary 2-vector mean can be
  // typed in directly. Malformed input falls back to the textbook's
  // original defaults.
  function parseMatrix2x2(text) {
    const rows = String(text)
      .trim()
      .split(/[\n;]+/)
      .map((r) => r.trim())
      .filter((r) => r.length > 0);
    if (rows.length !== 2) return null;
    const M = rows.map((r) =>
      r
        .split(/[,\s]+/)
        .filter((s) => s.length > 0)
        .map(Number)
    );
    if (M.some((row) => row.length !== 2 || row.some((v) => !Number.isFinite(v)))) return null;
    return M;
  }
  function parseVec2(text) {
    const v = String(text)
      .trim()
      .split(/[,\s]+/)
      .filter((s) => s.length > 0)
      .map(Number);
    if (v.length !== 2 || v.some((x) => !Number.isFinite(x))) return null;
    return v;
  }
  let Sigma0 = parseMatrix2x2(params.sigma0_text);
  if (!Sigma0) Sigma0 = parseMatrix2x2("3,0\n0,3");
  let Sigma10 = parseMatrix2x2(params.sigma10_text);
  if (!Sigma10) Sigma10 = parseMatrix2x2("2,0\n0,0.5");
  let mu0 = parseVec2(params.mu0_text);
  if (!mu0) mu0 = [0, 0];
  let mu10 = parseVec2(params.mu10_text);
  if (!mu10) mu10 = [0, 0];
  const kBar = 10;
  const kMid = Math.min(kBar - 1, Math.max(1, Math.round(params.k_mid)));
  const midCap = params.mid_cap;

  // ---- Variable layout (57 scalars): -----------------------------------
  // Sigma_k for k=1..9 (Sigma_0, Sigma_10 are fixed constants), 3 entries
  // each (11,12,22): 27 scalars. P_k for k=0..9, 2 entries each: 20
  // scalars. M_k for k=0..9, 1 entry each: 10 scalars.
  const NX = 9 * 3 + 10 * 2 + 10;
  const idxSigma = (k) => {
    const base = (k - 1) * 3;
    return [base, base + 1, base + 2];
  };
  const idxP = (k) => {
    const base = 27 + k * 2;
    return [base, base + 1];
  };
  const idxM = (k) => 27 + 20 + k;

  function sigmaOf(x, k) {
    if (k === 0) return Sigma0;
    if (k === kBar) return Sigma10;
    const [i11, i12, i22] = idxSigma(k);
    return [
      [x[i11], x[i12]],
      [x[i12], x[i22]],
    ];
  }
  function Pof(x, k) {
    const [i1, i2] = idxP(k);
    return [x[i1], x[i2]];
  }
  function Mof(x, k) {
    return x[idxM(k)];
  }

  // Residual of Sigma_{k+1} == A Sigma_k A' + A P_k B' + B P_k' A' + B M_k B' + N,
  // for k=0..9, 3 independent (symmetric) entries each -> 30-dim residual.
  // This function is affine in x, so A_eq/b_eq below are extracted exactly
  // (no truncation error) via basis-vector evaluations rather than a
  // hand-derived Jacobian.
  function equalityResidual(x) {
    const res = new Array(30).fill(0);
    for (let k = 0; k < kBar; k++) {
      const Sk = sigmaOf(x, k);
      const Sk1 = sigmaOf(x, k + 1);
      const Pk = Pof(x, k);
      const Mk = Mof(x, k);
      // pred = A Sk A' + A Pk B' + B Pk' A' + Mk * (B B') + N
      const ASkAt = L.matMul(L.matMul(A, Sk), L.transpose(A));
      // A Pk B' + (A Pk B')': A*Pk is a 2-vector, outer with B gives the
      // cross term, and its transpose gives the symmetric partner.
      const APk = [A[0][0] * Pk[0] + A[0][1] * Pk[1], A[1][0] * Pk[0] + A[1][1] * Pk[1]];
      const cross = [
        [APk[0] * Bvec[0], APk[0] * Bvec[1]],
        [APk[1] * Bvec[0], APk[1] * Bvec[1]],
      ];
      const crossT = L.transpose(cross);
      const pred = L.zeros(2, 2);
      for (let i = 0; i < 2; i++)
        for (let j = 0; j < 2; j++) pred[i][j] = ASkAt[i][j] + cross[i][j] + crossT[i][j] + Mk * Bvec[i] * Bvec[j] + N[i][j];
      const diff = L.sub(Sk1, pred);
      res[3 * k + 0] = diff[0][0];
      res[3 * k + 1] = diff[0][1];
      res[3 * k + 2] = diff[1][1];
    }
    return res;
  }

  const r0 = equalityResidual(new Array(NX).fill(0));
  const bEq = r0.map((v) => -v);
  const AEq = L.zeros(30, NX);
  for (let i = 0; i < NX; i++) {
    const e = new Array(NX).fill(0);
    e[i] = 1;
    const ri = equalityResidual(e);
    for (let r = 0; r < 30; r++) AEq[r][i] = ri[r] - r0[r];
  }
  const AEqT = L.transpose(AEq);

  // Objective: minimize sum_k M_k.
  const c = new Array(NX).fill(0);
  for (let k = 0; k < kBar; k++) c[idxM(k)] = 1;

  // Each LMI_k = [[Sigma_k, P_k],[P_k', M_k]] (3x3 symmetric), 6 independent
  // entries. blockVarMap says which x-index (or fixed constant) backs each
  // entry -- only k=0's Sigma part is constant (Sigma_0), everything else
  // is a free variable.
  const ENTRIES = ["11", "12", "22", "13", "23", "33"];
  const ENTRY_POS = { "11": [0, 0], "12": [0, 1], "22": [1, 1], "13": [0, 2], "23": [1, 2], "33": [2, 2] };
  function EMatrix(entry) {
    const [a, b] = ENTRY_POS[entry];
    const M = L.zeros(3, 3);
    if (a === b) M[a][b] = 1;
    else {
      M[a][b] = 1;
      M[b][a] = 1;
    }
    return M;
  }
  const EMATS = {};
  for (const e of ENTRIES) EMATS[e] = EMatrix(e);

  function blockVarMap(k) {
    let s11, s12, s22;
    if (k === 0) {
      s11 = s12 = s22 = ["const"];
    } else {
      const [i11, i12, i22] = idxSigma(k);
      s11 = ["var", i11];
      s12 = ["var", i12];
      s22 = ["var", i22];
    }
    const [p1i, p2i] = idxP(k);
    return { "11": s11, "12": s12, "22": s22, "13": ["var", p1i], "23": ["var", p2i], "33": ["var", idxM(k)] };
  }
  const blocks = [];
  for (let k = 0; k < kBar; k++) blocks.push(blockVarMap(k));

  function buildLMI(x, k) {
    const Sk = sigmaOf(x, k);
    const Pk = Pof(x, k);
    const Mk = Mof(x, k);
    return [
      [Sk[0][0], Sk[0][1], Pk[0]],
      [Sk[1][0], Sk[1][1], Pk[1]],
      [Pk[0], Pk[1], Mk],
    ];
  }

  // Strict-PD check via plain (no-jitter) Cholesky: fails as soon as a pivot
  // isn't safely positive, which is exactly the barrier's domain boundary.
  function isPD(M) {
    const n = M.length;
    const Lm = L.zeros(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = M[i][j];
        for (let p = 0; p < j; p++) sum -= Lm[i][p] * Lm[j][p];
        if (i === j) {
          if (sum <= 1e-12) return false;
          Lm[i][j] = Math.sqrt(sum);
        } else {
          Lm[i][j] = sum / Lm[j][j];
        }
      }
    }
    return true;
  }
  function inDomain(x, caseB, midIdx, cap) {
    for (let k = 0; k < kBar; k++) if (!isPD(buildLMI(x, k))) return false;
    if (caseB && cap - x[midIdx] <= 1e-12) return false;
    return true;
  }

  // Gradient/Hessian of phi(x) = c'x - mu * sum_k logdet(LMI_k(x)) [+ case b:
  // - mu*log(cap - Sigma_{kMid}(2,2))], restricted to the free (var-tagged)
  // coordinates. Uses the exact identities grad_e = -tr(Xinv @ E_e) and
  // Hess_{e,e'} = tr(Xinv E_e Xinv E_e') for a log-det barrier -- verified
  // against finite-difference gradients/Hessians while developing this.
  function barrierGradHess(x, mu, caseB, midIdx, cap) {
    const grad = new Array(NX).fill(0);
    const Hess = L.zeros(NX, NX);
    for (let k = 0; k < kBar; k++) {
      const M3 = buildLMI(x, k);
      const Xinv = L.inv(M3);
      const XinvE = {};
      for (const e of ENTRIES) XinvE[e] = L.matMul(Xinv, EMATS[e]);
      for (const e of ENTRIES) {
        const tag = blocks[k][e];
        if (tag[0] !== "var") continue;
        grad[tag[1]] += mu * -L.trace(XinvE[e]);
      }
      for (const e1 of ENTRIES) {
        const tag1 = blocks[k][e1];
        if (tag1[0] !== "var") continue;
        const i = tag1[1];
        for (const e2 of ENTRIES) {
          const tag2 = blocks[k][e2];
          if (tag2[0] !== "var") continue;
          const j = tag2[1];
          Hess[i][j] += mu * L.trace(L.matMul(XinvE[e1], XinvE[e2]));
        }
      }
    }
    if (caseB) {
      const s = cap - x[midIdx];
      grad[midIdx] += mu * (1 / s);
      Hess[midIdx][midIdx] += mu * (1 / (s * s));
    }
    for (let i = 0; i < NX; i++) grad[i] += c[i];
    return { grad, Hess };
  }

  function vecNorm(v) {
    return Math.sqrt(v.reduce((s, x) => s + x * x, 0));
  }

  // Infeasible-start Newton method (Boyd & Vandenberghe sec. 10.3.1): jointly
  // drives the KKT residual (dual stationarity + primal feasibility) to
  // zero, without needing an exactly-feasible starting point.
  function infeasibleNewton(x0, nu0, mu, caseB, midIdx, cap, maxIter) {
    let x = x0.slice();
    let nu = nu0.slice();
    for (let it = 0; it < maxIter; it++) {
      const { grad, Hess } = barrierGradHess(x, mu, caseB, midIdx, cap);
      const AeqTnu = L.matVec(AEqT, nu);
      const rDual = grad.map((g, i) => g + AeqTnu[i]);
      const rPrimal = L.matVec(AEq, x).map((v, i) => v - bEq[i]);
      const rnorm = Math.sqrt(rDual.reduce((s, v) => s + v * v, 0) + rPrimal.reduce((s, v) => s + v * v, 0));
      if (rnorm < 1e-10) break;

      const total = NX + 30;
      const KKT = L.zeros(total, total);
      for (let i = 0; i < NX; i++) for (let j = 0; j < NX; j++) KKT[i][j] = Hess[i][j];
      for (let i = 0; i < NX; i++)
        for (let j = 0; j < 30; j++) {
          KKT[i][NX + j] = AEqT[i][j];
          KKT[NX + j][i] = AEq[j][i];
        }
      const rhs = new Array(total).fill(0);
      for (let i = 0; i < NX; i++) rhs[i] = -rDual[i];
      for (let j = 0; j < 30; j++) rhs[NX + j] = -rPrimal[j];
      const sol = L.solve(KKT, rhs);
      const dx = sol.slice(0, NX);
      const dnu = sol.slice(NX);

      let t = 1.0;
      while (!inDomain(
        x.map((v, i) => v + t * dx[i]),
        caseB,
        midIdx,
        cap
      )) {
        t *= 0.5;
        if (t < 1e-14) break;
      }
      for (let tries = 0; tries < 60; tries++) {
        const xn = x.map((v, i) => v + t * dx[i]);
        const nun = nu.map((v, i) => v + t * dnu[i]);
        if (!inDomain(xn, caseB, midIdx, cap)) {
          t *= 0.5;
          continue;
        }
        const gh = barrierGradHess(xn, mu, caseB, midIdx, cap);
        const AeqTnun = L.matVec(AEqT, nun);
        const rd = gh.grad.map((g, i) => g + AeqTnun[i]);
        const rp = L.matVec(AEq, xn).map((v, i) => v - bEq[i]);
        const rn = Math.sqrt(rd.reduce((s, v) => s + v * v, 0) + rp.reduce((s, v) => s + v * v, 0));
        if (rn <= (1 - 0.01 * t) * rnorm || t < 1e-14) {
          x = xn;
          nu = nun;
          break;
        }
        t *= 0.5;
      }
    }
    return { x, nu };
  }

  // Solve the SDP via a barrier method: start from a comfortably-interior
  // (equality-infeasible is fine) point and geometrically shrink mu, doing a
  // handful of Newton steps at each level. 25 levels x <=30 Newton steps
  // converges to 5+ significant figures against the cvxpy/SCS reference.
  function solveCovarianceSteering(caseB) {
    const midIdx = idxSigma(kMid)[2];
    // A generic interior starting point: isotropic covariances around the
    // average of Sigma_0/Sigma_10's scale, comfortably PD regardless of what
    // was typed in, plus a strictly-feasible constant M_k.
    const initDiag = Math.max(0.5, (Sigma0[0][0] + Sigma0[1][1] + Sigma10[0][0] + Sigma10[1][1]) / 4);
    let x = new Array(NX).fill(0);
    for (let k = 1; k < kBar; k++) {
      const [i11, i12, i22] = idxSigma(k);
      x[i11] = initDiag;
      x[i12] = 0;
      x[i22] = initDiag;
    }
    for (let k = 0; k < kBar; k++) {
      x[idxM(k)] = Math.max(5, initDiag * 2);
    }
    if (caseB) x[midIdx] = Math.min(initDiag, midCap * 0.6); // strictly below the cap from the start
    let nu = new Array(30).fill(0);
    let mu = 5.0;
    for (let outer = 0; outer < 25; outer++) {
      const res = infeasibleNewton(x, nu, mu, caseB, midIdx, midCap, 30);
      x = res.x;
      nu = res.nu;
      mu *= 0.5;
    }
    const SigmaOpt = [];
    for (let k = 0; k <= kBar; k++) SigmaOpt.push(sigmaOf(x, k));
    const Popt = [];
    const Mopt = [];
    for (let k = 0; k < kBar; k++) {
      Popt.push(Pof(x, k));
      Mopt.push(Mof(x, k));
    }
    return { SigmaOpt, Popt, Mopt };
  }

  function eig2x2Sym(M) {
    const a = M[0][0],
      b = M[0][1],
      d = M[1][1];
    const mid = (a + d) / 2;
    const rad = Math.sqrt(((a - d) / 2) ** 2 + b * b);
    const eig1 = mid + rad,
      eig2 = mid - rad;
    let v1;
    if (Math.abs(b) > 1e-12) v1 = [b, eig1 - a];
    else v1 = a >= d ? [1, 0] : [0, 1];
    const norm = Math.hypot(v1[0], v1[1]) || 1;
    v1 = [v1[0] / norm, v1[1] / norm];
    return { eig1, eig2, v1 };
  }
  function ellipseXY(Sigma) {
    const { eig1, eig2, v1 } = eig2x2Sym(Sigma);
    const rot = Math.atan2(v1[1], v1[0]);
    const ax0 = 2 * Math.sqrt(Math.max(eig1, 0));
    const ax1 = 2 * Math.sqrt(Math.max(eig2, 0));
    const n = 60;
    const xs = new Array(n),
      ys = new Array(n);
    for (let i = 0; i < n; i++) {
      const th = (2 * Math.PI * i) / (n - 1);
      const ct = Math.cos(th),
        st = Math.sin(th);
      xs[i] = ax0 * ct * Math.cos(rot) - ax1 * st * Math.sin(rot);
      ys[i] = ax0 * ct * Math.sin(rot) + ax1 * st * Math.cos(rot);
    }
    return { xs, ys };
  }

  // matplotlib's default "tab10" color cycle -- matches what the original
  // Python gets for free (each ax.plot() call in the ellipse loop advances
  // to the next cycle color automatically); wraps if kBar+1 > 10.
  const TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];

  // Per Computer Exercise 11.1 (substituting Upsilon_k=diag(O,I) into eq.
  // 11.7): E[||u_k||^2] = Trace(M_k) + ||ubar_k||^2. The covariance-steering
  // SDP above already handles the Trace(M_k) part (assuming mean 0
  // throughout); ubar_k (the mean of u_k) is a COMPLETELY SEPARATE free
  // decision variable with no shared constraint against Sigma_k/M_k/K_k --
  // so xbar_k's steering from mu0 to mu10 is just a plain deterministic
  // minimum-energy control problem on xbar_{k+1}=A*xbar_k+Bvec*ubar_k, solved
  // via the same minimum-norm (pseudo-inverse) construction as the
  // covariance side, but for reachability rather than a barrier method.
  function computeMeanTrajectory(mu0In, mu10In) {
    const powA = [L.eye(2)];
    for (let k = 1; k <= kBar; k++) powA.push(L.matMul(powA[k - 1], A));
    // cols[k] = A^(kBar-1-k) @ Bvec: ubar_k's contribution to xbar_kBar.
    const cols = Array.from({ length: kBar }, (_, k) => L.matVec(powA[kBar - 1 - k], Bvec));
    const freeResp = L.matVec(powA[kBar], mu0In);
    const target = [mu10In[0] - freeResp[0], mu10In[1] - freeResp[1]];

    const MMt = L.zeros(2, 2);
    for (let r = 0; r < 2; r++)
      for (let c = 0; c < 2; c++) {
        let s = 0;
        for (let k = 0; k < kBar; k++) s += cols[k][r] * cols[k][c];
        MMt[r][c] = s;
      }
    const beta = L.matVec(L.inv(MMt), target);
    const ubar = cols.map((col) => col[0] * beta[0] + col[1] * beta[1]);

    let xbar = mu0In.slice();
    const meanTraj = [xbar.slice()];
    for (let k = 0; k < kBar; k++) {
      xbar = [A[0][0] * xbar[0] + A[0][1] * xbar[1] + Bvec[0] * ubar[k], A[1][0] * xbar[0] + A[1][1] * xbar[1] + Bvec[1] * ubar[k]];
      meanTraj.push(xbar.slice());
    }
    const meanCost = ubar.reduce((s, u) => s + u * u, 0); // sum ||ubar_k||^2 == sum gamma_k at optimum
    return { meanTraj, ubar, meanCost };
  }
  // Independent of caseB (the mid-horizon constraint only touches
  // Sigma_k), so compute it once and reuse for both panels.
  const { meanTraj, meanCost } = computeMeanTrajectory(mu0, mu10);

  function renderCase(title, caseB) {
    const { SigmaOpt, Popt, Mopt } = solveCovarianceSteering(caseB);
    // Total E[sum ||u_k||^2] = sum Trace(M_k) + sum ||ubar_k||^2 (Computer
    // Exercise 11.1's decomposition) -- the covariance-steering part plus
    // the mean-steering part, computed independently above.
    const optimalCost = Mopt.reduce((s, v) => s + v, 0) + meanCost;
    const Kgain = [];
    for (let k = 0; k < kBar; k++) {
      const SigInv = L.inv(SigmaOpt[k]);
      // K_k = (Sigma_k^{-1} P_k)' = P_k' Sigma_k^{-1} (P_k is a column, so
      // this is a 1x2 row vector).
      const t = L.matVec(SigInv, Popt[k]);
      Kgain.push(t);
    }

    const body = window.plotlib.makeCard(outputGrid, title);
    const kHalf = kBar; // matches python's xlim [10,0] essentially spanning 0..10
    const chart3d = window.plotlib.createChart3D(body, {
      width: 576,
      height: 350,
      margin: 28,
      xStretch: 1.2,
      elev: 15,
      azim: 35,
      xlabel: "$k$",
      ylabel: "$({\\rm x})_1$",
      zlabel: "$({\\rm x})_2$",
      xlim: [0, kHalf],
      ylim: [-10, 10],
      zlim: [-10, 10],
    });
    const kVals = Array.from({ length: kBar + 1 }, (_, i) => i);

    for (let s = 0; s < nSample; s++) {
      // x is the zero-mean FLUCTUATION around meanTraj -- covariance
      // steering (Kgain) only ever acts on this part; the mean itself is
      // added in separately (meanTraj already bakes in its own feedforward).
      let x = window.rnd.mvnSample([0, 0], Sigma0);
      const y1 = [meanTraj[0][0] + x[0]];
      const y2 = [meanTraj[0][1] + x[1]];
      for (let k = 0; k < kBar; k++) {
        const Kk = Kgain[k];
        const uk = Kk[0] * x[0] + Kk[1] * x[1];
        const v = window.rnd.mvnSample([0, 0], N);
        x = [A[0][0] * x[0] + A[0][1] * x[1] + Bvec[0] * uk + v[0], A[1][0] * x[0] + A[1][1] * x[1] + Bvec[1] * uk + v[1]];
        y1.push(meanTraj[k + 1][0] + x[0]);
        y2.push(meanTraj[k + 1][1] + x[1]);
      }
      chart3d.addLine(kVals, y1, y2, [0.6, 0.6, 0.6]);
    }
    for (let k = 0; k <= kBar; k++) {
      const { xs, ys } = ellipseXY(SigmaOpt[k]);
      const kConst = new Array(xs.length).fill(k);
      // Only the initial (k=0) and final (k=kBar) distributions are the ones
      // being steered between -- bold those two so they stand out from the
      // intermediate-time ellipses.
      const isEndpoint = k === 0 || k === kBar;
      const xsShift = xs.map((v) => v + meanTraj[k][0]);
      const ysShift = ys.map((v) => v + meanTraj[k][1]);
      chart3d.addLine(kConst, xsShift, ysShift, TAB10[k % TAB10.length], { lineWidth: isEndpoint ? 3 : 1.2 });
    }
    chart3d.finish();

    const costLine = document.createElement("p");
    costLine.className = "hint";
    costLine.style.margin = "8px 0 0";
    costLine.textContent = `最適コスト（Σ Trace(M_k) + Σ‖ū_k‖²、制御入力の平均2乗の総和）: ${optimalCost.toFixed(4)}`;
    body.appendChild(costLine);
  }

  renderCase("Figure 11.1(a) — 分散推移の最適化（中間時刻の制約なし）", false);
  renderCase(`Figure 11.1(b) — 分散推移の最適化（k=${kMid}でΣ(2,2)≤${midCap}の制約あり）`, true);
};
