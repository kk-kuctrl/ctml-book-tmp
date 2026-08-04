"use strict";

// Port of code_distribute/python/Figure11_4.py
window.figureLib = window.figureLib || {};
window.figureLib.figure11_4 = function (outputGrid, params) {
  const nK = Math.max(1, Math.round(params.n_k));

  const C = [params.c1, params.c2, params.c3, params.c4];
  const alpha = [params.alpha1, params.alpha2, params.alpha3, params.alpha4];
  const palette = ["blue", "orange", "green", "red"];
  const nSetting = C.length;
  const pIni = 1;
  const mean = params.noise_mean;
  const variance = params.noise_var;

  // z_k is the "observation" whose mean/quantile p_k is estimating.
  function sampleNoise() {
    if (params.noise_dist === "laplace") return mean + window.rnd.laplaceSample(Math.sqrt(variance / 2));
    return mean + Math.sqrt(variance) * window.rnd.randn();
  }

  // Mean-tracking (y = p-z) targets E[z] = mean directly. Exercise 11.5
  // (裾分布の推定) instead tracks the c-QUANTILE p* of z (P(z<=p*) = c; the
  // median is just c=0.5), via y_k = 1{z_k<=p_k} - c -- this only ever needs
  // an expectation of an INDICATOR, not of z itself, so it stays
  // well-defined even for distributions with no finite mean (not offered
  // here, but the same recursion works for e.g. Cauchy noise too).
  function stepError(pk, z) {
    if (params.est_target === "quantile") return (z <= pk ? 1 : 0) - params.quantile_c;
    return pk - z;
  }

  // Peter Acklam's rational approximation of the standard normal inverse
  // CDF (relative error ~1e-9), needed to get the TRUE quantile value for
  // the reference line below -- JS has no built-in erfinv/probit.
  function normInvCDF(p) {
    const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.383577518672690e2, -3.066479806614716e1, 2.506628277459239e0];
    const b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
    const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838e0, -2.549732539343734e0, 4.374664141464968e0, 2.938163982698783e0];
    const d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996e0, 3.754408661907416e0];
    const plow = 0.02425;
    if (p < plow) {
      const q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    if (p > 1 - plow) {
      const q = Math.sqrt(-2 * Math.log(1 - p));
      return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    const q = p - 0.5;
    const r = q * q;
    return ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }

  // The true value being estimated -- the plain mean, or (for "quantile")
  // the analytic c-quantile of the chosen (Gaussian/Laplace) distribution,
  // so it can be drawn as a reference line to check convergence against.
  function trueTargetValue() {
    if (params.est_target !== "quantile") return mean;
    const c = params.quantile_c;
    if (params.noise_dist === "laplace") {
      const b = Math.sqrt(variance / 2);
      return c <= 0.5 ? mean + b * Math.log(2 * c) : mean - b * Math.log(2 * (1 - c));
    }
    return mean + Math.sqrt(variance) * normInvCDF(c);
  }
  const trueVal = trueTargetValue();

  // Fixed [-2,2] no longer fits once mean/variance/quantile can push the
  // true target (and hence the converged curves) well outside that -- pad
  // around whichever of the initial value / true target is further out.
  const yPad = 1.5;
  const yMin = Math.min(-2, trueVal - yPad, pIni - yPad);
  const yMax = Math.max(2, trueVal + yPad, pIni + yPad);

  const body = window.plotlib.makeCard(outputGrid, "Figure 11.4 — 確率的勾配降下法による平均推定");
  const chart = window.plotlib.createChart(body, {
    xlim: [-5, nK],
    ylim: [yMin, yMax],
    xlabel: "$k$",
    ylabel: "$p_k$",
    legendLoc: "upper right",
  });

  const ks = Array.from({ length: nK + 1 }, (_, k) => k);

  for (let setting = 0; setting < nSetting; setting++) {
    const p = new Array(nK + 1).fill(0);
    p[0] = pIni;
    for (let k = 0; k < nK; k++) {
      const pk = p[k];
      const z = sampleNoise();
      const y = stepError(pk, z);
      p[k + 1] = pk - (C[setting] / Math.pow(k + 1, alpha[setting])) * y;
    }
    chart.line(ks, p, { color: palette[setting % palette.length], lineWidth: 1, label: `$C=${C[setting]}, \\alpha=${alpha[setting]}$` });
  }

  chart.line([0, nK], [trueVal, trueVal], { color: "black", lineWidth: 1, dash: [6, 4], label: "True Value" });
  chart.scatter([0], [pIni], { color: "black", marker: "circle", size: 4, filled: false, label: "Initial Value" });

  chart.finish();
};
