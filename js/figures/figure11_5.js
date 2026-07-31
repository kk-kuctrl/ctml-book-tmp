"use strict";

// Port of code_distribute/python/Figure11_5.py
//
// The Python differentiates L1/L2 symbolically with sympy, but both are
// simple closed-form expressions -- differentiating by hand avoids needing
// a symbolic-math engine entirely:
//   L1(p) = p^2/2 + p/2                          -> L1'(p) = p + 1/2
//   L2(p) = -p*cos(10p)/20 + sin(10p)/200 - p/2   -> L2'(p) = p*sin(10p)/2 - 1/2
//     (the -cos(10p)/20 terms from differentiating each piece cancel out)
//   L = (L1+L2)/2                                 -> L'  = (L1'+L2')/2
window.figureLib = window.figureLib || {};
window.figureLib.figure11_5 = function (outputGrid, params) {
  const nK = Math.max(10, Math.round(params.n_k));

  function L1(p) {
    return (p * p) / 2 + p / 2;
  }
  function L2(p) {
    return (-p * Math.cos(10 * p)) / 20 + Math.sin(10 * p) / 200 - p / 2;
  }
  function L(p) {
    return (L1(p) + L2(p)) / 2;
  }
  function gradL1(p) {
    return p + 0.5;
  }
  function gradL2(p) {
    return (p * Math.sin(10 * p)) / 2 - 0.5;
  }
  function gradL(p) {
    return (gradL1(p) + gradL2(p)) / 2;
  }

  // ---------------------------------------------------------------
  // Figure 11.5(a): L, grad L, grad L1, grad L2 over p in [-5, 5].
  // ---------------------------------------------------------------
  {
    const pVals = [];
    for (let p = -5; p <= 5; p += 0.01) pVals.push(p);

    const body = window.plotlib.makeCard(outputGrid, "Figure 11.5(a) — L(p) とその勾配");
    const chart = window.plotlib.createChart(body, { xlim: [-5, 5], ylim: [-15, 15], xlabel: "{\\rm p}" });
    chart.line(pVals, pVals.map(L), { color: "blue", lineWidth: 1.5, label: "L({\\rm p})" });
    chart.line(pVals, pVals.map(gradL), { color: "orange", lineWidth: 1.5, label: "\\nabla L({\\rm p})" });
    chart.line(pVals, pVals.map(gradL1), { color: "green", lineWidth: 1.5, label: "\\nabla L_1({\\rm p})" });
    chart.line(pVals, pVals.map(gradL2), { color: "red", lineWidth: 1.5, label: "\\nabla L_2({\\rm p})" });
    chart.finish();
  }

  // ---------------------------------------------------------------
  // Figure 11.5(b): SGD on L using a randomly-chosen sub-gradient
  // (grad_L1 or grad_L2, 50/50) at each step, for three step-size settings.
  // ---------------------------------------------------------------
  {
    const C = [1.0, 1.0, 1.0];
    const alpha = [0.4, 0.8, 1.2];
    const palette = ["blue", "orange", "green"];
    const nSetting = C.length;
    const pIni = 1;

    const body = window.plotlib.makeCard(outputGrid, "Figure 11.5(b) — SGDによる零点探索");
    const chart = window.plotlib.createChart(body, { xlim: [0, nK], ylim: [-2, 2], xlabel: "$k$", ylabel: "$p_k$" });
    const ks = Array.from({ length: nK + 1 }, (_, k) => k);

    for (let setting = 0; setting < nSetting; setting++) {
      const p = new Array(nK + 1).fill(0);
      p[0] = pIni;
      for (let k = 0; k < nK; k++) {
        const pk = p[k];
        const y = Math.random() < 0.5 ? gradL1(pk) : gradL2(pk);
        p[k + 1] = pk - C[setting] / Math.pow(k + 1, alpha[setting]) * y;
      }
      chart.line(ks, p, { color: palette[setting % palette.length], lineWidth: 2, label: `\\alpha=${alpha[setting]}` });
    }
    chart.scatter([0], [pIni], { color: "black", marker: "circle", size: 8, label: "\\text{Initial Value}" });
    chart.finish();
  }
};
