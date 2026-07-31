"use strict";

// Port of code_distribute/python/Figure4_1.py
window.figureLib = window.figureLib || {};
window.figureLib.figure4_1 = function (outputGrid, params) {
  const sysA = 1.1; // system gain (Python module-level `a`)
  const alpha = -0.2; // control law parameter
  const xMax = 5; // view range

  const kBar = Math.max(1, Math.round(params.k_bar));
  const nSample = Math.max(1, Math.round(params.n_sample));
  const ks = Array.from({ length: kBar + 1 }, (_, k) => k);

  {
    const x = new Array(kBar + 1).fill(0);
    x[0] = 1;
    for (let k = 0; k < kBar; k++) {
      const u = alpha * Math.pow(sysA + alpha, k);
      x[k + 1] = sysA * x[k] + u;
    }

    const body = window.plotlib.makeCard(outputGrid, "Figure 4.1(a) — Deterministic reference trajectory");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar],
      ylim: [-xMax, xMax],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });
    chart.line(ks, x, { color: "black", lineWidth: 1 });
    chart.finish();
  }

  {
    const x = new Array(kBar + 1).fill(0);
    x[0] = 1;
    const xV = Array.from({ length: nSample }, () => {
      const row = new Array(kBar + 1).fill(0);
      row[0] = 1;
      return row;
    });

    for (let k = 0; k < kBar; k++) {
      const v = window.rnd.randnVec(nSample);
      const u = alpha * Math.pow(sysA + alpha, k);
      x[k + 1] = sysA * x[k] + u;
      for (let s = 0; s < nSample; s++) {
        xV[s][k + 1] = sysA * xV[s][k] + u + 0.1 * v[s];
      }
    }

    const body = window.plotlib.makeCard(outputGrid, "Figure 4.1(b) — Trajectories under noise");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar],
      ylim: [-xMax, xMax],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });
    for (let s = 0; s < nSample; s++) {
      chart.line(ks, xV[s], { color: [0.7, 0.7, 0.7], lineWidth: 0.5 });
    }
    chart.line(ks, x, { color: "black", lineWidth: 1 });
    chart.finish();
  }
};
