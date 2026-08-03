"use strict";

// Port of code_distribute/python/Figure3_2.py
window.figureLib = window.figureLib || {};
window.figureLib.figure3_2 = function (outputGrid, params) {
  const kBar = Math.max(3, Math.round(params.k_bar));
  const nSample = Math.max(1, Math.round(params.n_sample));
  const ks = Array.from({ length: kBar }, (_, k) => k);

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.2(a) — 決定論的な遷移");
    const chart = window.plotlib.createChart(body, {
      xlim: [1, kBar - 1],
      ylim: [-1, 1],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });
    for (let i = 0; i < nSample; i++) {
      const x = new Array(kBar).fill(0);
      x[0] = -0.8 + 0.2 * i;
      for (let k = 0; k < kBar - 1; k++) x[k + 1] = x[k] + 0.1 * (x[k] - x[k] ** 3);
      chart.line(ks, x, { lineWidth: 1 });
    }
    chart.finish();
  }

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.2(b) — 一様雑音を伴う遷移");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar - 1],
      ylim: [-1, 1],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });
    for (let s = 0; s < nSample; s++) {
      const x = new Array(kBar).fill(0);
      x[0] = -0.5;
      for (let k = 0; k < kBar - 1; k++) {
        const rand = Math.random();
        x[k + 1] = x[k] + 0.1 * (x[k] - x[k] ** 3) + (rand - 0.5) * (1 - Math.abs(x[k]));
      }
      chart.line(ks, x, { lineWidth: 1 });
    }
    chart.finish();
  }
};
