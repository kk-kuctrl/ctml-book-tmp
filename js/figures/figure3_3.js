"use strict";

// Port of code_distribute/python/Figure3_3.py
window.figureLib = window.figureLib || {};
window.figureLib.figure3_3 = function (outputGrid, params) {
  const kBar = Math.max(1, Math.round(params.k_bar));
  const nSample = Math.max(1, Math.round(params.n_sample));
  const a = params.a;
  const xMax = 10;
  const ks = Array.from({ length: kBar + 1 }, (_, k) => k);

  // Figure 3.3(a): x[k+1] = a*x[k] + step, where `step` is drawn ONCE per
  // sample before the loop and reused (unchanged) at every time step --
  // this is the subtle detail that distinguishes (a) from (b) below.
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.3(a) — 固定ノイズを伴う遷移");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar],
      ylim: [-xMax, xMax],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });

    const x0 = window.rnd.randnVec(nSample);
    const step = window.rnd.randnVec(nSample);
    const x = Array.from({ length: nSample }, () => new Array(kBar + 1).fill(0));
    for (let s = 0; s < nSample; s++) x[s][0] = x0[s];
    for (let k = 0; k < kBar; k++) {
      for (let s = 0; s < nSample; s++) {
        x[s][k + 1] = a * x[s][k] + step[s];
      }
    }
    for (let s = 0; s < nSample; s++) {
      chart.line(ks, x[s], { lineWidth: 2, alpha: 0.8 });
    }
    chart.finish();
  }

  // Figure 3.3(b): x[k+1] = a*x[k] + noise*sqrt(3), where `noise` is
  // resampled fresh at EVERY step (not held fixed like (a) above).
  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 3.3(b) — 独立なノイズを伴う遷移");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar],
      ylim: [-xMax, xMax],
      xlabel: "$k$",
      ylabel: "$x_k$",
    });

    const x0 = window.rnd.randnVec(nSample);
    const x = Array.from({ length: nSample }, () => new Array(kBar + 1).fill(0));
    for (let s = 0; s < nSample; s++) x[s][0] = x0[s];
    for (let k = 0; k < kBar; k++) {
      const noise = window.rnd.randnVec(nSample);
      for (let s = 0; s < nSample; s++) {
        x[s][k + 1] = a * x[s][k] + noise[s] * Math.sqrt(3);
      }
    }
    for (let s = 0; s < nSample; s++) {
      chart.line(ks, x[s], { lineWidth: 2, alpha: 0.8 });
    }
    chart.finish();
  }
};
