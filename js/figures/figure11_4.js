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

  const body = window.plotlib.makeCard(outputGrid, "Figure 11.4 — 確率的勾配降下法による平均推定");
  const chart = window.plotlib.createChart(body, {
    xlim: [-5, nK],
    ylim: [-2, 2],
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
      const z = window.rnd.randn();
      const y = pk - z;
      p[k + 1] = pk - (C[setting] / Math.pow(k + 1, alpha[setting])) * y;
    }
    chart.line(ks, p, { color: palette[setting % palette.length], lineWidth: 1, label: `$C=${C[setting]}, \\alpha=${alpha[setting]}$` });
  }

  chart.scatter([0], [pIni], { color: "black", marker: "circle", size: 4, filled: false, label: "Initial Value" });

  chart.finish();
};
