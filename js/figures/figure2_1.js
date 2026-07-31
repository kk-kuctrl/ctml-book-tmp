"use strict";

// Port of code_distribute/python/Figure2_1.py
function normalPdf(x, mean, std) {
  const z = (x - mean) / std;
  return Math.exp(-0.5 * z * z) / (std * Math.sqrt(2 * Math.PI));
}

function laplacePdf(x, loc, scale) {
  return Math.exp(-Math.abs(x - loc) / scale) / (2 * scale);
}

function linspace(a, b, n) {
  return Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
}

window.figureLib = window.figureLib || {};
window.figureLib.figure2_1 = function (outputGrid, params) {
  const N = 600;

  // Figure 2.1(a) -- folded onto |x|, so the Normal's mean stays fixed at 0.
  {
    const aVar = params.a_var;
    const aLapScale = params.a_lap_scale;
    const x = linspace(0, 10, N);
    const pdfNormal = x.map((xi) => normalPdf(xi, 0, Math.sqrt(aVar)));
    const pdfLap = x.map((xi) => laplacePdf(xi, 0, aLapScale));

    const body = window.plotlib.makeCard(outputGrid, "Figure 2.1(a)");
    const chart = window.plotlib.createChart(body, {
      xlim: [0, 10],
      ylim: [1e-12, 1],
      yscale: "log",
      xlabel: "|{\\rm x}|",
      legendLoc: "upper right",
    });
    chart.line(x, pdfNormal, { color: "orange", lineWidth: 2, label: `\\mathcal{N}(0, ${aVar})` });
    chart.line(x, pdfLap, { color: "purple", lineWidth: 2, label: `{\\rm Lap}(0, ${aLapScale})` });
    chart.finish();
  }

  // Figure 2.1(b)
  {
    const mean1 = params.b_mean1,
      var1 = params.b_var1;
    const mean2 = params.b_mean2,
      var2 = params.b_var2;
    const w = params.b_mixw;
    const lapLoc = params.b_lap_loc,
      lapScale = params.b_lap_scale;

    const x = linspace(-10, 10, N);
    const std1 = Math.sqrt(var1),
      std2 = Math.sqrt(var2);
    const sumMean = mean1 + mean2,
      sumVar = var1 + var2;

    const pdf1 = x.map((xi) => normalPdf(xi, mean1, std1));
    const pdf2 = x.map((xi) => normalPdf(xi, mean2, std2));
    const pdfSum = x.map((xi) => normalPdf(xi, sumMean, Math.sqrt(sumVar)));
    const pdfMix = x.map((xi) => w * normalPdf(xi, mean1, std1) + (1 - w) * normalPdf(xi, mean2, std2));
    const pdfLap = x.map((xi) => laplacePdf(xi, lapLoc, lapScale));

    const body = window.plotlib.makeCard(outputGrid, "Figure 2.1(b)");
    const chart = window.plotlib.createChart(body, {
      xlim: [-5, 5],
      ylim: [0, 1],
      xlabel: "{\\rm x}",
      legendLoc: "upper right",
    });
    chart.line(x, pdf1, { color: "blue", lineWidth: 2, label: "x_1 \\sim \\mathcal{N}(\\mu_1, \\sigma_1^2)" });
    chart.line(x, pdf2, { color: "orange", lineWidth: 2, label: "x_2 \\sim \\mathcal{N}(\\mu_2, \\sigma_2^2)" });
    chart.line(x, pdfSum, { color: "green", lineWidth: 2, label: "(x_1+x_2) \\sim \\mathcal{N}(\\mu_1+\\mu_2, \\sigma_1^2+\\sigma_2^2)" });
    chart.line(x, pdfMix, { color: "red", lineWidth: 2, label: "(\\varphi_{x_1}+\\varphi_{x_2})/2" });
    chart.line(x, pdfLap, { color: "purple", lineWidth: 2, label: "x_3 \\sim {\\rm Lap}(\\mu_3, b_3)" });
    chart.finish();
  }
};
