"use strict";

// Port of code_distribute/python/Figure12.py
// Note: the Python reseeds with `np.random.seed(13)` inside figure12_3and4
// so the data points are reproducible across the three calls. This port
// deliberately does not replicate that seeding -- every Run uses fresh
// randomness, consistent with the rest of this site (see lib/random.js).
window.figureLib = window.figureLib || {};
window.figureLib.figure12 = function (outputGrid, params) {
  const c = params.c;
  const sBarSmall = Math.max(3, Math.round(params.s_bar));
  const sBarLarge = Math.max(3, Math.round(params.s_bar) * 5);

  function linspace(a, b, n) {
    if (n === 1) return [a];
    const out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = a + ((b - a) * i) / (n - 1);
    return out;
  }

  // K[i][j] = exp( -(x_i - y_j)^2 / c^2 )
  function gaussianKernel(x, y, cc) {
    const K = new Array(x.length);
    const c2 = cc * cc;
    for (let i = 0; i < x.length; i++) {
      const row = new Array(y.length);
      for (let j = 0; j < y.length; j++) {
        const d = x[i] - y[j];
        row[j] = Math.exp(-(d * d) / c2);
      }
      K[i] = row;
    }
    return K;
  }

  // K[i][j] = min(x_i, y_j)
  function minKernel(x, y) {
    const K = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const row = new Array(y.length);
      for (let j = 0; j < y.length; j++) row[j] = Math.min(x[i], y[j]);
      K[i] = row;
    }
    return K;
  }

  function matSub(A, B) {
    return A.map((row, i) => row.map((v, j) => v - B[i][j]));
  }

  function figure12_2(Nx, cWidth, cardTitle) {
    const x = linspace(0, 1, Nx);
    const K = gaussianKernel(x, x, cWidth);
    const priorMean = x.slice();

    const body = window.plotlib.makeCard(outputGrid, cardTitle);
    const chart = window.plotlib.createChart(body, {
      xlim: [0, 1],
      ylim: [-3, 3],
      xlabel: "{\\rm x}",
      legendLoc: "upper left",
    });

    for (let i = 0; i < 5; i++) {
      const sample = window.rnd.mvnSample(priorMean, K);
      chart.line(x, sample, { color: "gray", alpha: 0.8 });
    }
    chart.line(x, priorMean, { label: "\\mu({\\rm x})" });

    const xSD = window.linalg.diag(K).map((v) => Math.sqrt(Math.max(0, v)));
    const lo = priorMean.map((m, i) => m - xSD[i]);
    const hi = priorMean.map((m, i) => m + xSD[i]);
    chart.fillBetween(x, lo, hi, { color: "blue", alpha: 0.25 });

    chart.finish();
  }

  function figure12_3and4(Nx, sBar, kernel, cardTitle) {
    const xData = new Array(sBar);
    for (let i = 0; i < sBar; i++) xData[i] = Math.random() * Math.random();
    const x = linspace(0, 1, Nx);

    let K, Ktmp, kx;
    if (kernel === "Gaussian") {
      K = gaussianKernel(xData, xData, 0.1);
      Ktmp = gaussianKernel(x, x, 0.1);
      kx = gaussianKernel(x, xData, 0.1); // Nx x sBar
    } else {
      K = minKernel(xData, xData);
      Ktmp = minKernel(x, x);
      kx = minKernel(x, xData); // Nx x sBar
    }

    const muX = (v) => v; // mu(x) = x
    const fnX = (v) => Math.sin(4 * Math.PI * v);

    const fnList = x.map(fnX);

    const sigma = 0.1;
    const yData = xData.map((xi) => fnX(xi) + sigma * window.rnd.randn());
    const priorMean = xData.map(muX);

    // solved1 = solve(K + sigma^2 I, y_data - prior_mean)  -- length sBar
    const rhsVec = yData.map((v, i) => v - priorMean[i]);
    const Kreg = window.linalg.addDiag(K, sigma * sigma);
    const solved1 = window.linalg.solve(Kreg, rhsVec);

    // mean(x_i) = mu(x_i) + kx[i] . solved1
    const meanVec = window.linalg.matVec(kx, solved1);
    const mean = meanVec.map((v, i) => v + x[i]);

    // Kpost = Ktmp - kx @ solve(K + sigma^2 I, kx^T)
    const kxT = window.linalg.transpose(kx); // sBar x Nx
    const solvedMat = window.linalg.solve(Kreg, kxT); // sBar x Nx
    const Kpost = matSub(Ktmp, window.linalg.matMul(kx, solvedMat)); // Nx x Nx
    const xSD = window.linalg.diag(Kpost).map((v) => Math.sqrt(Math.max(0, v)));

    const body = window.plotlib.makeCard(outputGrid, cardTitle);
    const chart = window.plotlib.createChart(body, {
      xlim: [0, 1],
      ylim: [-1.5, 3.0],
      xlabel: "{\\rm x}",
      legendLoc: "upper left",
    });

    chart.line(x, x.map(muX), { color: "red", label: "\\mu({\\rm x})" });
    chart.line(x, fnList, { color: "black", dash: [6, 4], label: "\\sin(4\\pi {\\rm x})" });
    chart.scatter(xData, yData, { color: "blue", filled: false, label: "{\\rm y}_s" });

    const lo = mean.map((m, i) => m - xSD[i]);
    const hi = mean.map((m, i) => m + xSD[i]);
    chart.fillBetween(x, lo, hi, { color: "blue", alpha: 0.25 });
    chart.line(x, mean, { color: "blue", dash: [8, 3, 2, 3], label: "\\mu ({\\rm x}|\\mathcal{D})" });

    chart.finish();
  }

  figure12_2(100, c, "Figure 12.2(a)");
  figure12_2(100, 1.0, "Figure 12.2(b)");

  figure12_3and4(100, sBarSmall, "Gaussian", "Figure 12.3(a)");
  figure12_3and4(100, sBarLarge, "Gaussian", "Figure 12.3(b)");

  figure12_3and4(100, sBarSmall, "min", "Figure 12.4(a)");
  figure12_3and4(100, sBarLarge, "min", "Figure 12.4(b)");
};
