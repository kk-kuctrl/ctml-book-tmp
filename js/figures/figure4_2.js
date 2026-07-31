"use strict";

// Port of code_distribute/python/Figure4_2.py
// control.c2d(csys, T_c) with no method specified defaults to zero-order
// hold, which lib/linalg.js's c2dZOH implements via the augmented matrix
// exponential trick. C and D pass through unchanged under ZOH.
window.figureLib = window.figureLib || {};
window.figureLib.figure4_2 = function (outputGrid, params) {
  const Tc = params.t_c;
  // Steps are derived from the requested end TIME (seconds) rather than
  // exposed directly, since with Tc itself now adjustable, a fixed step
  // count would mean very different simulated durations depending on Tc.
  // No upper clamp: silently capping this would make the actual simulated
  // duration (kBar*Tc) quietly diverge from the requested end time for a
  // tiny Tc + large t_end combination, which is worse than the slider's
  // own min/max already bounding how large kBar can realistically get.
  const kBar = Math.max(10, Math.round(params.t_end / Tc));
  const d = params.d;
  const uMax = 3.2;

  const Ac = [
    [0, 4],
    [-3, 2],
  ];
  const Bc = [[0], [1]];
  const C = [[-0.3, -4]]; // unchanged by ZOH discretization
  const x0 = [[1], [0.5]];

  const { Ad, Bd } = window.linalg.c2dZOH(Ac, Bc, Tc);

  function simulate(uRule) {
    const y = new Array(kBar).fill(0);
    const u = new Array(kBar).fill(0);
    let xk = x0.map((row) => row.slice());
    for (let k = 0; k < kBar; k++) {
      y[k] = window.linalg.matMul(C, xk)[0][0];
      u[k] = uRule(y[k]);
      xk = window.linalg.add(window.linalg.matMul(Ad, xk), window.linalg.scale(Bd, u[k]));
    }
    return { y, u };
  }

  function stairsXY(values) {
    const xs = [0];
    const ys = [values[0]];
    for (let i = 0; i < values.length; i++) {
      xs.push(i + 1);
      ys.push(values[i]);
      if (i + 1 < values.length) {
        xs.push(i + 1);
        ys.push(values[i + 1]);
      }
    }
    return { xs, ys };
  }

  function drawPanel(title, y, u) {
    const body = window.plotlib.makeCard(outputGrid, title);
    const chart = window.plotlib.createChart(body, {
      xlim: [0, kBar],
      ylim: [-uMax, uMax],
      xlabel: "$k$",
      legendLoc: "upper right",
    });
    if (u) {
      const us = stairsXY(u);
      chart.line(us.xs, us.ys, { color: [0.7, 0.7, 0.7], lineWidth: 1, label: "u_k" });
    }
    const ys = stairsXY(y);
    chart.line(ys.xs, ys.ys, { color: "black", lineWidth: 1, label: "y_k" });
    chart.finish();
  }

  // (1) u_k = y_k
  {
    const { y } = simulate((yk) => yk);
    drawPanel("Figure 4.2(1) — u_k = y_k（そのままフィードバック）", y, null);
  }

  // (2) u_k = Q(y_k)  (nearest multiple of d)
  {
    const { y, u } = simulate((yk) => Math.floor((yk + d / 2) / d) * d);
    drawPanel("Figure 4.2(2) — u_k = Q(y_k)（量子化）", y, u);
  }

  // (3) u_k = Q(y_k + z_k), z_k ~ Uniform(-0.5, 0.5)  (dithered quantization)
  {
    const { y, u } = simulate((yk) => {
      const zk = Math.random() - 0.5;
      return Math.floor((yk + zk + d / 2) / d) * d;
    });
    drawPanel("Figure 4.2(3) — u_k = Q(y_k + z_k)（ディザ量子化）", y, u);
  }
};
