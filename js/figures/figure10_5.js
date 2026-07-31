"use strict";

// Port of code_distribute/python/Figure10_5.py
window.figureLib = window.figureLib || {};
window.figureLib.figure10_5 = function (outputGrid, params) {
  const idN = Math.max(1, Math.round(params.id_n));
  const sensorRange = Math.max(1, Math.round(params.sensor_range));
  const tMax = Math.max(1, Math.round(params.t_max));

  // Fixed constants (not parametrized), matching the Python source exactly.
  const gridSize = 100;
  const beta = 1;
  const iniState = Array.from({ length: idN }, () => [10, 10]);

  const directions = [
    [1, 0],
    [-1, 0],
    [0, 1],
    [0, -1],
  ];

  // value[I][J], I,J in 0..100 (101 x 101), matching numpy's value[I, J].
  function initValue() {
    const g = gridSize;
    const value = Array.from({ length: g + 1 }, () => new Array(g + 1).fill(0));
    for (let I = 0; I <= g; I++) {
      for (let J = 0; J <= g; J++) {
        const v =
          (1.5 * Math.exp(-((I - 20) ** 2) / 1000 - ((J - 20) ** 2) / 1000) +
            Math.exp(-((I - 20) ** 2) / 1000 - ((J - 90) ** 2) / 500) +
            Math.exp(-((I - 40) ** 2) / 1000 - ((J - 50) ** 2) / 500) +
            Math.exp(-((I - 70) ** 2) / 300 - ((J - 70) ** 2) / 500) +
            1.5 * Math.exp(-((I - 80) ** 2) / 300 - ((J - 40) ** 2) / 500) +
            Math.exp(-((I - 50) ** 2) / 800 - ((J - 50) ** 2) / 800) +
            1.2 * Math.exp(-((I - 80) ** 2) / 200 - ((J - 20) ** 2) / 200) +
            Math.exp(-((I - 90) ** 2) / 200 - ((J - 10) ** 2) / 200)) /
          10.0;
        value[I][J] = v;
      }
    }
    return value;
  }

  // Sum of value[i][j] over the sub-grid around (X, Y), restricted to cells
  // that have at least one agent within sensor range (near_mask semantics:
  // a cell counts if the MINIMUM squared distance across all agents is < r2).
  function sumValueNearAgents(state, sRange, X, Y, value) {
    const r2 = sRange * sRange;
    const g = gridSize;
    const evaluationRange = sRange + 3;
    const Xi = Math.trunc(X);
    const Yi = Math.trunc(Y);
    const i0 = Math.max(0, Xi - evaluationRange);
    const i1 = Math.min(g, Xi + evaluationRange);
    const j0 = Math.max(0, Yi - evaluationRange);
    const j1 = Math.min(g, Yi + evaluationRange);

    let sum = 0;
    const nAgents = state.length;
    for (let i = i0; i <= i1; i++) {
      for (let j = j0; j <= j1; j++) {
        let minD2 = Infinity;
        for (let a = 0; a < nAgents; a++) {
          const dx = state[a][0] - i;
          const dy = state[a][1] - j;
          const d2 = dx * dx + dy * dy;
          if (d2 < minD2) minD2 = d2;
        }
        if (minD2 < r2) sum += value[i][j];
      }
    }
    return sum;
  }

  function allAgentsLocallyOptimal(state, value, sRange) {
    const idn = state.length;
    for (let ID = 0; ID < idn; ID++) {
      const x0 = state[ID][0];
      const y0 = state[ID][1];
      const baseValue = sumValueNearAgents(state, sRange, x0, y0, value);

      for (let d = 0; d < directions.length; d++) {
        const dx = directions[d][0];
        const dy = directions[d][1];
        const proposed = state.map((row) => row.slice());
        proposed[ID][0] = Math.min(gridSize, Math.max(0, x0 + dx));
        proposed[ID][1] = Math.min(gridSize, Math.max(0, y0 + dy));

        // NOTE: evaluation center is the ORIGINAL (x0, y0), not the moved
        // position, even though `proposed` (the moved state) is what's summed.
        const newValue = sumValueNearAgents(proposed, sRange, x0, y0, value);
        if (newValue > baseValue) return false;
      }
    }
    return true;
  }

  function simulation(value, iniStateArr, Tmax, deterministic) {
    const idn = iniStateArr.length;
    let state = iniStateArr.map((row) => row.slice());
    const xList = Array.from({ length: idn }, () => new Array(Tmax).fill(0));
    const yList = Array.from({ length: idn }, () => new Array(Tmax).fill(0));
    let lastK = Tmax - 1;

    for (let k = 0; k < Tmax; k++) {
      for (let a = 0; a < idn; a++) {
        xList[a][k] = state[a][0];
        yList[a][k] = state[a][1];
      }

      const agentId = Math.floor(Math.random() * idn);
      const direction = directions[Math.floor(Math.random() * 4)];
      const candidateState = state.map((row) => row.slice());
      candidateState[agentId][0] = Math.min(gridSize, Math.max(0, candidateState[agentId][0] + direction[0]));
      candidateState[agentId][1] = Math.min(gridSize, Math.max(0, candidateState[agentId][1] + direction[1]));

      // NOTE: x, y are the CURRENT (pre-move) position of the chosen agent,
      // used as the evaluation center for both current and candidate sums.
      const x = state[agentId][0];
      const y = state[agentId][1];
      const currentValue = sumValueNearAgents(state, sensorRange, x, y, value);
      const candidateValue = sumValueNearAgents(candidateState, sensorRange, x, y, value);

      lastK = k;

      if (deterministic) {
        if (candidateValue > currentValue) {
          state = candidateState;
          if (allAgentsLocallyOptimal(state, value, sensorRange)) {
            break;
          }
        }
      } else {
        const probMove = 1 / (1 + Math.exp(beta * (currentValue - candidateValue)));
        if (Math.random() < probMove) {
          state = candidateState;
        }
      }
    }

    return {
      state,
      xList: xList.map((row) => row.slice(0, lastK + 1)),
      yList: yList.map((row) => row.slice(0, lastK + 1)),
    };
  }

  function plotSimResult(body, value, Tmax, deterministic) {
    const { state, xList, yList } = simulation(value, iniState, Tmax, deterministic);
    const chart = window.plotlib.createChart(body, {
      xlim: [0, gridSize],
      ylim: [0, gridSize],
      xlabel: "$x$",
      ylabel: "$y$",
    });

    const iniX = iniState.map((p) => p[0]);
    const iniY = iniState.map((p) => p[1]);
    chart.scatter(iniX, iniY, { marker: "star", size: 8, color: "black", label: "初期位置" });

    const finalX = state.map((p) => p[0]);
    const finalY = state.map((p) => p[1]);
    chart.scatter(finalX, finalY, { size: 8, color: "red", label: "最終位置" });

    const cycle = ["blue", "orange", "green", "purple", "gray"];
    const nT = 100;
    const tArr = Array.from({ length: nT }, (_, i) => (2 * Math.PI * i) / (nT - 1));
    for (let i = 0; i < idN; i++) {
      const color = cycle[i % cycle.length];
      chart.line(xList[i], yList[i], { color, lineWidth: 1 });
      const cx = tArr.map((t) => 10 * Math.sin(t) + state[i][0]);
      const cy = tArr.map((t) => 10 * Math.cos(t) + state[i][1]);
      chart.line(cx, cy, { color, lineWidth: 1 });
    }
    chart.finish();
  }

  const value = initValue();

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.5(a) — 価値関数のヒートマップ");
    window.plotlib.createHeatmap(body, value, { width: 420, height: 380 });
  }

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.5(b) — 確率的方策によるシミュレーション結果");
    plotSimResult(body, value, tMax, false);
  }

  {
    const body = window.plotlib.makeCard(outputGrid, "Figure 10.5(c) — 決定論的方策によるシミュレーション結果");
    plotSimResult(body, value, tMax, true);
  }
};
