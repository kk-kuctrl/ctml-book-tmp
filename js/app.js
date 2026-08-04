"use strict";

// Figures ported to plain JS from code_distribute/python/. Each entry's
// `render(container, params)` (see figures/*.js) draws directly and
// synchronously -- no loading step, unlike the Pyodide site.
const FIGURES = [
  {
    key: "figure2_1",
    chapter: 2,
    title: "Figure 2.1",
    note: "正規分布とラプラス分布の比較",
    params: [
      { key: "a_var", label: "2.1(a) N(0, σ²) の分散 σ²", default: 2, min: 0.1, max: 10, step: 0.1 },
      { key: "a_lap_scale", label: "2.1(a) Laplace のスケール b", default: 1, min: 0.1, max: 5, step: 0.1 },
      { key: "b_mean1", label: "2.1(b) x₁ の平均 μ₁", default: -2, min: -5, max: 5, step: 0.1 },
      { key: "b_var1", label: "2.1(b) x₁ の分散 σ₁²", default: 1, min: 0.1, max: 5, step: 0.1 },
      { key: "b_mean2", label: "2.1(b) x₂ の平均 μ₂", default: 2, min: -5, max: 5, step: 0.1 },
      { key: "b_var2", label: "2.1(b) x₂ の分散 σ₂²", default: 1, min: 0.1, max: 5, step: 0.1 },
      { key: "b_mixw", label: "2.1(b) 混合比 w（x₁側の重み）", default: 0.5, min: 0, max: 1, step: 0.05 },
      { key: "b_lap_loc", label: "2.1(b) Laplace の位置 μ", default: 0, min: -5, max: 5, step: 0.1 },
      { key: "b_lap_scale", label: "2.1(b) Laplace のスケール b", default: 1, min: 0.1, max: 5, step: 0.1 },
    ],
  },
  {
    key: "figure3_1",
    chapter: 3,
    title: "Figure 3.1",
    note: "分布の時間発展（3D）とサンプル軌道",
    params: [
      { key: "k_bar_a", label: "3.1(a) ステップ数", default: 9, min: 1, max: 30, step: 1 },
      { key: "k_bar_b", label: "3.1(b) ステップ数", default: 30, min: 5, max: 100, step: 1 },
      { key: "n_sample_b", label: "3.1(b) サンプル数", default: 10, min: 1, max: 30, step: 1 },
      {
        type: "range2",
        keyLow: "x0_min",
        keyHigh: "x0_max",
        label: "初期値の範囲（一様分布、(a)(b)共通）",
        default: [-0.5, 0.5],
        min: -1,
        max: 1,
        step: 0.05,
        newRow: true,
      },
    ],
  },
  {
    key: "figure3_2",
    chapter: 3,
    title: "Figure 3.2",
    note: "決定論的・確率的な状態遷移",
    params: [
      { key: "k_bar", label: "ステップ数", default: 50, min: 5, max: 200, step: 1 },
      { key: "n_sample", label: "サンプル数", default: 10, min: 1, max: 30, step: 1 },
      {
        type: "range2",
        keyLow: "x0_min_b",
        keyHigh: "x0_max_b",
        label: "(b) 初期値の範囲（一様分布）",
        default: [-0.5, 0.5],
        min: -2,
        max: 2,
        step: 0.05,
        newRow: true,
      },
    ],
  },
  {
    key: "figure3_3",
    chapter: 3,
    title: "Figure 3.3",
    note: "線形フィードバックのサンプル軌道と収束性",
    method:
      "ゲイン a は線形遷移 x_{k+1}=a*x_k+雑音 の係数（フィードバック強度）です。|a|<1 なら軌道は原点付近に収束・分布が定常化し、|a|≥1 なら発散します。",
    params: [
      { key: "k_bar", label: "ステップ数", default: 10, min: 5, max: 200, step: 1 },
      { key: "n_sample", label: "サンプル数", default: 20, min: 1, max: 50, step: 1 },
      { key: "a", label: "ゲイン a", default: 0.5, min: -1.5, max: 1.5, step: 0.05 },
    ],
  },
  {
    key: "figure3_4",
    chapter: 3,
    title: "Figure 3.4",
    note: "周波数重みフィルタと有色雑音生成（1次系のTustin離散化＋離散リアプノフ方程式）",
    method:
      "離散リアプノフ方程式 P = APA' + Q は vec(P) = (I - A⊗A)^{-1} vec(Q) とベクトル化し、線形連立方程式として解いています（scipyのsolve_discrete_lyapunovの代わり）。Tustin（双一次変換）離散化も行列演算で直接計算しています。",
    params: [
      { key: "a", label: "雑音モデル F(s)=1/(s+a) の a", default: 0.3, min: 0.05, max: 2, step: 0.05 },
      { key: "ts", label: "サンプル時間 Ts", default: 0.1, min: 0.01, max: 0.5, step: 0.01 },
      { key: "k_bar", label: "3.4(b) ステップ数", default: 200, min: 20, max: 500, step: 10 },
    ],
  },
  {
    key: "figure4_1",
    chapter: 4,
    title: "Figure 4.1",
    note: "フィードフォワード制御のもとでの不安定システムの状態推移",
    params: [
      { key: "k_bar", label: "ステップ数", default: 50, min: 5, max: 200, step: 1 },
      { key: "n_sample", label: "4.1(b) サンプル数", default: 20, min: 1, max: 50, step: 1 },
    ],
  },
  {
    key: "figure4_2",
    chapter: 4,
    title: "Figure 4.2",
    note: "量子化フィードバックと2状態系のZOH離散化（行列指数関数）",
    method:
      "行列指数関数 exp(At) はスケーリング＆スクエアリング法（行列を2の冪で縮小してから冪級数展開し、2乗を繰り返して元のスケールに戻す）で計算しています（scipyのexpmの代わり）。",
    params: [
      { key: "t_end", label: "終了時刻（秒）", default: 8, min: 0.1, max: 10, step: 0.1 },
      { key: "t_c", label: "サンプル時間 T_c", default: 0.01, min: 0.0001, max: 0.5, log: true },
      { key: "d", label: "量子化幅 d", default: 1.0, min: 0.2, max: 3, step: 0.1 },
    ],
  },
  {
    key: "figure5_3",
    chapter: 5,
    title: "Figure 5.3",
    note: "LQR/LQG制御と有色雑音下でのカルマンフィルタ推定",
    method: "Figure 3.4と同じ離散リアプノフ方程式（Kronecker積によるベクトル化）とTustin離散化を、カルマンフィルタの定常共分散の計算にも使っています。",
    params: [
      { key: "a", label: "雑音モデル F(s)=1/(s+a) の a", default: 0.3, min: 0.05, max: 2, step: 0.05 },
      { key: "k_bar", label: "ステップ数", default: 300, min: 50, max: 600, step: 10 },
    ],
  },
  {
    key: "figure5_4",
    chapter: 5,
    title: "Figure 5.4",
    note: "外れ値に頑健な状態推定",
    method:
      "L1正則化された最小二乗を、射影付きAdamとε-アニーリング（平滑化パラメータεを1e-1から1e-10まで徐々に小さくする継続法／graduated non-convexity）で解いています（cvxpyの代わり）。固定εのままだと真のL1解のようなスパースな残差構造が再現できませんでした。",
    params: [
      { key: "k_bar", label: "ステップ数", default: 60, min: 20, max: 150, step: 5 },
      {
        key: "init_dist",
        type: "select",
        label: "初期状態 x_0 の分布",
        default: "gaussian",
        newRow: true,
        options: [
          { value: "gaussian", label: "正規分布" },
          { value: "laplace", label: "ラプラス分布" },
          { value: "uniform", label: "一様分布" },
        ],
      },
      { key: "init_mean", label: "初期状態 x_0 の平均（各成分共通）", default: 1, min: -3, max: 3, step: 0.1 },
      { key: "init_scale", label: "初期状態 x_0 のばらつき（一様:台の半幅／それ以外:分散）", default: 1, min: 0.05, max: 5, step: 0.05 },
      {
        key: "proc_dist",
        type: "select",
        label: "外乱 v の分布",
        default: "gaussian",
        newRow: true,
        options: [
          { value: "gaussian", label: "正規分布" },
          { value: "laplace", label: "ラプラス分布" },
          { value: "uniform", label: "一様分布" },
        ],
      },
      { key: "proc_scale", label: "外乱 v の分散（一様のときは無視）", default: 1, min: 0.05, max: 5, step: 0.05 },
      { key: "proc_support", label: "外乱 v の台（±、分布によらず指定・一様のときはそのまま範囲）", default: 1, min: 0.1, max: 5, step: 0.1 },
      { key: "proc_unbounded", type: "checkbox", label: "外乱 v を非有界に（台を無視、一様のときは無効）", default: false },
      {
        key: "noise_dist",
        type: "select",
        label: "観測雑音 w の分布",
        default: "laplace",
        newRow: true,
        options: [
          { value: "gaussian", label: "正規分布" },
          { value: "laplace", label: "ラプラス分布" },
          { value: "uniform", label: "一様分布" },
        ],
      },
      { key: "noise_scale", label: "観測雑音 w のばらつき（一様:台の半幅／それ以外:分散）", default: 2, min: 0.05, max: 10, step: 0.05 },
    ],
  },
  {
    key: "figure6_1",
    chapter: 6,
    title: "Figure 6.1",
    note: "割引LQRの価値反復・方策反復による収束比較",
    method: "割引LQRの最適値 P_opt は、リカッチ再帰式を（十分な回数）収束するまで反復して求めています（python-controlのdlqrの代わり）。",
    params: [
      { key: "beta", label: "割引率 β", default: 0.95, min: 0.5, max: 0.999, step: 0.005 },
      { key: "n_iter", label: "反復回数", default: 11, min: 5, max: 30, step: 1 },
      {
        key: "A_text",
        type: "matrix",
        rows: 3,
        label: "A行列（NxN、次元は自由・行は改行区切り、値はカンマ区切り）",
        default: "0.8,0.9,0.86\n0.3,0.25,1\n0.1,0.55,0.5",
        newRow: true,
      },
      { key: "Q_text", type: "matrix", rows: 3, label: "Q行列（Aの次元に合わせる）", default: "1,0,0\n0,1,0\n0,0,1" },
      { key: "B_text", type: "matrix", rows: 3, label: "B ベクトル（縦、Aの次元に合わせる）", default: "1\n0\n0" },
    ],
  },
  {
    key: "figure7_1",
    chapter: 7,
    title: "Figure 7.1",
    note: "多項式基底によるベイズ回帰",
    params: [
      { key: "n_sample_a", label: "7.1(a) サンプル関数の数", default: 20, min: 1, max: 60, step: 1 },
      { key: "n_sample_b", label: "7.1(b) 事後サンプル数", default: 20, min: 1, max: 60, step: 1 },
      { key: "s_bar", label: "7.1(b) データ数", default: 20, min: 3, max: 30, step: 1 },
    ],
  },
  {
    key: "figure7_4",
    chapter: 7,
    title: "Figure 7.4",
    note: "最小二乗・Ridge・Lasso回帰の比較",
    method:
      "Lasso回帰は座標降下法（各係数を閉形式のsoft-threshold更新で1つずつ最適化）で解いています。（この多項式基底はヴァンデルモンド型で悪条件のため、近接勾配法（FISTA）は収束が遅く不安定でした。）",
    params: [
      { key: "sigma_sq", label: "正則化の重み σ²", default: 0.01, min: 0.001, max: 0.2, step: 0.001 },
      { key: "s_bar", label: "データ数", default: 30, min: 5, max: 60, step: 1 },
    ],
  },
  {
    key: "figure8_2",
    chapter: 8,
    title: "Figure 8.2",
    note: "逐次最小二乗法によるシステム同定",
    params: [
      { key: "k_bar_a", label: "8.2(a) ステップ数", default: 1000, min: 100, max: 5000, step: 100 },
      { key: "k_bar_b", label: "8.2(b) ステップ数", default: 5000, min: 500, max: 20000, step: 500 },
      { key: "k_bar_cdef", label: "8.2(c-f) ステップ数", default: 20000, min: 1000, max: 100000, step: 1000 },
    ],
  },
  {
    key: "figure8_8",
    chapter: 8,
    title: "Figure 8.8",
    note: "マルコフモデルのサンプル軌道（離散リアプノフ方程式で初期分布を設定）",
    method:
      "初期分布の共分散は離散リアプノフ方程式（Kronecker積によるベクトル化）で、系の安定化係数はべき乗法によるスペクトル半径推定で計算しています（scipy/numpy.linalg.eigの代わり）。",
    params: [
      { key: "k_bar", label: "ステップ数", default: 200, min: 20, max: 500, step: 10 },
      { key: "n_sample", label: "サンプル数", default: 5, min: 1, max: 10, step: 1 },
      {
        key: "A_text",
        type: "matrix",
        rows: 5,
        label: "A行列（NxN、次元は自由・不安定な場合のみスペクトル半径0.95に自動正規化）",
        default:
          "0.42,0.72,0.00,0.30,0.15\n0.09,0.19,0.35,0.40,0.54\n0.42,0.69,0.20,0.88,0.03\n0.67,0.42,0.56,0.14,0.20\n0.80,0.97,0.31,0.69,0.88",
        newRow: true,
      },
      {
        key: "C_text",
        type: "matrix",
        rows: 2,
        label: "C行列（出力 y=Cx、行数=2固定・列数はAの次元に合わせる）",
        default: "1,0,0,0,0\n0,1,0,0,0",
      },
    ],
  },
  {
    key: "figure9_3and4",
    chapter: 9,
    title: "Figure 9.3 & 9.4",
    note: "最小二乗TD学習によるLQRの強化学習（RLSはFigure8.2と同じ更新式）",
    params: [
      { key: "k_update", label: "ゲイン更新あたりのステップ数", default: 50, min: 10, max: 150, step: 5 },
      { key: "n_path", label: "サンプルパス数", default: 20, min: 5, max: 50, step: 1 },
      { key: "sigma_a", label: "9.3(a) 探索雑音レベル σ", default: 2, min: 0.1, max: 20, step: 0.1 },
      { key: "sigma_b", label: "9.3(b)/9.4 探索雑音レベル σ", default: 10, min: 0.1, max: 20, step: 0.1 },
      { key: "iter_gain", label: "9.4 ゲイン更新回数", default: 5, min: 1, max: 20, step: 1 },
    ],
  },
  {
    key: "figure10_2and4",
    chapter: 10,
    title: "Figure 10.2(b) & 10.4",
    note: "KL制御による最適方策とオンライン逆強化学習での状態コスト推定。P(i,j) は状態jから状態iへの遷移重みで、同じ列（同じ遷移元j）内で自動的に合計1へ正規化されます（0にすればその遷移を削除でき、トポロジーも自由に変更可）。状態数はPの行数（正方行列である必要あり）で自由に決まり、コストベクトルはその次元に合わせて下さい。",
    method:
      "L_KL・L_IRLはどちらも凸関数（アフィン項＋log-sum-exp項）なので、scipy.optimize.minimizeの代わりに解析的勾配を手計算してAdamで最小化しています。",
    params: [
      { key: "k_bar_2", label: "10.2(b) ステップ数", default: 50, min: 10, max: 300, step: 5 },
      { key: "k_bar_4", label: "10.4 ステップ数", default: 1000, min: 100, max: 3000, step: 100 },
      {
        key: "P_text",
        type: "matrix",
        rows: 4,
        label: "遷移重み行列 P（NxN、次元は自由・行=遷移先i／列=遷移元j、各列は自動正規化）",
        default: "0.3333,0.3333,0,0\n0,0.3333,0.3333,0\n0,0.3333,0.3333,0.3333\n0.6667,0,0.3333,0.6667",
        newRow: true,
      },
      {
        key: "cost_text",
        type: "text",
        label: "コストベクトル ℓ（Pの次元に合わせる、カンマ区切り）",
        default: "1,2,3,4",
      },
    ],
  },
  {
    key: "figure10_5",
    chapter: 10,
    title: "Figure 10.5",
    note: "マルチエージェントの分散最適化（重め: 数秒かかる場合あり）。センサー範囲は各エージェントがその半径内にいる（自分または他の）エージェントの位置を基準に価値場を感知・合算できる距離で、これが各エージェントの移動判断（現在位置と移動先で感知できる価値の合計を比較）に使われます。",
    params: [
      { key: "id_n", label: "エージェント数", default: 12, min: 2, max: 20, step: 1 },
      { key: "sensor_range", label: "センサー範囲", default: 8, min: 2, max: 15, step: 1 },
      { key: "t_max", label: "シミュレーション長", default: 30000, min: 500, max: 50000, step: 500 },
      { key: "init_x", label: "初期位置 x（全エージェント共通）", default: 10, min: 0, max: 100, step: 1, newRow: true },
      { key: "init_y", label: "初期位置 y（全エージェント共通）", default: 10, min: 0, max: 100, step: 1 },
      {
        key: "value_field",
        type: "select",
        label: "価値場の種類",
        default: "original",
        newRow: true,
        options: [
          { value: "original", label: "教科書の例" },
          { value: "two_peaks", label: "2峰（離れた2箇所）" },
          { value: "ripple", label: "同心円状のさざ波" },
        ],
      },
      { key: "show_sensor_circles", type: "checkbox", label: "アニメ中もセンサー範囲の丸を表示（終了後は常に表示）", default: true },
    ],
  },
  {
    key: "figure11_1",
    chapter: 11,
    title: "Figure 11.1",
    note: "状態分布を目標分布へ操舵する最適制御（最適輸送制御）",
    method:
      "半正定値計画問題（SDP）を内点法（infeasible-start Newton法＋対数バリア）で解いています。等式制約（共分散の遷移関係）はKKT連立方程式で厳密に扱い、半正定値制約は-log detバリアで滑らかな目的関数に変換し、バリア係数μを5から3e-7まで25段階で減衰させながらNewton法を反復します（cvxpy/SCSの代わり）。",
    params: [
      { key: "n_sample", label: "サンプル軌道数", default: 20, min: 1, max: 50, step: 1 },
      { key: "k_mid", label: "制約を課す時刻 k", default: 5, min: 1, max: 9, step: 1, newRow: true },
      { key: "mid_cap", label: "制約幅（Σ_k(2,2) の上限）", default: 0.5, min: 0.05, max: 5, step: 0.05 },
      {
        key: "sigma0_text",
        type: "matrix",
        rows: 2,
        label: "初期分散 Σ_0（2x2、対称正定値行列）",
        default: "3,0\n0,3",
        newRow: true,
      },
      {
        key: "sigma10_text",
        type: "matrix",
        rows: 2,
        label: "終端分散 Σ_10（2x2、対称正定値行列）",
        default: "2,0\n0,0.5",
      },
      { key: "mu0_text", type: "text", label: "初期分布の平均 mu_0（カンマ区切り）", default: "0,0" },
      { key: "mu10_text", type: "text", label: "終端分布の平均 mu_10（カンマ区切り）", default: "0,0" },
    ],
  },
  {
    key: "figure11_4",
    chapter: 11,
    title: "Figure 11.4",
    note: "確率的勾配降下法による平均推定",
    params: [
      { key: "n_k", label: "反復回数", default: 100, min: 10, max: 5000, step: 10 },
      { key: "c1", label: "設定1: C", default: 0.6, min: 0.1, max: 2, step: 0.05 },
      { key: "alpha1", label: "設定1: α", default: 1.0, min: 0.1, max: 2, step: 0.05 },
      { key: "c2", label: "設定2: C", default: 0.6, min: 0.1, max: 2, step: 0.05 },
      { key: "alpha2", label: "設定2: α", default: 0.3, min: 0.1, max: 2, step: 0.05 },
      { key: "c3", label: "設定3: C", default: 1.0, min: 0.1, max: 2, step: 0.05 },
      { key: "alpha3", label: "設定3: α", default: 1.0, min: 0.1, max: 2, step: 0.05 },
      { key: "c4", label: "設定4: C", default: 1.0, min: 0.1, max: 2, step: 0.05 },
      { key: "alpha4", label: "設定4: α", default: 1.5, min: 0.1, max: 2, step: 0.05 },
      {
        key: "noise_dist",
        type: "select",
        label: "確率変数 z の分布",
        default: "gaussian",
        newRow: true,
        options: [
          { value: "gaussian", label: "正規分布" },
          { value: "laplace", label: "ラプラス分布（裾がやや重い）" },
        ],
      },
      { key: "noise_mean", label: "z の平均", default: 0, min: -3, max: 3, step: 0.1 },
      { key: "noise_var", label: "z の分散", default: 1, min: 0.05, max: 5, step: 0.05 },
      {
        key: "est_target",
        type: "select",
        label: "推定対象（y_kの定義）",
        default: "mean",
        newRow: true,
        options: [
          { value: "mean", label: "平均（y=p-z）" },
          { value: "quantile", label: "分位点（y=𝟙[z≤p]-c）" },
        ],
      },
      { key: "quantile_c", label: "分位点 c（P(z≤p*)=c・裾の端まで指定可）", default: 0.5, min: 0.001, max: 0.999, step: 0.001 },
    ],
  },
  {
    key: "figure11_5",
    chapter: 11,
    title: "Figure 11.5",
    note: "ミニバッチSGDによる零点探索（手計算の微分でsympy不要）",
    params: [{ key: "n_k", label: "反復回数", default: 2000, min: 100, max: 10000, step: 100 }],
  },
  {
    key: "figure12",
    chapter: 12,
    title: "Figure 12",
    note: "ガウス過程回帰（事前分布・事後分布）",
    params: [
      { key: "c", label: "12.2 カーネル幅 c", default: 0.1, min: 0.01, max: 2, step: 0.01 },
      { key: "s_bar", label: "12.3/12.4 データ数", default: 10, min: 3, max: 60, step: 1 },
    ],
  },
];

const els = {
  figureList: document.getElementById("figure-list"),
  figTitle: document.getElementById("fig-title"),
  figNote: document.getElementById("fig-note"),
  status: document.getElementById("status"),
  runBtn: document.getElementById("run-btn"),
  paramsPanel: document.getElementById("params-panel"),
  outputGrid: document.getElementById("output-grid"),
  outputHint: document.getElementById("output-hint"),
  methodNote: document.getElementById("method-note"),
  methodNoteText: document.getElementById("method-note-text"),
};

let currentEntry = null;
let currentParamValues = {};

function setStatus(text, kind) {
  els.status.textContent = text;
  els.status.className = "status" + (kind ? " " + kind : "");
}

// Chapter titles from the textbook's table of contents (parenthetical
// qualifiers in the original titles, e.g. "(主に自律的な)確率システム", are
// dropped here per request -- just the core topic name).
const CHAPTER_NAMES = {
  1: "はじめに",
  2: "確率と統計",
  3: "確率システム",
  4: "最適制御",
  5: "状態推定",
  6: "漸近挙動",
  7: "機械学習",
  8: "システム同定",
  9: "強化学習",
  10: "マルコフ連鎖と定常状態",
  11: "数理最適化",
  12: "ガウス過程回帰",
};

function buildSidebar() {
  const groups = new Map();
  for (const entry of FIGURES) {
    if (!groups.has(entry.chapter)) groups.set(entry.chapter, []);
    groups.get(entry.chapter).push(entry);
  }
  els.figureList.innerHTML = "";
  for (const chapter of [...groups.keys()].sort((a, b) => a - b)) {
    const label = document.createElement("div");
    label.className = "chapter-label";
    label.textContent = `第${chapter}章 ${CHAPTER_NAMES[chapter] || ""}`;
    els.figureList.appendChild(label);
    for (const entry of groups.get(chapter)) {
      const btn = document.createElement("button");
      btn.className = "figure-item";
      btn.dataset.key = entry.key;
      btn.textContent = entry.title;
      btn.addEventListener("click", () => selectEntry(entry));
      els.figureList.appendChild(btn);
    }
  }
}

// Debounced so dragging a slider doesn't re-run the (synchronous) figure
// function on every "input" tick -- cheap figures wouldn't notice, but
// Figure 10.5 can take ~100-300ms and would visibly stutter otherwise.
let autoRunTimer = null;
function scheduleAutoRun() {
  clearTimeout(autoRunTimer);
  autoRunTimer = setTimeout(runCurrent, 80);
}

function renderParamsPanel(entry) {
  els.paramsPanel.style.display = "block";
  els.paramsPanel.innerHTML = "";
  const title = document.createElement("div");
  title.className = "params-title";
  title.textContent = "パラメータ（ドラッグすると自動で再実行されます）";
  els.paramsPanel.appendChild(title);

  // A param with `newRow: true` starts a fresh grid container instead of
  // joining the current one -- a real block-level break, so the group
  // always starts on its own line regardless of the auto-fill column count
  // (unlike e.g. `grid-column: 1` on an item, which only forces a new row
  // when the columns happen to line up that way).
  let grid = null;
  for (const p of entry.params) {
    if (!grid || p.newRow) {
      grid = document.createElement("div");
      grid.className = "params-grid";
      els.paramsPanel.appendChild(grid);
    }
    const field = document.createElement("div");
    field.className = "params-field";

    const header = document.createElement("div");
    header.className = "params-field-header";
    const labelSpan = document.createElement("span");
    labelSpan.className = "params-label";
    labelSpan.textContent = p.label;
    header.appendChild(labelSpan);

    let input;
    if (p.type === "matrix" || p.type === "text") {
      // Free-text entry for things a slider can't reasonably express (a
      // matrix, a vector) -- the figure itself parses the string, with a
      // fallback to its default on malformed input.
      input = document.createElement(p.type === "matrix" ? "textarea" : "input");
      if (p.type === "text") input.type = "text";
      else input.rows = p.rows || 3;
      input.className = "params-text-input";
      input.value = currentParamValues[p.key];
      // Unlike sliders (where auto-run-while-dragging is the point), typing
      // a matrix/vector is many keystrokes -- re-solving on every single one
      // is wasted work at best and, for an expensive figure (e.g. 11.1's
      // ~1-3s SDP solve), makes the whole page stutter while you type. Only
      // re-run once the field loses focus (or Enter is pressed for a
      // single-line input), not on every "input" event.
      input.addEventListener("change", () => {
        currentParamValues[p.key] = input.value;
        scheduleAutoRun();
      });
    } else if (p.type === "range2") {
      // A single visual track with two independently-draggable handles
      // (native <input type=range> only ever has one) -- two transparent-
      // track range inputs are stacked via CSS so only their thumbs are
      // interactive, with a separate rail/fill pair drawn behind them for
      // the single-track look. Backed by two plain scalar params (keyLow/
      // keyHigh), same as if they were two ordinary sliders.
      const valueSpan = document.createElement("span");
      valueSpan.className = "params-value";
      const fmtPair = () => `${currentParamValues[p.keyLow]} .. ${currentParamValues[p.keyHigh]}`;
      valueSpan.textContent = fmtPair();
      header.appendChild(valueSpan);

      const track = document.createElement("div");
      track.className = "params-range2-track";
      const rail = document.createElement("div");
      rail.className = "range2-rail";
      const fill = document.createElement("div");
      fill.className = "range2-fill";
      track.appendChild(rail);
      track.appendChild(fill);

      const lowInput = document.createElement("input");
      lowInput.type = "range";
      lowInput.min = p.min;
      lowInput.max = p.max;
      lowInput.step = p.step;
      lowInput.value = currentParamValues[p.keyLow];

      const highInput = document.createElement("input");
      highInput.type = "range";
      highInput.min = p.min;
      highInput.max = p.max;
      highInput.step = p.step;
      highInput.value = currentParamValues[p.keyHigh];

      function updateFill() {
        const lo = parseFloat(lowInput.value);
        const hi = parseFloat(highInput.value);
        const pctLo = ((lo - p.min) / (p.max - p.min)) * 100;
        const pctHi = ((hi - p.min) / (p.max - p.min)) * 100;
        fill.style.left = pctLo + "%";
        fill.style.width = Math.max(0, pctHi - pctLo) + "%";
      }
      updateFill();

      lowInput.addEventListener("input", () => {
        let lo = parseFloat(lowInput.value);
        const hi = parseFloat(highInput.value);
        if (lo > hi) {
          lo = hi;
          lowInput.value = lo;
        }
        currentParamValues[p.keyLow] = lo;
        valueSpan.textContent = fmtPair();
        updateFill();
        scheduleAutoRun();
      });
      highInput.addEventListener("input", () => {
        let hi = parseFloat(highInput.value);
        const lo = parseFloat(lowInput.value);
        if (hi < lo) {
          hi = lo;
          highInput.value = hi;
        }
        currentParamValues[p.keyHigh] = hi;
        valueSpan.textContent = fmtPair();
        updateFill();
        scheduleAutoRun();
      });

      track.appendChild(lowInput);
      track.appendChild(highInput);
      input = track;
    } else if (p.type === "select") {
      // A fixed menu of named choices (e.g. "which value field") instead of
      // a numeric range -- the figure switches behavior based on the chosen
      // option's `value`, not a continuously-varying number.
      input = document.createElement("select");
      input.className = "params-select";
      for (const opt of p.options) {
        const optionEl = document.createElement("option");
        optionEl.value = opt.value;
        optionEl.textContent = opt.label;
        if (opt.value === currentParamValues[p.key]) optionEl.selected = true;
        input.appendChild(optionEl);
      }
      input.addEventListener("change", () => {
        currentParamValues[p.key] = input.value;
        scheduleAutoRun();
      });
    } else if (p.type === "checkbox") {
      // A plain on/off toggle -- no numeric value readout, just the label
      // (above, via the shared header) and the checkbox itself (below).
      input = document.createElement("input");
      input.type = "checkbox";
      input.checked = !!currentParamValues[p.key];
      input.addEventListener("change", () => {
        currentParamValues[p.key] = input.checked;
        scheduleAutoRun();
      });
    } else if (p.log) {
      // Log-scaled range: the <input> itself is a plain linear slider over
      // an internal 0..LOG_STEPS position, mapped exponentially onto the
      // real [p.min, p.max] value (both must be > 0) -- so e.g. 0.0001 and
      // 0.1 get the same drag precision, instead of small values being
      // squeezed into an unusable sliver of a linear slider's track.
      const valueSpan = document.createElement("span");
      valueSpan.className = "params-value";
      header.appendChild(valueSpan);

      const LOG_STEPS = 1000;
      const logMin = Math.log(p.min);
      const logMax = Math.log(p.max);
      const posToValue = (pos) => Math.exp(logMin + (pos / LOG_STEPS) * (logMax - logMin));
      const valueToPos = (v) => Math.round(((Math.log(v) - logMin) / (logMax - logMin)) * LOG_STEPS);
      const fmt = (v) => Number(v.toPrecision(3)).toString();

      valueSpan.textContent = fmt(currentParamValues[p.key]);
      input = document.createElement("input");
      input.type = "range";
      input.min = 0;
      input.max = LOG_STEPS;
      input.step = 1;
      input.value = valueToPos(currentParamValues[p.key]);
      input.addEventListener("input", () => {
        const v = posToValue(parseFloat(input.value));
        currentParamValues[p.key] = v;
        valueSpan.textContent = fmt(v);
        scheduleAutoRun();
      });
    } else {
      const valueSpan = document.createElement("span");
      valueSpan.className = "params-value";
      valueSpan.textContent = currentParamValues[p.key];
      header.appendChild(valueSpan);

      input = document.createElement("input");
      input.type = "range";
      input.min = p.min;
      input.max = p.max;
      input.step = p.step;
      input.value = currentParamValues[p.key];
      input.addEventListener("input", () => {
        const v = parseFloat(input.value);
        currentParamValues[p.key] = Number.isFinite(v) ? v : p.default;
        valueSpan.textContent = currentParamValues[p.key];
        scheduleAutoRun();
      });
    }

    field.appendChild(header);
    field.appendChild(input);
    grid.appendChild(field);
  }
}

function selectEntry(entry) {
  currentEntry = entry;
  document.querySelectorAll(".figure-item").forEach((el) => el.classList.toggle("active", el.dataset.key === entry.key));
  els.figTitle.textContent = entry.title;
  els.figNote.textContent = entry.note || "";
  if (entry.method) {
    els.methodNoteText.textContent = entry.method;
    els.methodNote.style.display = "";
  } else {
    els.methodNote.style.display = "none";
  }
  els.runBtn.disabled = false;
  setStatus("未実行");
  els.outputGrid.innerHTML = "";
  els.outputHint.style.display = "";

  currentParamValues = {};
  if (entry.params) {
    entry.params.forEach((p) => {
      if (p.type === "range2") {
        currentParamValues[p.keyLow] = p.default[0];
        currentParamValues[p.keyHigh] = p.default[1];
      } else {
        currentParamValues[p.key] = p.default;
      }
    });
    renderParamsPanel(entry);
  } else {
    els.paramsPanel.style.display = "none";
    els.paramsPanel.innerHTML = "";
  }

  // Execution is instant here (unlike the Pyodide site), so running once
  // immediately on selection -- with the default parameters -- reads as
  // natural rather than as an unwanted background load.
  runCurrent();
}

function runCurrent() {
  if (!currentEntry) return;
  const entry = currentEntry;
  const fn = window.figureLib && window.figureLib[entry.key];
  if (!fn) {
    setStatus("未実装です", "error");
    return;
  }
  els.outputGrid.innerHTML = "";
  els.outputHint.style.display = "none";
  try {
    const t0 = performance.now();
    fn(els.outputGrid, currentParamValues);
    window.plotlib.finalizeResponsiveLayout();
    const ms = Math.round(performance.now() - t0);
    setStatus(`完了（${ms}ms）`, "ok");
  } catch (err) {
    console.error(err);
    setStatus("エラー: " + err.message, "error");
  }
}

els.runBtn.addEventListener("click", runCurrent);

buildSidebar();
if (FIGURES.length) selectEntry(FIGURES[0]);
