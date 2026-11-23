/* Block-by-block section: one run of the model version chosen in the model tab above (window.FS_MODEL, event
   'modelchange'), as functions with inputs, outputs, constants and state, and a chart of the hovered block.
   Data: BK_DATA (docs/blocks/traces.json, made by docs/blocks/run.py). */
(() => {
"use strict";
const DATA = JSON.parse(document.getElementById('bk-traces').textContent);
const root = document.getElementById('blocks');
const $b = s => root.querySelector(s);

const A3 = 'A Table 3', A1 = 'A Table 1a', A1b = 'A Table 1b', A2 = 'A Table 2', B1 = 'B Table 1';
const R3 = 'Randall 2019 Table 3', OURS = 'ours', U98 = 'Ursino 1998', UT1 = 'Ursino 1998 Table 1', UT2 = 'Ursino 1998 Table 2',
  UT3 = 'Ursino 1998 Table 3', CK = 'Clark 2018 code';

/* ---------- pipeline f1-f3 ---------- */
const RUN = {
  0: ['run(p_ut, w_um=1, g=1)', 'Jongen 2016', 'maternal blocks M1 and fetal blocks 1–10'],
  1: ['run(p_ut, w_um=1, g=1, mhr=2, vaso=1, autotx=1)', 'mhr=2 closed loop, vaso=1, autotx=1', 'maternal blocks M1–M8 and fetal blocks 1–10'],
  2: ['run(p_ut, w_um=1, g=1, circ=1, uterus=2)', 'circ=1 Ursino 1998, uterus=2 Clark 2018 with L and C', 'maternal blocks M1–M5 and fetal blocks 1–10'],
};
const pipe = m => [
  {k: 'toco', n: 'f1', name: 'Bin the TOCO', fn: 'toco_bins', file: 'input/ctg_signals.py',
   in: [['toco_4hz', 'CTG recording']], out: ['bins'], par: [['samples per bin', '8 (2 s)', 'ours']], st: [],
   eq: 'bins[i] = mean(toco_4hz[8i … 8i+7])'},
  {k: 'put', n: 'f2', name: 'Make uterine pressure', fn: 'uterine_pressure', file: 'input/toco_to_pressure.py',
   in: [['bins', 'f1']], out: ['p_ut'],
   par: [['a', '1 (fitted)', 'fit'], ['height', '70 mmHg', 'A §2.4'], ['tone', '15 mmHg', 'A §2.4'], ['σ', '15 s', 'ours'],
         ['peak rules', 'prominence 15, ≥80 s apart, ≥16 s wide', 'ours']], st: [],
   eq: 't_k = find_peaks(medfilt(bins))\np_ut(t) = 15 + 70a · Σ_k exp(−(t − t_k)² / 2σ²)'},
  {k: 'out', n: 'f3', name: 'Run the model', fn: 'run', file: 'model/simulate.py → core.c',
   in: [['p_ut', 'f2']], out: ['FHR', 'MHR'],
   par: [['w_um', '1 (fitted)', 'fit'], ['g', '1 (fitted)', 'fit'], ['mode', RUN[m][1], OURS]],
   st: ['all block states below, from rest'],
   eq: `o = ${RUN[m][0]}\none run, every 0.02 s: ${RUN[m][2]}\nFHR = o["fhr"],  MHR = o["mhr"]     (averaged into 2-s bins)\no also holds SpO2 and StO2 (block 11)`},
];

/* ---------- fetal blocks (the same in all versions except block 1) ---------- */
const PLACENTA_J = (pa, pv, from) => ({k: 'placenta', n: '1', name: 'Uterine blood flow (links mother and fetus)', fn: 'placenta_inflow', file: 'core.c, A(3) A(6) A(7)',
  in: [['p_ut', 'f2'], [pa, from, 1], [pv, from, 1]], out: ['q_uta', 'q_utv'],
  par: [['R_uta', '5.71 mmHg·s/ml', A1b], ['R_utv law', 'R0 0.358, p0 −2.5, p1 25', A1b], ['IVS law', 'V_max 350 ml, p0 5, p1 46.4', 'derived']],
  st: ['V_ivs'],
  eq: `p_ivs = p0 + p1·tan(π(V_ivs/V_max − ½)) + p_ut\nq_uta = (${pa} − p_ivs) / R_uta\nq_utv = (p_ivs − ${pv}) / R_utv(½(p_ivs + ${pv}) − p_ut)\ndV_ivs/dt = q_uta − q_utv`});
const PLACENTA_CK = {k: 'placenta', n: '1', name: 'Uterine network (links mother and fetus)', fn: 'uterine_network', file: 'core.c, ck_rlc_step, A(3) A(6) A(7)',
  in: [['p_ut', 'f2'], ['P_sa', 'M1', 1], ['P_ev', 'M1', 1]], out: ['q_ut,art', 'q_uta', 'q_an', 'q_utv'],
  par: [['uterine artery R, L, C', '0.257, 0.0705, 0.00335 (5 segments)', CK], ['arcuate R, L, C', '0.0462, 0.0127, 0.000603', CK],
        ['radial R, L, C', '3.08, 0.00846, 0.000101', CK], ['anastomoses R, L, C', '4.62, 0.0127, 0.000151', CK], ['R_D', '1.69 mmHg·s/ml', CK],
        ['R_spiral', '1.13 mmHg·s/ml', 'calibrated'], ['k_a', '0.7', OURS],
        ['R_utv law', 'R0 0.358, p0 −2.5, p1 25', A1b], ['IVS law', 'V_max 350 ml, p0 5, p1 46.4', 'derived']],
  st: ['q_k, p_k of the 8 elements', 'V_ivs'],
  eq: 'elements k: uterine artery ×5, arcuate, radial (in series from P_sa), anastomoses (from the radial end)\nL_k dq_k/dt = p_k−1 − p_k − R_k q_k,   C_k dp_k/dt = q_k − Σ q_out\nq_ut,art = q_0\nq_uta = (p_radial − p_ivs) / R_spiral                (spiral arteries → IVS)\nq_an = (p_av − p_out) / R_D,  p_out = max(P_ev, k_a·p_ut),  0 if p_av ≤ p_out\np_ivs = p0 + p1·tan(π(V_ivs/V_max − ½)) + p_ut\nq_utv = (p_ivs − P_ev) / R_utv(½(p_ivs + P_ev) − p_ut)\ndV_ivs/dt = q_uta − q_utv\nbackward Euler in each Ursino substep'};
const FETAL = [
  {k: 'cord', n: '2', name: 'Cord compression', fn: 'cord', file: 'core.c, B(1) B(2)',
   in: [['p_ut', 'f2'], ['p_a', '3', 1], ['p_vil', '3', 1], ['p_v', '3', 1]], out: ['R_uma', 'R_umv'],
   par: [['w_um', '1 (fitted)', 'fit'], ['p_rest', '15 mmHg', 'A §2.4'], ['artery law', 'R0 2.34, p0 −27, p1 16, R_ref 3.84', B1],
         ['vein law', 'R0 2.05, p0 −17, p1 10, R_ref 2.56', B1]], st: [],
   eq: 'p_um = w_um · (p_ut − 15)\np_tm = ½(p_vil + p_a) − p_um        (vein: p_v)\nR_uma, R_umv = max(R_ref, R0 / f²),  f = ½ + atan((p_tm − p0)/p1)/π'},
  {k: 'circ', n: '3', name: 'Fetal circulation', fn: 'fetal_circulation', file: 'core.c, A(9)–A(11)',
   in: [['R_uma', '2'], ['R_umv', '2'], ['R_cera', '7', 1], ['T', '10', 1], ['z_E', '9', 1], ['z_R', '9', 1], ['z_V', '9', 1]],
   out: ['p_a', 'p_v', 'p_vil', 'q (fetal flows)', 'V (fetal volumes)'],
   par: [['C, V0 per compartment', '', A1], ['R_tis, R_cer,v', '', A1], ['E_max, E_min, V0_max, V0_min', '', A1]],
   st: ['V_art', 'V_ven', 'V_tis', 'V_cer', 'V_vil'],
   eq: 'p_i = (V_i − V0_i) / C_i                 (V0_ven × (1 + z_V))\nq = Δp / R                               (R_tis,a × (1 + z_R))\nq_co = (V_ed − V_es) / T,  V_es = V0_max + p_a / (E_max(1 + z_E))\ndV_i/dt = Σ q_in − Σ q_out\nreturns p_a, p_v, p_vil, the flows q and the volumes V'},
  {k: 'oxygen', n: '4', name: 'Oxygen', fn: 'oxygen', file: 'core.c, A(16) A(17)',
   in: [['q_uta', '1'], ['q_utv', '1'], ['q (fetal flows)', '3'], ['V (fetal volumes)', '3']], out: ['pO2_a', 'c_a', 'pO2_tis', 'pO2_cer', 'pO2_a,m', 'pO2_tis,m'],
   par: [['D_plac', '0.082 ml O2/(s·mmHg)', A2], ['Q_tis, Q_cer', '0.336, 0.131 ml O2/s', A2], ['Hb curve α, Hb, c1, c2, β', '', A2]],
   st: ['M_i (O2 in each compartment, mother and fetus)'],
   eq: 'q_O2 = D_plac · (pO2_ivs − pO2_vil)\ndM_ivs/dt = q_uta · c_a,m − q_utv · c_ivs − q_O2     (uterine blood flow carries the O2 in)\ndM_i/dt = Σ q · c_upstream − consumption_i\nc_i = M_i / V_i\npO2_i = c⁻¹(c_i)   c(p) = αHb·g/(g + c1) + βp, g = p³ + c2·p (Newton)\nc_a = M_art / V_art,   pO2_a = c⁻¹(c_a)'},
  {k: 'chemo', n: '5', name: 'Chemoreceptor', fn: 'chemoreceptor', file: 'core.c, A(20) A(21)',
   in: [['pO2_a', '4']], out: ['r_c'],
   par: [['η1', '−6.909', 'B(3)'], ['η2', '0.52 /mmHg', A3], ['r_min, r_max', '−0.001, 1', A3], ['y_ref', '18.1 mmHg', A3], ['τ', '2 s', A3]],
   st: ['r_c'], eq: 'r_c* = r_min + (r_max − r_min) · exp(η1 · exp(η2 (pO2_a − y_ref)))\nτ · dr_c/dt = r_c* − r_c'},
  {k: 'baro', n: '6', name: 'Baroreceptor', fn: 'baroreceptor', file: 'core.c, A(20) A(21)',
   in: [['p_a', '3']], out: ['r_b'],
   par: [['η1', '−0.693', 'B(3)'], ['η2', '−0.06 /mmHg', A3], ['r_min, r_max', '−1, 1', A3], ['y_ref', '45 mmHg', A3], ['τ', '8 s', A3]],
   st: ['r_b'], eq: 'r_b* = r_min + (r_max − r_min) · exp(η1 · exp(η2 (p_a − y_ref)))\nτ · dr_b/dt = r_b* − r_b'},
  {k: 'brain', n: '7', name: 'Brain autoregulation', fn: 'brain_autoregulation', file: 'core.c, A(18) A(19)',
   in: [['c_a', '4']], out: ['R_cera'],
   par: [['R_cer,a ref', '7.24 mmHg·s/ml', A1], ['γ', '0.5', A3], ['c_ref', '0.107 ml O2/ml', A3], ['τ', '10 s', A3]],
   st: ['R_cera'], eq: 'R* = R_ref / ((1 − γ) + γ · c_ref / c_a)     (≥ R_ref/2)\nτ · dR_cera/dt = R* − R_cera'},
  {k: 'reflex', n: '8', name: 'Vagal and sympathetic', fn: 'autonomic', file: 'core.c, A(22)',
   in: [['r_b', '6'], ['r_c', '5']], out: ['e_v', 'e_s'],
   par: [['w_bv', '0.25', A3], ['w_cv', '0.75', A3], ['w_bs', '0.25', A3], ['w_cs', '0.75', A3]], st: [],
   eq: 'e_v =  w_bv · r_b + w_cv · r_c\ne_s = −w_bs · r_b + w_cs · r_c'},
  {k: 'effectors', n: '9', name: 'Effectors', fn: 'effectors', file: 'core.c, A(23)',
   in: [['e_v', '8'], ['e_s', '8']], out: ['z_Tv', 'z_Ts', 'z_E', 'z_R', 'z_V'],
   par: [['g', '1 (fitted)', 'fit'], ['Tv: D, τ, k', '0.2 s, 1.5 s, 2.4', A3], ['Ts: D, τ, k', '2 s, 2 s, 1.2', A3],
         ['E: D, τ, k', '2 s, 8 s, 1.2', A3], ['R: D, τ, k', '2 s, 6 s, 4.2', A3], ['V: D, τ, k', '5 s, 20 s, −1.2', A3]],
   st: ['x_Tv … x_V', 'delay buffers of e_v, e_s'],
   eq: 'τ_i · dx_i/dt = e(t − D_i) − x_i     e = e_v for Tv, e_s for the other four\nz_i = g · k_i · x_i      → z_Tv, z_Ts, z_E, z_R, z_V'},
  {k: 'heart', n: '10', name: 'Heart period', fn: 'heart', file: 'core.c, A(25)',
   in: [['z_Tv', '9'], ['z_Ts', '9']], out: ['T', 'FHR'], par: [['T_ref', '0.444 s (135 bpm)', A1]], st: [],
   eq: 'T = T_ref · (1 + z_Tv − z_Ts)\nFHR = 60 / T'},
  {k: 'sat', n: '11', name: 'Oxygen saturation', fn: 'add_saturations', file: 'simulate.py (run_core), hemoglobin term of A(13)',
   in: [['pO2_a', '4'], ['pO2_tis', '4'], ['pO2_cer', '4'], ['pO2_a,m', '4'], ['pO2_tis,m', '4']],
   out: ['SpO2_f', 'StO2_f_tis', 'StO2_f_cer', 'SpO2_m', 'StO2_m'],
   par: [['fetal c1, c2', '10400, 150', A2], ['maternal c1, c2', '23400, 150', A2], ['STO2_ART', '0.25', 'ours']], st: [],
   eq: 'S(p) = 100 · g / (g + c1),  g = p³ + c2·p     (fetal S_f, maternal S_m)\nSpO2_f = S_f(pO2_a)          SpO2_m = S_m(pO2_a,m)\nStO2_f_tis = 0.25 · S_f(pO2_a) + 0.75 · S_f(pO2_tis)\nStO2_f_cer = 0.25 · S_f(pO2_a) + 0.75 · S_f(pO2_cer)\nStO2_m = 0.25 · S_m(pO2_a,m) + 0.75 · S_m(pO2_tis,m)'},
];

/* ---------- maternal blocks per version ---------- */
const MCIRC_EQ = (tm, vu) => `p_a,m = (V_art,m − V0_art) / C_art\np_v,m = (V_ven,m − (V0_ven${vu})) / C_ven\nq_co,m = (V_ed − V_es) / ${tm},  V_ed = V0_min + p_v,m / E_min,  V_es = p_a,m / E_max\nq_tis = (p_a,m − p_tis,m) / R_tis,a,m\ndV_art,m/dt = q_co,m − q_tis − q_uta,   dV_ven,m/dt = q_tis,v + q_utv − q_co,m`;
const MPAR = [['maternal C, V0', 'art 3.9 ml/mmHg 721 ml, ven 308.9 / 3604, tis 6.4 / 504', A1b], ['R_tis,v', '0.04 mmHg·s/ml', A1b],
  ['maternal E_max, E_min', '1.83, 0.0571 mmHg/ml', A1b]];
const MOTHER = {
  0: [
    {k: 'mcirc', n: 'M1', name: 'Maternal circulation', fn: 'maternal_circulation', file: 'core.c, A(9)–A(11) for the mother',
     in: [['q_uta', '1'], ['q_utv', '1']], out: ['p_a,m', 'p_v,m'],
     par: MPAR.concat([['T_m', '0.75 s (80 bpm)', A1b], ['R_tis,a', '0.67 mmHg·s/ml', A1b]]),
     st: ['V_art,m', 'V_ven,m', 'V_tis,m'], eq: MCIRC_EQ('T_m', '') + '\nMHR = 60 / T_m = 80 bpm, not regulated'},
  ],
  1: [
    {k: 'autotx', n: 'M1', name: 'Autotransfusion', fn: 'autotransfusion', file: 'core.c, MA_MODE',
     in: [['p_ut', 'f2']], out: ['v_sq'],
     par: [['V_sq', '625 ml', OURS], ['p_ref', '70 mmHg', OURS], ['p_rest', '15 mmHg', 'A §2.4']], st: [],
     eq: 'x_c = max(0, p_ut − 15)\nv_sq = V_sq · x_c / p_ref'},
    {k: 'mcirc', n: 'M2', name: 'Maternal circulation', fn: 'maternal_circulation', file: 'core.c, A(9)–A(11) for the mother',
     in: [['q_uta', '1'], ['q_utv', '1'], ['v_sq', 'M1'], ['T_m', 'M7', 1], ['R_tis,a,m', 'M8', 1], ['v_vaso', 'M8', 1]],
     out: ['p_a,m', 'p_v,m'], par: MPAR, st: ['V_art,m', 'V_ven,m', 'V_tis,m'], eq: MCIRC_EQ('T_m', ' − v_vaso − v_sq')},
    {k: 'mbaro', n: 'M3', name: 'Baroreceptor', fn: 'baroreceptor', file: 'core.c, Randall 2019',
     in: [['p_a,m', 'M2']], out: ['n_b'],
     par: [['ΔSBP', '34.66 mmHg', OURS], ['P_th', '4.5 mmHg', OURS], ['A', '5', R3], ['q_w', '0.04 /mmHg', R3], ['s_w', '115 mmHg', R3],
           ['K_b', '0.1', R3], ['τ_b', '0.9 s', R3], ['B', '0.5', R3]],
     st: ['ε_b,c', 'ε_b,a'],
     eq: 'SBP = p_a,m + ΔSBP\nε_w(P) = 1 − √((1 + e^(−q_w(P − s_w))) / (A + e^(−q_w(P − s_w))))\nε_w,c = ε_w(SBP),   ε_w,a = ε_w(SBP − P_th)\nτ_b · dε_b/dt = −ε_b + K_b · ε_w          (carotid c, aortic a)\nn_b = B(ε_w,c − ε_b,c) + (1 − B)(ε_w,a − ε_b,a)'},
    {k: 'resp', n: 'M4', name: 'Respiratory parasympathetic', fn: 'respiratory', file: 'core.c, Randall 2019',
     in: [['P_th', 'constant']], out: ['T_pr'],
     par: [['K_pr', '1', R3], ['τ_pr', '6 s', R3], ['q_pr', '1 /mmHg', R3], ['s_pr', '4.5 mmHg', R3]], st: ['T_pr'],
     eq: 'G_pr = 1 / (1 + e^(q_pr(P_th − s_pr)))\nτ_pr · dT_pr/dt = −T_pr + K_pr · G_pr\nP_th is fixed at 4.5 mmHg here, so T_pr stays at 0.5'},
    {k: 'para', n: 'M5', name: 'Parasympathetic', fn: 'parasympathetic', file: 'core.c, Randall 2019',
     in: [['n_b', 'M3']], out: ['T_pb'],
     par: [['K_pb', '5', R3], ['τ_pb', '1.8 s', R3], ['q_pb', '10', R3], ['s_pb', '0.538', R3]], st: ['T_pb'],
     eq: 'G_pb = 1 / (1 + e^(−q_pb(n_b − s_pb)))\nτ_pb · dT_pb/dt = −T_pb + K_pb · G_pb'},
    {k: 'symp', n: 'M6', name: 'Sympathetic', fn: 'sympathetic', file: 'core.c, Randall 2019',
     in: [['n_b', 'M3']], out: ['T_s'],
     par: [['K_s', '5', R3], ['τ_s', '10 s', R3], ['q_s', '10', R3], ['s_s', '0.0544', R3], ['D_s', '3 s', R3]],
     st: ['T_s', 'delay buffer of T_s'],
     eq: 'G_s = 1 / (1 + e^(q_s(n_b − s_s)))\nτ_s · dT_s/dt = −T_s(t − D_s) + K_s · G_s'},
    {k: 'mheart', n: 'M7', name: 'Heart rate', fn: 'heart_rate', file: 'core.c, Randall 2019',
     in: [['T_pb', 'M5'], ['T_pr', 'M4'], ['T_s', 'M6']], out: ['MHR', 'T_m'],
     par: [['H_I', '100.9 bpm', R3], ['H_pb', '0.425', R3], ['H_pr', '0.198', R3], ['H_s', '0.171', R3], ['τ_H', '0.5 s', R3]],
     st: ['MHR'],
     eq: 'H̃ = H_I · (1 − H_pb·T_pb + H_pr·T_pr + H_s·T_s)\nτ_H · dMHR/dt = −MHR + H̃\nT_m = 60 / MHR'},
    {k: 'vaso', n: 'M8', name: 'Vasomotor', fn: 'vasomotor', file: 'core.c, MV_MODE',
     in: [['T_s', 'M6']], out: ['R_tis,a,m', 'v_vaso'],
     par: [['T_s0', '0.2', OURS], ['k_R', '0.671', OURS], ['k_V', '3850 ml', OURS], ['r_min', '0.3', OURS],
           ['τ_R', '6 s', U98], ['τ_V', '20 s', U98], ['R_tis,a', '0.67 mmHg·s/ml', A1b]],
     st: ['x_R', 'x_V'],
     eq: 'τ_R · dx_R/dt = (T_s − T_s0) − x_R\nτ_V · dx_V/dt = (T_s − T_s0) − x_V\nR_tis,a,m = R_tis,a · max(r_min, 1 + k_R · x_R)\nv_vaso = k_V · x_V'},
  ],
  2: [
    {k: 'uhv', n: 'M1', name: 'Heart and vessels', fn: 'ursino_circulation', file: 'ursino1998.c, U(1)–U(29)',
     in: [['T_m', 'M5', 1], ['E_max,lv', 'M4', 1], ['E_max,rv', 'M4', 1], ['R_sp', 'M4', 1], ['R_ep', 'M4', 1], ['V_u,sv', 'M4', 1], ['V_u,ev', 'M4', 1],
          ['q_ut,art', '1'], ['q_utv', '1'], ['q_an', '1']],
     out: ['P_sa', 'P_ev'],
     par: [['C_sa, R_sa, L_sa', '0.28 ml/mmHg, 0.06 mmHg·s/ml, 0.00022 mmHg·s²/ml', UT1], ['C_sp, C_ep', '2.05, 1.67 ml/mmHg', UT1],
           ['C_sv, C_ev', '61.11, 50 ml/mmHg', UT1], ['R_sv, R_ev', '0.038, 0.016 mmHg·s/ml', UT1], ['R_sp,0, R_ep,0', '2.73, 0.855 mmHg·s/ml', 'calibrated'],
           ['V_t', '5300 ml', UT1], ['left heart P0, k_E, k_R', '1.5 mmHg, 0.014 /ml, 0.000375 s/ml', UT2], ['T_sys,0, k_sys', '0.5 s, 0.075 s²', UT2],
           ['n_sub', '100 RK4 substeps per 0.02 s', OURS]],
     st: ['P_sa, F_sa, P_sp, P_ev', 'P_la, V_lv, P_ra, V_rv', 'P_pa, F_pa, P_pp, P_pv', 'ξ (heart phase)', 'V_t'],
     eq: 'φ = sin²(π·T_m·u / T_sys) for u ≤ T_sys/T_m, else 0;   u = frac(ξ),  dξ/dt = 1/T_m,  T_sys = T_sys,0 − k_sys/T_m\nP_max,lv = φ·E_max,lv·(V_lv − V_u,lv) + (1 − φ)·P0,lv·(e^(k_E,lv·V_lv) − 1)\nF_o,l = (P_max,lv − P_sa) / (k_R,lv·P_max,lv)  if P_max,lv > P_sa, else 0       (right heart likewise)\nC_sa · dP_sa/dt = F_o,l − F_sa − q_ut,art\nL_sa · dF_sa/dt = P_sa − P_sp − R_sa·F_sa\n(C_sp + C_ep) · dP_sp/dt = F_sa − (P_sp − P_sv)/R_sp − (P_sp − P_ev)/R_ep\nC_ev · dP_ev/dt = (P_sp − P_ev)/R_ep − (P_ev − P_ra)/R_ev + q_utv + q_an − dV_u,ev/dt\nP_sv from the blood volume:  dV_t/dt = q_utv + q_an − q_ut,art'},
    {k: 'ucs', n: 'M2', name: 'Carotid baroreceptor', fn: 'carotid_afferent', file: 'ursino1998.c, U(30) U(31)',
     in: [['P_sa', 'M1']], out: ['f_cs'],
     par: [['τ_p, τ_z', '2.076, 6.37 s', UT3], ['P_n', '92 mmHg', UT3], ['k_a', '11.758 mmHg', UT3], ['f_min, f_max', '2.52, 47.78 spikes/s', UT3]],
     st: ['P̃'],
     eq: 'τ_p · dP̃/dt = P_sa + τ_z · dP_sa/dt − P̃\nf_cs = (f_min + f_max · e^((P̃ − P_n)/k_a)) / (1 + e^((P̃ − P_n)/k_a))'},
    {k: 'ueff', n: 'M3', name: 'Sympathetic and vagal', fn: 'efferent', file: 'ursino1998.c, U(33) U(34)',
     in: [['f_cs', 'M2']], out: ['f_es', 'f_ev'],
     par: [['f_es,∞, f_es,0, k_es', '2.1, 16.11 spikes/s, 0.0675 s', UT3], ['f_ev,0, f_ev,∞', '3.2, 6.3 spikes/s', UT3], ['f_cs,0, k_ev', '25, 7.06 spikes/s', UT3]],
     st: [],
     eq: 'f_es = f_es,∞ + (f_es,0 − f_es,∞) · e^(−k_es · f_cs)\nf_ev = (f_ev,0 + f_ev,∞ · e^((f_cs − f_cs,0)/k_ev)) / (1 + e^((f_cs − f_cs,0)/k_ev))'},
    {k: 'ueffec', n: 'M4', name: 'Effectors', fn: 'effectors', file: 'ursino1998.c, U(35)–U(41)',
     in: [['f_es', 'M3'], ['f_ev', 'M3']], out: ['E_max,lv', 'E_max,rv', 'R_sp', 'R_ep', 'V_u,sv', 'V_u,ev', 'ΔT_s', 'ΔT_v'],
     par: [['E_max,lv: G, τ, D', '0.475, 8 s, 2 s', UT3], ['E_max,rv: G, τ, D', '0.282, 8 s, 2 s', UT3], ['R_sp: G, τ, D', '0.695, 6 s, 2 s', UT3],
           ['R_ep: G, τ, D', '0.53, 6 s, 2 s', UT3], ['V_u,sv: G, τ, D', '−265.4 ml, 20 s, 5 s', UT3], ['V_u,ev: G, τ, D', '−132.5 ml, 20 s, 5 s', UT3],
           ['T_s: G, τ, D', '−0.13 s, 2 s, 2 s', UT3], ['T_v: G, τ, D', '0.09 s, 1.5 s, 0.2 s', UT3], ['f_es,min', '2.66 spikes/s', UT3]],
     st: ['Δθ for the eight effectors', 'delay buffers of f_es, f_ev'],
     eq: 'σ_θ = G_θ · ln(f_es(t − D_θ) − f_es,min + 1)      (0 below f_es,min)\nτ_θ · dΔθ/dt = σ_θ − Δθ       θ = E_max,lv, E_max,rv, R_sp, R_ep, V_u,sv, V_u,ev, T_s\nτ_Tv · dΔT_v/dt = G_Tv · f_ev(t − D_Tv) − ΔT_v\nθ = θ_0 + Δθ'},
    {k: 'uhp', n: 'M5', name: 'Heart period', fn: 'heart_period', file: 'ursino1998.c, U(42)',
     in: [['ΔT_s', 'M4'], ['ΔT_v', 'M4']], out: ['T_m', 'MHR'], par: [['T_0', '0.58 s', UT3]], st: [],
     eq: 'T_m = T_0 + ΔT_s + ΔT_v\nMHR = 60 / T_m'},
  ],
};
const KEYS = {0: ['placenta', 'oxygen', 'sat', 'heart', 'mcirc'], 1: ['placenta', 'oxygen', 'sat', 'heart', 'mcirc'], 2: ['placenta', 'oxygen', 'sat', 'heart', 'uhv']};

/* ---------- diagram: mother band per version, fetal band shared ---------- */
const PUT_F = `<polyline class="w" points="130,405 150,405 150,882 178,882" marker-end="url(#bk-fa)"/>
  <line class="w" x1="150" y1="635" x2="168" y2="635" marker-end="url(#bk-fa)"/><text class="wl" x="156" y="520">p_ut</text>
  <g class="nd src" data-k="put" tabindex="0"><rect x="20" y="370" width="110" height="70" rx="6"/>
    <text class="nt" x="75" y="400" text-anchor="middle">p_ut(t)</text><text class="nf" x="75" y="420" text-anchor="middle">from f2</text></g>`;
const LINK = (up1, up2, dn1, dn2, up3) => `<polyline class="w" points="220,600 220,285 440,285 440,242" marker-end="url(#bk-fa)"/>
  ${up3 ? `<text class="wl" x="214" y="386" text-anchor="end">${up1}</text><text class="wl" x="214" y="400" text-anchor="end">${up2}</text><text class="wl" x="214" y="414" text-anchor="end">${up3}</text>`
    : `<text class="wl" x="226" y="400">${up1}</text><text class="wl" x="226" y="414">${up2}</text>`}
  <polyline class="fbw" points="480,240 480,320 290,320 290,598" marker-end="url(#bk-fb)"/>
  <text class="fbl" x="296" y="400">${dn1}</text><text class="fbl" x="296" y="414">${dn2}</text>`;
const ND = (k, x, y, w, t, f, key) => `<g class="nd${key ? ' bkkey' : ''}" data-k="${k}" tabindex="0"><rect x="${x}" y="${y}" width="${w}" height="70" rx="6"/>
  <text class="nt" x="${x + 12}" y="${y + 28}">${t}</text><text class="nf" x="${x + 12}" y="${y + 50}">${f}</text></g>`;
const OUTN = (k, x, y, t) => `<g class="nd" data-k="${k}" tabindex="0"><rect x="${x}" y="${y}" width="${t.length * 10 + 14}" height="40" fill="transparent" stroke="none"/>
  <text class="outl" x="${x + 4}" y="${y + 25}">${t}</text></g>`;
const MBAND = {
  0: PUT_F + LINK('q_uta,', 'q_utv', 'p_a,m,', 'p_v,m') +
    `<line class="w" x1="570" y1="205" x2="1280" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="925" y="196" text-anchor="middle">T_m fixed</text>` +
    ND('mcirc', 400, 170, 170, 'M1 · Maternal circ.', 'maternal_circulation', 1) + OUTN('out', 1282, 185, 'MHR'),
  1: `<polyline class="w" points="130,405 150,405 150,75 178,75" marker-end="url(#bk-fa)"/>
  <polyline class="w" points="150,405 150,882 178,882" marker-end="url(#bk-fa)"/>
  <line class="w" x1="150" y1="635" x2="168" y2="635" marker-end="url(#bk-fa)"/><text class="wl" x="156" y="250">p_ut</text>
  <polyline class="w" points="340,75 372,75 372,205 398,205" marker-end="url(#bk-fa)"/><text class="wl" x="378" y="140">v_sq</text>
  <line class="w" x1="570" y1="205" x2="628" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="599" y="196" text-anchor="middle">p_a,m</text>
  <line class="w" x1="790" y1="205" x2="848" y2="205" marker-end="url(#bk-fa)"/>
  <polyline class="w" points="820,205 820,335 848,335" marker-end="url(#bk-fa)"/><text class="wl" x="819" y="196" text-anchor="middle">n_b</text>
  <polyline class="w" points="1010,75 1045,75 1045,190 1068,190" marker-end="url(#bk-fa)"/><text class="wl" x="1051" y="130">T_pr</text>
  <line class="w" x1="1010" y1="205" x2="1068" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="1027" y="199" text-anchor="middle">T_pb</text>
  <polyline class="w" points="1010,322 1045,322 1045,220 1068,220" marker-end="url(#bk-fa)"/><text class="wl" x="1051" y="280">T_s</text>
  <line class="w" x1="1010" y1="350" x2="1068" y2="350" marker-end="url(#bk-fa)"/><text class="wl" x="1039" y="366" text-anchor="middle">T_s</text>
  <line class="w" x1="1230" y1="205" x2="1280" y2="205" marker-end="url(#bk-fa)"/>` + LINK('q_uta,', 'q_utv', 'p_a,m,', 'p_v,m') +
  `<polyline class="fbw" points="1150,170 1150,24 520,24 520,168" marker-end="url(#bk-fb)"/><text class="fbl" x="835" y="16" text-anchor="middle">T_m</text>
  <polyline class="fbw" points="1150,370 1150,420 540,420 540,242" marker-end="url(#bk-fb)"/><text class="fbl" x="845" y="414" text-anchor="middle">R_tis,a,m, v_vaso</text>
  <g class="nd src" data-k="put" tabindex="0"><rect x="20" y="370" width="110" height="70" rx="6"/>
    <text class="nt" x="75" y="400" text-anchor="middle">p_ut(t)</text><text class="nf" x="75" y="420" text-anchor="middle">from f2</text></g>` +
    ND('autotx', 180, 40, 160, 'M1 · Autotransfusion', 'autotransfusion') + ND('mcirc', 400, 170, 170, 'M2 · Maternal circ.', 'maternal_circulation', 1) +
    ND('mbaro', 630, 170, 160, 'M3 · Baroreceptor', 'baroreceptor') + ND('resp', 850, 40, 160, 'M4 · Respiratory', 'P_th = 4.5 (fixed)') +
    ND('para', 850, 170, 160, 'M5 · Parasympathetic', 'parasympathetic') + ND('symp', 850, 300, 160, 'M6 · Sympathetic', 'sympathetic') +
    ND('mheart', 1070, 170, 160, 'M7 · Heart rate', 'heart_rate') + ND('vaso', 1070, 300, 160, 'M8 · Vasomotor', 'vasomotor') + OUTN('out', 1282, 185, 'MHR'),
  2: PUT_F + LINK('q_ut,art,', 'q_utv,', 'P_sa,', 'P_ev', 'q_an') +
    `<line class="w" x1="570" y1="205" x2="628" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="599" y="196" text-anchor="middle">P_sa</text>
  <line class="w" x1="790" y1="205" x2="848" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="819" y="196" text-anchor="middle">f_cs</text>
  <polyline class="w" points="1010,220 1040,220 1040,335 1068,335" marker-end="url(#bk-fa)"/><text class="wl" x="1046" y="290">f_es, f_ev</text>
  <line class="w" x1="1150" y1="300" x2="1150" y2="242" marker-end="url(#bk-fa)"/><text class="wl" x="1158" y="276">ΔT_s, ΔT_v</text>
  <line class="w" x1="1230" y1="205" x2="1280" y2="205" marker-end="url(#bk-fa)"/>
  <polyline class="fbw" points="1150,170 1150,24 520,24 520,168" marker-end="url(#bk-fb)"/><text class="fbl" x="835" y="16" text-anchor="middle">T_m</text>
  <polyline class="fbw" points="1150,370 1150,420 540,420 540,242" marker-end="url(#bk-fb)"/><text class="fbl" x="845" y="414" text-anchor="middle">E_max,lv, E_max,rv, R_sp, R_ep, V_u,sv, V_u,ev</text>` +
    ND('uhv', 400, 170, 170, 'M1 · Heart, vessels', 'ursino_circulation', 1) + ND('ucs', 630, 170, 160, 'M2 · Carotid sinus', 'carotid_afferent') +
    ND('ueff', 850, 170, 160, 'M3 · Symp., vagal', 'efferent') + ND('uhp', 1070, 170, 160, 'M5 · Heart period', 'heart_period') +
    ND('ueffec', 1070, 300, 160, 'M4 · Effectors', 'effectors') + OUTN('out', 1282, 185, 'MHR'),
};
const FBAND = m => `<g transform="translate(0,560)">
  <line class="w" x1="340" y1="75" x2="398" y2="75" marker-end="url(#bk-fa)"/>
  <text class="wl" x="369" y="64" text-anchor="middle">q_uta,</text><text class="wl" x="369" y="92" text-anchor="middle">q_utv</text>
  <line class="w" x1="340" y1="322" x2="398" y2="322" marker-end="url(#bk-fa)"/><text class="wl" x="369" y="311" text-anchor="middle">R_um</text>
  <line class="w" x1="485" y1="300" x2="485" y2="112" marker-end="url(#bk-fa)"/><text class="wl" x="492" y="205">q, V</text>
  <line class="w" x1="570" y1="65" x2="628" y2="65" marker-end="url(#bk-fa)"/><text class="wl" x="599" y="56" text-anchor="middle">pO2_a</text>
  <polyline class="w" points="570,95 600,95 600,205 628,205" marker-end="url(#bk-fa)"/><text class="wl" x="606" y="150">c_a</text>
  <line class="w" x1="570" y1="335" x2="628" y2="335" marker-end="url(#bk-fa)"/><text class="wl" x="599" y="326" text-anchor="middle">p_a</text>
  <polyline class="w" points="790,75 820,75 820,195 848,195" marker-end="url(#bk-fa)"/><text class="wl" x="826" y="140">r_c</text>
  <polyline class="w" points="790,335 820,335 820,215 848,215" marker-end="url(#bk-fa)"/><text class="wl" x="826" y="280">r_b</text>
  <line class="w" x1="1010" y1="205" x2="1068" y2="205" marker-end="url(#bk-fa)"/><text class="wl" x="1039" y="196" text-anchor="middle">e_v, e_s</text>
  <line class="w" x1="1150" y1="170" x2="1150" y2="112" marker-end="url(#bk-fa)"/><text class="wl" x="1158" y="145">z_Tv, z_Ts</text>
  <line class="w" x1="1230" y1="65" x2="1280" y2="65" marker-end="url(#bk-fa)"/>
  <polyline class="w" points="485,40 485,-40 628,-40" marker-end="url(#bk-fa)"/><text class="wl" x="556" y="-49" text-anchor="middle">pO2_i</text>
  <line class="ow" x1="790" y1="-40" x2="848" y2="-40" marker-end="url(#bk-fb)"/>
  <line class="fbw" x1="400" y1="352" x2="342" y2="352" marker-end="url(#bk-fb)"/><text class="fbl" x="370" y="392" text-anchor="middle">p_a, p_vil, p_v</text>
  <polyline class="fbw" points="710,240 710,272 530,272 530,298" marker-end="url(#bk-fb)"/><text class="fbl" x="620" y="265" text-anchor="middle">R_cera</text>
  <polyline class="fbw" points="1150,240 1150,420 510,420 510,372" marker-end="url(#bk-fb)"/><text class="fbl" x="830" y="414" text-anchor="middle">z_E, z_R, z_V</text>
  <polyline class="fbw" points="1230,95 1260,95 1260,450 470,450 470,372" marker-end="url(#bk-fb)"/><text class="fbl" x="830" y="444" text-anchor="middle">T</text>` +
  (m === 2 ? ND('placenta', 170, 40, 170, '1 · Uterine network', 'uterine_network', 1) : ND('placenta', 170, 40, 170, '1 · Uterine blood flow', 'placenta_inflow', 1)) +
  ND('cord', 180, 300, 160, '2 · Cord', 'cord') + ND('circ', 400, 300, 170, '3 · Fetal circulation', 'fetal_circulation') +
  ND('oxygen', 400, 40, 170, '4 · Oxygen', 'oxygen', 1) + ND('chemo', 630, 40, 160, '5 · Chemoreceptor', 'chemoreceptor') +
  ND('baro', 630, 300, 160, '6 · Baroreceptor', 'baroreceptor') + ND('brain', 630, 170, 160, '7 · Brain autoreg.', 'brain_autoregulation') +
  ND('reflex', 850, 170, 160, '8 · Vagal, sympathetic', 'autonomic') + ND('effectors', 1070, 170, 160, '9 · Effectors', 'effectors') +
  ND('heart', 1070, 40, 160, '10 · Heart period', 'heart', 1) + ND('sat', 630, -75, 160, '11 · O2 saturation', 'add_saturations', 1) +
  OUTN('sat', 852, -60, 'SpO2, StO2') + OUTN('out', 1282, 45, 'FHR') + '</g>';
const diagram = m => `<svg class="bkflow" viewBox="0 0 1340 1030" role="img" aria-label="Maternal blocks (top) and fetal blocks (bottom) of one model step, linked by block 1">
  <defs><marker id="bk-fa" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" class="ah"/></marker>
  <marker id="bk-fb" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" class="ahf"/></marker></defs>
  <line class="sep" x1="170" y1="458" x2="1330" y2="458"/><text class="band" x="20" y="60">MOTHER</text><text class="band" x="20" y="620">FETUS</text>
  ${MBAND[m]}${FBAND(m)}</svg>`;

/* ---------- charts ---------- */
const range = (a, d = 1) => `${Math.min(...a).toFixed(d)} to ${Math.max(...a).toFixed(d)}`;
const RECOMP = 'Recorded by the instrumented copy of the model that docs/blocks/run.py builds (output lines only).';
function charts(m, D) {
  const T1 = D.t, TT = D.t.slice(0, D.toco.length);
  const C = {
    toco: {title: 'f1 · TOCO in, bins out', unit: 'TOCO units', s: [['toco', TT, 'toco_4hz (shown every 1 s)'], ['bins', D.tb, 'bins']],
      note: 'Check here that the contractions you expect are present and the dotted peak lines sit on them.'},
    put: {title: 'f2 · p_ut out', unit: 'mmHg', s: [['p_ut', T1, 'p_ut']], note: 'Every contraction gets the same height; only peak times come from TOCO.'},
    out: {title: 'f3 · FHR and MHR out', unit: 'bpm', s: [['fhr', T1, 'FHR'], ['mhr', T1, 'MHR']],
      note: `Both from the same run. FHR rests at ${D.fhr[0].toFixed(1)} bpm and falls to ${Math.min(...D.fhr).toFixed(1)}; ` +
        (m === 0 ? 'MHR is fixed at 80 bpm.' : `MHR rests at ${D.mhr[0].toFixed(1)} and ranges from ${range(D.mhr)}.`)},
    placenta: m === 2 ? {title: '1 · Uterine network', unit: 'ml/min', s: [['q_ut_art', T1, 'q_ut,art'], ['q_uta', T1, 'q_uta'], ['q_anast', T1, 'q_an']],
        note: 'Beat means. q_utv is not plotted.'}
      : {title: '1 · Uterine blood flow', unit: 'ml/min', s: [['q_uta', T1, 'q_uta'], ['q_utv', T1, 'q_utv']]},
    cord: {title: '2 · Cord compression', unit: 'mmHg·s/ml', s: [['r_uma', T1, 'R_uma']], note: 'R_umv is not plotted (needs p_v, which the model does not record).'},
    circ: {title: '3 · Fetal circulation', unit: 'mmHg', s: [['ptm_a', T1, 'p_a'], ['p_vil', T1, 'p_vil']]},
    oxygen: {title: '4 · Oxygen', unit: 'mmHg', s: [['po2_a', T1, 'pO2_a']]},
    chemo: {title: '5 · Chemoreceptor', unit: '−1 to 1', s: [['r_c', T1, 'r_c']]},
    baro: {title: '6 · Baroreceptor', unit: '−1 to 1', s: [['r_b', T1, 'r_b']]},
    brain: {title: '7 · Brain autoregulation', unit: 'ml/min', s: [['q_cera', T1, 'brain arterial flow']], note: 'R_cera itself is not recorded; the brain flow it sets is shown instead.'},
    reflex: {title: '8 · Vagal and sympathetic', unit: '', s: [['e_v', T1, 'e_v'], ['e_s', T1, 'e_s']]},
    effectors: {title: '9 · Effectors', unit: '', s: [['z_tv', T1, 'z_Tv'], ['z_ts', T1, 'z_Ts']], note: 'Only z_Tv and z_Ts are recorded; z_E, z_R, z_V are not in the model output.'},
    heart: {title: '10 · Fetal heart rate', unit: 'bpm', s: [['fhr', T1, 'FHR (1 s)'], ['fhr_bins', D.tb, 'FHR (2-s bins)']], note: 'FHR = 60 / T.'},
    sat: {title: '11 · Oxygen saturation (fetus)', unit: '%', s: [['spo2_f', T1, 'SpO2_f'], ['sto2_f_tis', T1, 'StO2_f_tis'], ['sto2_f_cer', T1, 'StO2_f_cer']],
      note: `Maternal values are not plotted: SpO2_m ${range(D.spo2_m)} %, StO2_m ${range(D.sto2_m)} %.`},
    mcirc: {title: `${m === 0 ? 'M1' : 'M2'} · Maternal circulation`, unit: 'mmHg', s: [['pa_m', T1, 'p_a,m']],
      note: `p_v,m (${range(D.pv_m)} mmHg) is on a different scale and not plotted.`},
  };
  if (m === 1) Object.assign(C, {
    autotx: {title: 'M1 · Autotransfusion', unit: 'ml', s: [['v_sq', T1, 'v_sq']]},
    mbaro: {title: 'M3 · Baroreceptor', unit: '', s: [['ewc', T1, 'ε_w,c (carotid wall strain)'], ['nb', T1, 'n_b']], note: RECOMP},
    resp: {title: 'M4 · Respiratory parasympathetic', unit: '', s: [['tpr', T1, 'T_pr']], note: 'Flat: there is no breathing input. ' + RECOMP},
    para: {title: 'M5 · Parasympathetic', unit: '', s: [['tpb', T1, 'T_pb']], note: RECOMP},
    symp: {title: 'M6 · Sympathetic', unit: '', s: [['ts_m', T1, 'T_s']]},
    mheart: {title: 'M7 · Maternal heart rate', unit: 'bpm', s: [['mhr', T1, 'MHR (1 s)'], ['mhr_bins', D.tb, 'MHR (2-s bins)']]},
    vaso: {title: 'M8 · Vasomotor', unit: 'mmHg·s/ml', s: [['rtisa_m', T1, 'R_tis,a,m']], note: `v_vaso (${range(D.v_vaso)} ml) is on a different scale and not plotted.`},
  });
  if (m === 2) Object.assign(C, {
    uhv: {title: 'M1 · Heart and vessels', unit: 'mmHg', s: [['sbp', T1, 'SBP'], ['pa_m', T1, 'P_sa'], ['dbp', T1, 'DBP']],
      note: `Systolic, mean and diastolic: maximum, mean and minimum of P_sa over each beat. P_ev (${range(D.pv_m)} mmHg) is not plotted.`},
    ucs: {title: 'M2 · Carotid baroreceptor', unit: 'spikes/s', s: [['u_fcs', T1, 'f_cs']], note: 'Beat mean. ' + RECOMP},
    ueff: {title: 'M3 · Sympathetic and vagal', unit: 'spikes/s', s: [['u_fes', T1, 'f_es'], ['u_fev', T1, 'f_ev']], note: 'Beat means. ' + RECOMP},
    ueffec: {title: 'M4 · Effectors', unit: 'mmHg·s/ml', s: [['u_rsp', T1, 'R_sp'], ['u_rep', T1, 'R_ep']],
      note: `Not plotted: E_max,lv ${range(D.u_elv, 3)} mmHg/ml, V_u,sv ${range(D.u_vusv, 0)} ml, V_u,ev ${range(D.u_vuev, 0)} ml. ` + RECOMP},
    uhp: {title: 'M5 · Maternal heart rate', unit: 'bpm', s: [['mhr', T1, 'MHR (1 s)'], ['mhr_bins', D.tb, 'MHR (2-s bins)']]},
  });
  return C;
}

/* ---------- glossary tooltips: [full name, meaning] ---------- */
const G = {
  toco_4hz: ['Raw TOCO, 4 Hz', 'Tocodynamometer reading of uterine activity from the CTG belt, arbitrary units, 4 samples per second.'],
  bins: ['TOCO in 2-s bins', 'Mean of each 8 consecutive TOCO samples.'],
  p_ut: ['Uterine pressure (mmHg)', 'Pressure inside the uterus. 15 mmHg at rest, rising with each contraction. The model\'s only driving input.'],
  FHR: ['Fetal heart rate (bpm)', 'Beats per minute, 60 / T.'],
  MHR: ['Maternal heart rate (bpm)', 'Beats per minute, 60 / T_m.'],
  q_uta: ['Spiral-artery flow into the IVS (ml/s)', 'Maternal blood entering the intervillous space, the pool that bathes the placenta. Drops during a contraction.'],
  q_utv: ['Uterine vein flow (ml/s)', 'Maternal blood leaving the intervillous space, back into the maternal veins.'],
  'q_ut,art': ['Uterine artery flow (ml/s)', 'Flow into the uterine artery from the systemic arteries: the IVS inflow plus the anastomotic bypass, both sides.'],
  q_an: ['Anastomotic bypass flow (ml/s)', 'Flow from the radial arteries through the arterio-venous anastomoses into the maternal veins, past the IVS.'],
  V_ivs: ['Intervillous space volume (ml)', 'Maternal blood held in the placental pool.'],
  maternal: ['Maternal blood volumes', 'They set maternal arterial and venous pressure.'],
  v_sq: ['Squeezed-out blood volume (ml)', 'Blood the contracting uterine muscle pushes into the maternal veins (autotransfusion). Up to 625 ml at a 70 mmHg contraction.'],
  'p_a,m': ['Maternal arterial pressure (mmHg)', 'Mean pressure in the mother\'s arteries; one compartment, beat-averaged. About 80 mmHg at rest. It drives uterine artery flow.'],
  'p_v,m': ['Maternal venous pressure (mmHg)', 'Pressure in the mother\'s veins; sets how much the heart fills each beat.'],
  P_sa: ['Systemic arterial pressure (mmHg)', 'Ursino 1998: pressure in the large systemic arteries, pulsatile. The uterine network leaves from here.'],
  P_ev: ['Extrasplanchnic venous pressure (mmHg)', 'Ursino 1998: the venous pool the uterine veins and the bypass drain into.'],
  SBP: ['Systolic pressure (mmHg)', 'Highest P_sa of the beat.'], DBP: ['Diastolic pressure (mmHg)', 'Lowest P_sa of the beat.'],
  T_m: ['Maternal heart period (s)', 'Time per beat = 60 / MHR. Closes the loop: a slower heart pumps less, lowering pressure.'],
  'R_tis,a,m': ['Maternal body artery resistance', 'Resistance of the mother\'s body arteries; the vasomotor reflex raises or lowers it.'],
  v_vaso: ['Venous volume moved by the vasomotor reflex (ml)', 'Blood shifted from unstressed to stressed venous volume (venoconstriction).'],
  'V_art,m': ['Maternal arterial volume (ml)', ''], 'V_ven,m': ['Maternal venous volume (ml)', ''], 'V_tis,m': ['Maternal tissue blood volume (ml)', ''],
  n_b: ['Baroreceptor afferent signal', 'Combined carotid and aortic stretch signal. Rises with blood pressure.'],
  'ε_b,c': ['Carotid baroreceptor strain state', 'Slowly adapting part of the carotid wall strain.'],
  'ε_b,a': ['Aortic baroreceptor strain state', 'Slowly adapting part of the aortic wall strain.'],
  'ε_w,c': ['Carotid wall strain', ''],
  P_th: ['Thoracic pressure (mmHg)', 'Pressure in the chest. Fixed at 4.5 mmHg here (no breathing or pushing input).'],
  T_pr: ['Respiratory parasympathetic tone', 'Vagal tone driven by breathing. Constant here.'],
  T_pb: ['Baroreflex parasympathetic tone', 'Vagal tone driven by the baroreceptors. Rises with pressure and slows the heart.'],
  T_s: ['Sympathetic tone', 'Randall 2019: sympathetic drive. Falls when pressure rises.'],
  x_R: ['Resistance effector state', 'Smoothed sympathetic excess driving body artery resistance.'],
  x_V: ['Venous effector state', 'Smoothed sympathetic excess driving venous volume.'],
  f_cs: ['Carotid sinus afferent firing (spikes/s)', 'Rises with carotid pressure and its rate of rise.'],
  'P̃': ['Filtered carotid pressure (mmHg)', 'P_sa passed through the receptor\'s lead-lag filter.'],
  f_es: ['Sympathetic efferent firing (spikes/s)', 'Falls as f_cs rises.'], f_ev: ['Vagal efferent firing (spikes/s)', 'Rises with f_cs.'],
  R_sp: ['Splanchnic peripheral resistance (mmHg·s/ml)', 'Raised by sympathetic firing.'], R_ep: ['Extrasplanchnic peripheral resistance (mmHg·s/ml)', 'Raised by sympathetic firing.'],
  'V_u,sv': ['Splanchnic venous unstressed volume (ml)', 'Lowered by sympathetic firing (venoconstriction).'], 'V_u,ev': ['Extrasplanchnic venous unstressed volume (ml)', 'Lowered by sympathetic firing.'],
  'E_max,lv': ['Left ventricle end-systolic elastance (mmHg/ml)', 'Contractility, raised by sympathetic firing.'], 'E_max,rv': ['Right ventricle end-systolic elastance (mmHg/ml)', ''],
  'ΔT_s': ['Sympathetic change of heart period (s)', 'Negative: speeds the heart.'], 'ΔT_v': ['Vagal change of heart period (s)', 'Positive: slows the heart.'],
  'Δθ': ['Effector changes', 'One first-order state per effector.'],
  'ξ': ['Heart phase', 'Counts beats; its fractional part is the position inside the current beat.'],
  'P_sa, F_sa, P_sp, P_ev': ['Systemic states', 'Arterial pressure and flow, peripheral and extrasplanchnic venous pressure.'],
  'P_la, V_lv, P_ra, V_rv': ['Heart states', 'Atrial pressures and ventricular volumes.'],
  'P_pa, F_pa, P_pp, P_pv': ['Pulmonary states', 'Pulmonary arterial pressure and flow, peripheral and venous pressure.'],
  V_t: ['Total blood volume (ml)', 'Changes only with the intervillous volume.'],
  'q_k, p_k of the 8 elements': ['Network states', 'Flow through and pressure at the end of each of the 8 network elements.'],
  R_spiral: ['Spiral-artery resistance', 'Calibrated so that resting IVS inflow and pressure are those of the single-R_uta model.'],
  R_D: ['Anastomotic vein resistance', 'Drains the bypass into the maternal veins.'], k_a: ['Bypass collapse factor', 'The bypass veins drain against max(P_ev, k_a·p_ut); ours, not in Clark 2018.'],
  R_um: ['Umbilical resistances', 'Resistance of the umbilical artery and vein; rises steeply when the cord is squeezed.'],
  R_uma: ['Umbilical artery resistance (mmHg·s/ml)', 'Rises as the cord is squeezed.'], R_umv: ['Umbilical vein resistance (mmHg·s/ml)', 'Collapses first under compression.'],
  p_a: ['Fetal arterial pressure (mmHg)', 'Transmural. Rest about 45 mmHg.'], p_v: ['Fetal venous pressure (mmHg)', ''],
  p_vil: ['Placental villi pressure (mmHg)', 'Fetal blood inside the placental villi.'],
  q: ['Fetal blood flows (ml/s)', 'Heart output, body, brain, and through the cord to the villi and back.'], V: ['Fetal blood volumes (ml)', ''],
  R_cera: ['Fetal brain artery resistance (mmHg·s/ml)', 'Falls when oxygen is low (brain sparing).'],
  T: ['Fetal heart period (s)', 'FHR = 60 / T. 0.444 s at 135 bpm.'],
  z_Tv: ['Vagal effect on heart period', 'Lengthens T. Main driver of decelerations.'], z_Ts: ['Sympathetic effect on heart period', 'Shortens T.'],
  z_E: ['Contractility effect', ''], z_R: ['Peripheral resistance effect', ''], z_V: ['Venous volume effect', ''],
  pO2_a: ['Fetal arterial pO2 (mmHg)', 'Rest about 18.5 mmHg (Jongen mother).'], c_a: ['Fetal arterial O2 content (ml O2 / ml blood)', ''],
  pO2_i: ['pO2 of each blood pool (mmHg)', 'Fetal and maternal pools.'], pO2_tis: ['Fetal body-tissue pO2 (mmHg)', ''], pO2_cer: ['Fetal brain pO2 (mmHg)', ''],
  'pO2_a,m': ['Maternal arterial pO2 (mmHg)', ''], 'pO2_tis,m': ['Maternal tissue pO2 (mmHg)', ''],
  SpO2_f: ['Fetal arterial saturation (%)', 'What a pulse oximeter would read.'], SpO2_m: ['Maternal arterial saturation (%)', ''],
  StO2_f_tis: ['Fetal body-tissue saturation (%)', '25 % arterial, 75 % body-tissue pool.'], StO2_f_cer: ['Fetal brain saturation (%)', '25 % arterial, 75 % brain pool.'],
  StO2_m: ['Maternal tissue saturation (%)', ''], M_i: ['O2 amount per compartment (ml O2)', ''],
  r_c: ['Chemoreceptor output (−1 to 1)', 'Rises when fetal pO2 falls below its set point.'], r_b: ['Baroreceptor output (−1 to 1)', ''],
  e_v: ['Vagal (parasympathetic) signal', ''], e_s: ['Sympathetic signal', ''],
  delay: ['Delay buffer', 'Recent history of the signal, so the block can use its value from D seconds ago.'],
  all: ['Every block state', 'The run starts all of them from the steady state at rest.'],
  a: ['Contraction strength', 'Pulse height as a multiple of 70 mmHg.'], w_um: ['Cord-compression weight', 'Share of the contraction pressure that squeezes the cord.'],
  g: ['Reflex gain', 'Multiplies all five fetal effector gains. 1 = published values.'],
  mode: ['Model version', 'The switches of run() for this version.'],
  R_uta: ['Uterine artery resistance', 'One resistance for the uterine and spiral arteries; not squeezed by contractions.'],
  'R_utv law': ['Uterine vein collapse law', 'R = R0 / f² with f an S-curve of transmural pressure: the vein collapses as the uterus squeezes it.'],
  'IVS law': ['Intervillous space pressure–volume curve', ''],
};
const look = k => G[k] || G[k.split(' ')[0]] || G[k.replace(/[,(]$/, '')] || G[k.split(' ')[0].replace(/[,(]$/, '')];

/* ---------- rendering ---------- */
const mk = (tag, cls, text) => { const e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; };
const base = v => v.split(' ')[0];
function chip(kind, label, v, title) {
  const c = mk('span', 'bkchip ' + kind); c.dataset.v = base(v); c.dataset.g = v; c.dataset.from = title || '';
  if (typeof label === 'string') c.textContent = label; else c.append(...label);
  return c;
}
function row(label, chips) { const r = mk('div', 'bkrow'); r.append(mk('b', null, label)); const c = mk('div', 'bkchips'); c.append(...chips); r.append(c); return r; }
function eqBlock(text, outs) {
  const pre = mk('pre', 'bkeqs');
  const esc = s => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const re = new RegExp(`(?<![\\w,])(${[...outs].sort((a, b) => b.length - a.length).map(esc).join('|')})(?![\\w]|,[a-z])`, 'g');
  let last = 0;
  for (const mm of text.matchAll(re)) { pre.append(text.slice(last, mm.index), mk('span', 'ov', mm[0])); last = mm.index + mm[0].length; }
  pre.append(text.slice(last));
  return pre;
}
function card(b) {
  const el = mk('article', 'bkcard'); el.tabIndex = 0; el.dataset.k = b.k;
  const h = mk('h3'); h.append(mk('span', 'n', b.n), document.createTextNode(b.name)); el.append(h);
  const sig = mk('div', 'sig'); sig.append(mk('span', 'fnm', b.fn), '(');
  b.in.forEach(([v, , fb], i) => { if (i) sig.append(', '); sig.append(chip('in' + (fb ? ' fb' : ''), v, v)); });
  sig.append(') → ');
  b.out.forEach((v, i) => { if (i) sig.append(', '); sig.append(chip('out', v, v)); });
  el.append(sig);
  el.append(row('inputs', b.in.map(([v, from, fb]) => chip('in' + (fb ? ' fb' : ''), `${v} ← ${from}`, v,
    fb ? `from block ${from}, read from the state at the start of the step` : `from ${/^M?\d/.test(from) ? 'block ' : ''}${from}`))));
  el.append(row('outputs', b.out.map(v => chip('out', v, v))));
  if (b.par.length) el.append(row('constants', b.par.map(([s, val, src]) => {
    const lab = [document.createTextNode(val ? `${s} = ${val} ` : `${s} `), mk('span', 'src', `(${src})`)];
    return chip('par', lab, s, `source: ${src}`);
  })));
  if (b.st.length) el.append(row('state', b.st.map(v => chip('st', v, v))));
  el.append(eqBlock(b.eq, b.out.map(base)));
  el.append(row('code', [mk('span', 'bknote', b.file)]));
  return el;
}
const tip = mk('div', 'bk-tip'); tip.hidden = true; document.body.append(tip);
function tipOn(target, keys, extra) {
  const move = e => { const r = tip.getBoundingClientRect();
    let x = e.clientX + 14, y = e.clientY + 16;
    if (x + r.width > innerWidth - 8) x = Math.max(8, innerWidth - r.width - 8);
    if (y + r.height > innerHeight - 8) y = Math.max(8, e.clientY - r.height - 12);
    tip.style.left = x + 'px'; tip.style.top = y + 'px'; };
  target.addEventListener('pointerenter', e => {
    const ents = keys().map(k => [k, look(k)]).filter(([, g]) => g);
    if (!ents.length) return;
    tip.replaceChildren(...ents.map(([k, [nm, mn]]) => { const d = mk('div'); d.append(mk('b', null, k), ' ', mk('span', 'nm', nm));
      if (mn) { d.append(mk('br')); d.append(mk('span', 'mn', mn)); } return d; }));
    const ex = extra && extra(); if (ex) tip.append(mk('div', 'mn', ex));
    tip.hidden = false; move(e);
  });
  target.addEventListener('pointermove', e => { if (!tip.hidden) move(e); });
  target.addEventListener('pointerleave', () => { tip.hidden = true; });
}

const R = 12, TP = 10, tmax = 20, NS = 'http://www.w3.org/2000/svg';
const sv = (tag, attrs, parent) => { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); parent && parent.appendChild(e); return e; };
function niceTicks(lo, hi, n = 5) {
  const raw = (hi - lo || 1) / n, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map(x => x * mag).find(s => s >= raw);
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(6));
  return out;
}
let D = null, CH = null, pinned = 'out', model = -1;
const MAIN = () => ({title: $b('#bk-ptitle'), note: $b('#bk-pnote'), ro: $b('#bk-readout'), leg: $b('#bk-plegend'), chart: $b('#bk-chart')});
function draw(k, box = MAIN()) {
  const b = CH[k] || CH.out, W = box.W || 1000, H = box.H || 220, L = box.L || 56, sw = box.sw || 1.8, BT = box.BT || 30;
  box.title.textContent = b.title; box.note.textContent = b.note || '';
  const ro = box.ro; ro.textContent = '';
  const leg = box.leg; leg.innerHTML = '';
  b.s.forEach(([, , lab], i) => { const g = look(lab), long = g && g[0].replace(/\s*\([^)]*\)$/, '');
    const s = mk('span', ['', 'd', 'dt'][i], long ? `${long} · ${lab}` : lab); s.style.setProperty('--c', `var(--bk-s${i + 1})`); leg.append(s); });
  let lo = Infinity, hi = -Infinity;
  b.s.forEach(([key]) => D[key].forEach(v => { if (v < lo) lo = v; if (v > hi) hi = v; }));
  const pad = (hi - lo) * 0.08 || 1; lo -= pad; hi += pad;
  const X = mm => L + (mm / tmax) * (W - L - R), Y = v => TP + (1 - (v - lo) / (hi - lo)) * (H - TP - BT);
  const chart = box.chart; chart.innerHTML = '';
  const svg = sv('svg', {viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': `${b.title} over 20 minutes`}, chart);
  niceTicks(lo, hi).forEach(v => { sv('line', {class: 'gl', x1: L, x2: W - R, y1: Y(v), y2: Y(v)}, svg);
    sv('text', {class: 'ax', x: L - 6, y: Y(v) + 4, 'text-anchor': 'end'}, svg).textContent = v; });
  for (let mm = 0; mm <= tmax; mm += 2.5) sv('text', {class: 'ax', x: X(mm), y: H - BT / 3, 'text-anchor': 'middle'}, svg).textContent = mm;
  sv('text', {class: 'ax', x: W - R, y: H - BT / 3, 'text-anchor': 'end'}, svg).textContent = 'min';
  if (b.unit) sv('text', {class: 'ax', x: 4, y: TP + 2, 'dominant-baseline': 'hanging'}, svg).textContent = b.unit;
  D.peaks.forEach(mm => sv('line', {class: 'pk', x1: X(mm), x2: X(mm), y1: TP, y2: H - BT}, svg));
  b.s.forEach(([key, ts], i) => {
    const pts = D[key].map((v, j) => `${X(ts[j] / 60).toFixed(1)},${Y(v).toFixed(1)}`).join(' ');
    sv('polyline', {points: pts, fill: 'none', stroke: `var(--bk-s${i + 1})`, 'stroke-width': sw,
      'stroke-dasharray': ['none', '6 4', '2 3'][i], 'stroke-linejoin': 'round'}, svg);
  });
  const xh = sv('line', {class: 'xh', y1: TP, y2: H - BT, x1: -10, x2: -10}, svg);
  svg.addEventListener('pointermove', e => {
    const r = svg.getBoundingClientRect(), x = (e.clientX - r.left) / r.width * W;
    const mm = Math.max(0, Math.min(tmax, (x - L) / (W - L - R) * tmax)), s = mm * 60;
    xh.setAttribute('x1', X(mm)); xh.setAttribute('x2', X(mm));
    const vals = b.s.map(([key, ts]) => { let j = 0, best = Infinity;
      ts.forEach((tt, i) => { const d = Math.abs(tt - s); if (d < best) { best = d; j = i; } }); return D[key][j]; });
    ro.textContent = `${Math.floor(s / 60)}:${String(Math.round(s % 60)).padStart(2, '0')}  ` +
      vals.map((v, i) => `${b.s[i][2].split(' ')[0]} ${v}`).join('   ') + (b.unit ? `  ${b.unit}` : '');
  });
  svg.addEventListener('pointerleave', () => { xh.setAttribute('x1', -10); xh.setAttribute('x2', -10); ro.textContent = ''; });
}
function render(m) {
  if (m === model) return; model = m;
  D = Object.assign({}, DATA, DATA['cm' + m]); CH = charts(m, D);
  $b('#bk-flow').innerHTML = diagram(m);
  $b('#bk-pipe').replaceChildren(...pipe(m).map(card));
  $b('#bk-mother-h').textContent = `Mother: block${MOTHER[m].length > 1 ? 's M1 to M' + MOTHER[m].length : ' M1'}`;
  $b('#bk-mother').replaceChildren(...MOTHER[m].map(card));
  $b('#bk-fetus').replaceChildren(...[m === 2 ? PLACENTA_CK : PLACENTA_J('p_a,m', 'p_v,m', m === 0 ? 'M1' : 'M2')].concat(FETAL).map(card));
  root.querySelectorAll('.bkcard').forEach(c => c.classList.toggle('bkkey', KEYS[m].includes(c.dataset.k)));
  $b('#bk-version').textContent = ['Jongen 2016', 'fetosim: + maternal layers', 'fetosim: Ursino + Clark'][m];
  const chips = [...root.querySelectorAll('.bkchip[data-v]')];
  chips.forEach(c => {
    c.addEventListener('pointerenter', () => chips.forEach(o => o.classList.toggle('hl', o.dataset.v === c.dataset.v)));
    c.addEventListener('pointerleave', () => chips.forEach(o => o.classList.remove('hl')));
    tipOn(c, () => [c.dataset.g], () => c.dataset.from);
  });
  root.querySelectorAll('.bkflow .wl, .bkflow .fbl').forEach(t => tipOn(t, () => t.textContent.split(/,?\s+/).filter(Boolean)));
  cards = [...root.querySelectorAll('.bkcard, .bkflow .nd')];
  if (!CH[pinned]) pinned = 'out';
  const hov = k => { if (typeof hoverPart === 'function') hoverPart(k ? GROUP_OF_BLOCK(k) : null); };
  cards.forEach(c => {
    const k = c.dataset.k;
    c.addEventListener('pointerenter', () => { show(k); hov(k); });
    c.addEventListener('pointerleave', () => { show(pinned); hov(null); });
    c.addEventListener('click', () => { pinned = k; show(k); pickPart(k);
      if (c.closest('.bkflow')) { const d = root.querySelector(`.bkcard[data-k="${k}"]`); if (d) d.scrollIntoView({behavior: 'smooth', block: 'start'}); } });
    c.addEventListener('focus', () => { pinned = k; show(k); });
  });
  show(pinned); markPart(typeof current !== 'undefined' ? current : null); mini();
}
let cards = [];
function show(k) { draw(k); cards.forEach(c => c.classList.toggle('on', c.dataset.k === k)); }

/* ---------- link with the cross-section above (numbered parts: GROUP keys of the explorer) ---------- */
const PART_BLOCKS = m => ({
  mother: MOTHER[m].map(b => b.k), uterus: ['toco', 'put'], spiral: ['placenta'], placenta: ['placenta', 'oxygen', 'sat'],
  cord: ['cord'], heart: ['circ', 'brain', 'reflex', 'effectors', 'heart'], baro: ['baro'], chemo: ['chemo']});
const PART_CHART = m => ({mother: ['mcirc', 'mcirc', 'uhv'][m], uterus: 'put', spiral: 'placenta', placenta: 'oxygen', cord: 'cord',
  heart: 'heart', baro: 'baro', chemo: 'chemo'});
function GROUP_OF_BLOCK(k) {
  if (k === 'placenta') return 'spiral';
  if (k === 'oxygen' || k === 'sat') return 'placenta';
  const pb = PART_BLOCKS(model); return Object.keys(pb).find(g => pb[g].includes(k)) || null;
}
let picking = false;
function pickPart(k) {                       // a block click selects its part above (and zooms the cross-section there)
  const g = GROUP_OF_BLOCK(k); if (!g || typeof select !== 'function') return;
  picking = true; select(g, false); picking = false;
}
function markPart(g) {                       // outline the blocks of the selected part
  const ks = g ? PART_BLOCKS(model)[g] || [] : [];
  cards.forEach(c => c.classList.toggle('bkpart', ks.includes(c.dataset.k)));
}
/* the explorer's part panel gets the signal of the selected part (FHR and MHR when nothing is selected) */
function mini() {
  const panel = document.getElementById('panel'); if (!panel || !CH) return;
  let box = panel.querySelector('.bk-mini');
  if (!box) {
    box = mk('div', 'bk-mini');
    box.innerHTML = '<div class="bk-phead"><h4></h4><span class="bk-readout"></span></div><div class="bk-plegend"></div><div class="bk-mchart"></div><p class="bknote"></p><button type="button" class="bk-go">Blocks and signals ↓</button>';
    const h = panel.querySelector('.ptitle'); h ? h.after(box) : panel.prepend(box);
    box.querySelector('.bk-go').addEventListener('click', () => {
      const g = typeof current !== 'undefined' ? current : null, ks = g ? PART_BLOCKS(model)[g] : null;
      const d = ks && root.querySelector(`.bkcard[data-k="${ks[0]}"]`);
      (d || root).scrollIntoView({behavior: 'smooth', block: 'start'});
    });
  }
  const g = typeof current !== 'undefined' ? current : null;
  draw(g ? PART_CHART(model)[g] : 'out', {title: box.querySelector('h4'), note: box.querySelector('.bknote'), ro: box.querySelector('.bk-readout'),
    leg: box.querySelector('.bk-plegend'), chart: box.querySelector('.bk-mchart'), W: 1000, H: 420, L: 120, sw: 4, BT: 56});
}
if (typeof hooks !== 'undefined') hooks.select.push(g => {
  markPart(g);
  if (!picking && g) { const k = PART_CHART(model)[g]; if (CH[k]) { pinned = k; show(k); } }
});
const pnl = document.getElementById('panel');
if (pnl) new MutationObserver(() => { if (!pnl.querySelector('.bk-mini')) mini(); }).observe(pnl, {childList: true});
document.addEventListener('modelchange', e => render(+e.detail.cm));
render(+(window.FS_MODEL || 0));
})();
