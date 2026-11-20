"""Run the maternal-fetal model from Python: uterine pressure p_ut(t) in -> FHR, pressures, flows and pO2 out.

The ODEs are solved by core.c (forward Euler, compiled with the system C compiler on first use into build/ next to
this file; build/ is not tracked by git). The parameter vector comes from a JSON parameter set in parameters/
(default parameters/jongen2016.json), whose entry order is the index order of core.c.

The equations are those of Jongen et al. 2016 (J Biomech 49(12), Part A: 2466-2473, Part B: 2474-2480) and the FHR
variability noise of Jongen's 2016 thesis Ch. 4. Tags: A(8) = Part A Eq. 8, B(2) = Part B Eq. 2, A-T1a = Part A
Table 1a. Two additions that are not in the papers: S_GAIN (g, a common multiplier of the five effector gains k of
A(23)) and the inputs pum (extra umbilical pressure) and ruta (R_uta multiplier).

    from fetosim.model import simulate as S          # with the folder that holds fetosim/ on sys.path
    out, state, ok = S.run(S.contraction(300, 0.02), w_um=1.0, du=0.02)   # Part B reference contraction
"""
import ctypes
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PARAMETER_SET = HERE / 'parameters/jongen2016.json'
MATERNAL_HR_SET = HERE / 'parameters/maternal_hr_randall2019.json'   # appended after Jongen 2016 (core.c MR_*)
MATERNAL_THORAX_SET = HERE / 'parameters/maternal_thorax_lu2001.json'  # appended after that (core.c MT_*)
MATERNAL_VASO_SET = HERE / 'parameters/maternal_vasomotor.json'        # appended after that (core.c MV_*, MA_*)
MATERNAL_CIRC_SETS = (HERE / 'parameters/maternal_circulation_coupling.json',     # core.c MU_MODE, MU_NSUB
                      HERE / 'parameters/maternal_circulation_ursino1998.json',   # then Ursino 1998 (ursino1998.c U_*)
                      HERE / 'parameters/placental_o2_hill1972.json',             # then Hill 1972 (hill1972.c H_*)
                      HERE / 'parameters/uterine_network_clark2018.json')         # then Clark 2018 network (core.c CK_*)


def load_parameter_set(path=PARAMETER_SET):
    """List of (name, value, source) in core.c order, plus the full entries. JSON floats round-trip exactly (repr),
    so the values are bit-identical to the ones written."""
    entries = json.loads(Path(path).read_text())['parameters']
    return [(e['name'], e['value'], e['source_detail']) for e in entries], entries


PARAMS, ENTRIES, NAMES, IDX = [], [], [], {}


def reload_parameters():
    """(Re)read the parameter files into PARAMS, ENTRIES, NAMES, IDX and drop the cached steady state."""
    PARAMS.clear(); ENTRIES.clear(); NAMES.clear(); IDX.clear()
    for path in (PARAMETER_SET, MATERNAL_HR_SET, MATERNAL_THORAX_SET, MATERNAL_VASO_SET) + MATERNAL_CIRC_SETS:
        p, e = load_parameter_set(path); PARAMS.extend(p); ENTRIES.extend(e)
    NAMES.extend(p[0] for p in PARAMS); IDX.update({n: i for i, n in enumerate(NAMES)})
    _SS.clear()


_SS = {}
reload_parameters()
STATE = ['VM_ART', 'VM_VEN', 'VM_TIS', 'VM_IVS', 'MM_ART', 'MM_VEN', 'MM_TIS', 'MM_IVS', 'MM_H',
         'VF_ART', 'VF_VEN', 'VF_TIS', 'VF_CER', 'VF_VIL', 'MF_ART', 'MF_VEN', 'MF_TIS', 'MF_CER', 'MF_VIL', 'MF_H',
         'HM_PREV', 'HF_PREV', 'RCERA', 'RB', 'RC', 'XTV', 'XTS', 'XE', 'XR', 'XV'] + [f'PG{i}' for i in range(11)] + \
        ['R_EBC', 'R_EBA', 'R_TPB', 'R_TPR', 'R_TS', 'R_H', 'VM_TH', 'MM_TH', 'MV_XR', 'MV_XV'] + \
        ['MU_' + n for n in ('PPA', 'FPA', 'PPP', 'PPV', 'PSA', 'FSA', 'PSP', 'PEV', 'PLA', 'VLV', 'PRA', 'VRV', 'XI',
                             'PT', 'DELV', 'DERV', 'DRSP', 'DREP', 'DVSV', 'DVEV', 'DTS', 'DTV', 'VT')]  # ursino1998.c
SI = {n: i for i, n in enumerate(STATE)}
OUT = ['fhr', 'ptm_a', 'po2_a', 'q_uta', 'q_utv', 'v_ivs', 'po2_tis', 'po2_cer', 'q_tisa', 'q_cera', 'q_uma', 'q_umv',
       'v_villi', 'r_b', 'r_c', 'p_ut', 'e_s', 'e_v', 'z_tv', 'z_ts', 'q_co',
       'mhr', 'pa_m', 'q_co_m', 'sbp_m', 'pth_m', 'pv_m', 'q_vr_m', 'pth_tm_m', 'rtisa_m', 'v_sh_m', 'ts_m',
       'po2_m_a', 'po2_m_tis', 'po2_v',
       'q_ut_art', 'q_anast', 'p_node']
# q_ut_art, q_anast (ml/s), p_node (mmHg): uterine artery flow, anastomotic bypass flow and network node pressure of the
# Clark 2018 uterine network (CK_MODE 1); q_uta is then the spiral-artery flow into the intervillous space. With
# CK_MODE 0: q_ut_art = q_uta, q_anast = 0, p_node = pa_m.
# Oxygen saturations, %, added to every output dict by run_core (not computed in core.c):
#   spo2_m, spo2_f   arterial saturation of mother and fetus (what a pulse oximeter reads)
#   sto2_m, sto2_f_tis, sto2_f_cer   tissue saturation as NIRS reads it: STO2_ART x arterial + (1 - STO2_ART) x the
#                    tissue (cerebral) pool, which holds blood that has given up its O2 (venous side).
# Saturation is the hemoglobin term of A(13), g / (g + c1) with g = p^3 + c2 p, using the same curves as the model.
STO2_ART = 0.25     # arterial share of the NIRS volume (our choice; the common 25:75 arterial:venous assumption)
MR_RANDALL = ['MR_A', 'MR_B', 'MR_KB', 'MR_KPB', 'MR_KPR', 'MR_KS', 'MR_TAUB', 'MR_TAUPB', 'MR_TAUPR', 'MR_TAUS',
              'MR_TAUH', 'MR_QW', 'MR_QPB', 'MR_QPR', 'MR_QS', 'MR_SW', 'MR_SPB', 'MR_SPR', 'MR_SS', 'MR_HI', 'MR_HPB',
              'MR_HPR', 'MR_HS', 'MR_DS']                         # same order as randall2019.NAMES

_lib = None
_R = None


def _randall():
    """model/randall2019.py, loaded by path so this file works both as fetosim.model.simulate and on its own."""
    global _R
    if _R is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location('fetosim_randall2019', HERE / 'randall2019.py')
        _R = importlib.util.module_from_spec(spec); spec.loader.exec_module(_R)
    return _R


def lib():
    global _lib
    if _lib is None:
        src = HERE / 'core.c'
        tag = hashlib.sha256(src.read_bytes() + (HERE / 'ursino1998.c').read_bytes()
                             + (HERE / 'hill1972.c').read_bytes()).hexdigest()[:12]
        so = HERE / f'build/core_{tag}.so'
        if not so.exists():
            so.parent.mkdir(exist_ok=True)
            subprocess.run(['cc', '-O3', '-shared', '-fPIC', '-o', str(so), str(src), '-lm'], check=True)
        L = ctypes.CDLL(str(so)); d = ctypes.POINTER(ctypes.c_double)
        L.simulate.argtypes = [d, d, ctypes.c_int, ctypes.c_double, d, d, d, d, d, ctypes.c_int, ctypes.c_double, ctypes.c_int, d]
        L.simulate.restype = ctypes.c_int
        L.hill_exchange_rate.argtypes = [d] + [ctypes.c_double] * 4 + [d]; L.hill_exchange_rate.restype = ctypes.c_double
        assert L.n_par() == len(NAMES) and L.n_state() == len(STATE) and L.n_out() == len(OUT)
        _lib = L
    return _lib


def params(**kw):
    P = np.array([p[1] for p in PARAMS], float)
    for k, v in kw.items():
        P[IDX[k]] = v
    return P


def _ptr(a):
    return None if a is None else a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def cO2(p, fetal=True):
    ahb = 1.34 * (0.17 if fetal else 0.12); c1 = 1.04e4 if fetal else 2.34e4; g = p ** 3 + 150 * p
    return ahb * g / (g + c1) + 3.1e-5 * p


def initial_state(P):
    """Derived steady-state guess from the papers' stated targets (B §2.2.1-2.2.2); exact state from steady_state()."""
    g = lambda n: P[IDX[n]]
    S = np.zeros(len(STATE))
    pa_m, pv_m, pt_m = 80.0, 5.0, 9.25
    S[SI['VM_ART']] = g('M_VART0') + g('M_CART') * pa_m; S[SI['VM_VEN']] = g('M_VVEN0') + g('M_CVEN') * pv_m
    S[SI['VM_TIS']] = g('M_VTIS0') + g('M_CTIS') * pt_m; S[SI['VM_IVS']] = 175.0
    pa, pv, pt, pc, pl = 45.0, 3.0, 5.13, 5.09, 19.8
    S[SI['VF_ART']] = g('F_VART0') + g('F_CART') * pa; S[SI['VF_VEN']] = g('F_VVEN0') + g('F_CVEN') * pv
    S[SI['VF_TIS']] = g('F_VTIS0') + g('F_CTIS') * pt; S[SI['VF_CER']] = g('F_VCER0') + g('F_CCER') * pc
    S[SI['VF_VIL']] = g('F_VVIL0') + g('F_CVIL') * pl
    hm = 0.5 * (g('M_V0MIN') + pv_m / g('M_EMIN') + g('M_V0MAX') + pa_m / g('M_EMAX'))
    hf = 0.5 * (g('F_V0MIN') + pv / g('F_EMIN') + g('F_V0MAX') + pa / g('F_EMAX'))
    S[SI['HM_PREV']] = hm; S[SI['HF_PREV']] = hf
    po2 = {'M_ART': 98, 'M_VEN': 98, 'M_TIS': 40, 'M_IVS': 40, 'M_H': 98, 'F_ART': 18.5, 'F_VEN': 16, 'F_TIS': 15.4,
           'F_CER': 15.4, 'F_VIL': 22, 'F_H': 16}
    vol = {'M_ART': S[SI['VM_ART']], 'M_VEN': S[SI['VM_VEN']], 'M_TIS': S[SI['VM_TIS']], 'M_IVS': 175.0, 'M_H': hm,
           'F_ART': S[SI['VF_ART']], 'F_VEN': S[SI['VF_VEN']], 'F_TIS': S[SI['VF_TIS']], 'F_CER': S[SI['VF_CER']],
           'F_VIL': S[SI['VF_VIL']], 'F_H': hf}
    for i, (k, p) in enumerate(po2.items()):
        S[SI[('MM_' if k[0] == 'M' else 'MF_') + {'ART': 'ART', 'VEN': 'VEN', 'TIS': 'TIS', 'IVS': 'IVS', 'H': 'H', 'CER': 'CER', 'VIL': 'VIL'}[k[2:]]]] = cO2(p, k[0] == 'F') * vol[k]
        S[SI['PG0'] + i] = p
    S[SI['RCERA']] = g('F_RCERA')
    set_reflex_rest(P, S)
    return S


def set_reflex_rest(P, S):
    """Put the maternal HR reflex states (Randall 2019) at rest: SBP = MR_SW (the resting SBP), Pth = MR_PTH_REST and
    MHR = 60 / M_T. The nominal Randall parameters are built so that this is a steady state."""
    R = _randall()
    p = np.array([P[IDX[n]] for n in MR_RANDALL])
    y0 = R.initial_conditions(p, P[IDX['MR_SW']], P[IDX['MR_PTH_REST']], 60.0 / P[IDX['M_T']])
    for n, v in zip(['R_EBC', 'R_EBA', 'R_TPB', 'R_TPR', 'R_TS', 'R_H'], y0):
        S[SI[n]] = v
    return S


def set_thorax_rest(P, S):
    """Maternal thorax (MT_MODE 1) at Jongen's rest: move MT_VTH ml (and its O2) from the venous pool into the
    intrathoracic reservoir, set MT_DPF = MT_RVR x resting CO and MT_V0TH so that the heart's filling pressure, and
    therefore the whole resting state, equals Jongen's. Changes P and S in place."""
    g = lambda n: P[IDX[n]]
    pa = maternal_pa(S, P); pv = (S[SI['VM_VEN']] - g('M_VVEN0')) / g('M_CVEN')
    q = ((g('M_V0MIN') + pv / g('M_EMIN')) - (g('M_V0MAX') + pa / g('M_EMAX'))) / g('M_T')
    P[IDX['MT_DPF']] = g('MT_RVR') * q
    P[IDX['MT_V0TH']] = g('MT_VTH') - g('MT_CTH') * (pv - g('MT_DPF'))
    c = S[SI['MM_VEN']] / S[SI['VM_VEN']]
    S[SI['VM_VEN']] -= g('MT_VTH'); S[SI['MM_VEN']] -= c * g('MT_VTH')
    S[SI['VM_TH']] = g('MT_VTH'); S[SI['MM_TH']] = c * g('MT_VTH')
    return P, S


def run_core(P, S, put, dt=0.02, du=None, pum=None, ruta=None, noise=None, pth=None, out_every=None, copy_state=True):
    """Low-level run with a full parameter vector P: len(put)*du/dt steps from state S (copied unless copy_state=False).
    pth is the maternal thoracic pressure for the HR reflex and, with MT_MODE 1, the thorax (None = MR_PTH_REST). Returns (outputs dict, final state, ok)."""
    du = dt if du is None else du
    put = np.ascontiguousarray(put, float); n_in = len(put)
    n_steps = int(round((n_in - 1) * du / dt)) + 1
    out_every = out_every or max(1, int(round(0.5 / dt)))
    nrec = (n_steps + out_every - 1) // out_every
    out = np.zeros(nrec * len(OUT)); S = np.array(S, float) if copy_state else S
    arr = [None if a is None else np.ascontiguousarray(a, float) for a in (pum, ruta, noise, pth)]
    for a in arr:
        assert a is None or len(a) == n_in
    P = np.ascontiguousarray(P, float)
    bad = lib().simulate(_ptr(P), _ptr(S), n_steps, dt, _ptr(put), _ptr(arr[0]), _ptr(arr[1]), _ptr(arr[2]),
                         _ptr(arr[3]), n_in, du, out_every, _ptr(out))
    o = out.reshape(nrec, len(OUT))
    res = {k: o[:, i] for i, k in enumerate(OUT)}; res['t'] = np.arange(nrec) * out_every * dt
    add_saturations(res, P)
    return res, S, bad == 0


def saturation(p, P, fetal=True):
    """Hemoglobin O2 saturation, %, for pO2 p (mmHg) on the model's dissociation curve (hemoglobin term of A(13))."""
    c1, c2 = (P[IDX['O_C1F']], P[IDX['O_C2F']]) if fetal else (P[IDX['O_C1M']], P[IDX['O_C2M']])
    p = np.maximum(p, 0.0); g = p ** 3 + c2 * p
    return 100 * g / (g + c1)


def add_saturations(res, P):
    """Add spo2_m, spo2_f, sto2_m, sto2_f_tis, sto2_f_cer (%) to an output dict (see OUT)."""
    res['spo2_m'] = saturation(res['po2_m_a'], P, fetal=False); res['spo2_f'] = saturation(res['po2_a'], P)
    mix = lambda a, v: STO2_ART * a + (1 - STO2_ART) * v
    res['sto2_m'] = mix(res['spo2_m'], saturation(res['po2_m_tis'], P, fetal=False))
    res['sto2_f_tis'] = mix(res['spo2_f'], saturation(res['po2_tis'], P))
    res['sto2_f_cer'] = mix(res['spo2_f'], saturation(res['po2_cer'], P))


def steady_state(P, dt=0.02, t_max=3000.0, tol=1e-4):
    """Run at rest until |d p_tm,a/dt| and |d pO2,a/dt| < tol mmHg/s (A s2.4) over the last 10 s."""
    S = initial_state(P); rest = np.full(int(100 / dt) + 1, P[IDX['S_PREST']]); t = 0.0
    while t < t_max:
        o, S, ok = run_core(P, S, rest, dt=dt, out_every=int(round(1 / dt)))
        assert ok, 'steady-state run became non-finite'
        t += 100
        rate = max(abs(o['ptm_a'][-1] - o['ptm_a'][-11]) / 10, abs(o['po2_a'][-1] - o['po2_a'][-11]) / 10)
        if rate < tol:
            return S
    raise RuntimeError(f'no steady state within {t_max} s (rate {rate:.2e})')


def contraction(t_total, dt, t_con=30.0, T_con=60.0, p_con=70.0, p_rest=15.0):
    """Jongen 2016 A(8): p_ut = p_rest + p_con sin^2(pi (t - t_con) / T_con) on (t_con, t_con + T_con), else p_rest.
    Sampled every dt; the reference contraction is p_con 70 mmHg, T_con 60 s (A §2.4)."""
    t = np.arange(int(round(t_total / dt)) + 1) * dt
    p = np.full(len(t), p_rest); m = (t > t_con) & (t < t_con + T_con)
    p[m] += p_con * np.sin(np.pi * (t[m] - t_con) / T_con) ** 2
    return p


# ---- convenience layer used by the input and fit code -----------------------------------------------------------
DTS = 0.02                      # model time step, s (Jongen 2016 A §2.4: forward Euler, dt 0.02 s)
DU = 0.5                        # input grid for p_ut, s (linearly interpolated inside the core)
OE = 5                          # core output every 5 steps = 0.1 s
TONE = 15.0                     # resting uterine pressure p_rest, mmHg (A §2.4)
P_CON = 70.0                    # reference contraction amplitude above rest, mmHg (A §2.4; peak 85 mmHg)


def steady():
    """Steady state of the default parameter set at rest (A §2.4: run until |dp_tm,a/dt| and |dpO2,a/dt| < 1e-4
    mmHg/s). Cached. w_um and g do not change the resting state (p_um = w_um (p_ut - p_rest) = 0 and all effector
    inputs are 0 at rest), so every fit starts from this one state."""
    if 'S' not in _SS:
        _SS['S'] = steady_state(params())
    return _SS['S']


def hill_exchange(pm_in, pf_in, qm, qf, **kw):
    """Hill 1972 exchanger alone: inflowing maternal and fetal pO2 (mmHg), capillary flows (ml/s); kw overrides H_*
    parameters. Returns (O2 transfer ml/s, maternal end-capillary pO2, fetal end-capillary pO2, mean Dp)."""
    P = params(**kw); H = np.ascontiguousarray(P[IDX['H_MODE']:IDX['CK_MODE']]); end = np.zeros(3)
    r = lib().hill_exchange_rate(_ptr(H), pm_in, pf_in, qm, qf, _ptr(end))
    return r, end[0], end[1], end[2]


def steady_hill():
    """Resting state with the Hill 1972 placental exchanger (H_MODE 1) and Jongen's mother. Cached."""
    if 'SH' not in _SS:
        _SS['SH'] = steady_state(params(H_MODE=1))
    return _SS['SH']


def steady_ursino(t_settle=600.0, placenta=0):
    """Resting state with the Ursino 1998 maternal circulation (MU_MODE 1): Jongen's resting state with the Ursino
    states at their own basal state, then t_settle s at rest so that the uterine branch, intervillous space and both
    O2 systems settle together. placenta=1: with the Hill 1972 exchanger. Cached."""
    key = 'SU%d' % placenta
    if key not in _SS:
        import importlib.util
        spec = importlib.util.spec_from_file_location('fetosim_ursino1998', HERE / 'ursino1998.py')
        U = importlib.util.module_from_spec(spec); spec.loader.exec_module(U)
        P = params(MU_MODE=1, H_MODE=placenta); S = (steady_hill() if placenta else steady()).copy()
        i0 = SI['MU_PPA']; S[i0:i0 + len(U.USTATE)] = U.steady_state(P[IDX['U_CSA']:IDX['U_CSA'] + len(U.NAMES)])
        rest = np.full(int(t_settle / DU) + 1, P[IDX['S_PREST']])
        _, S, ok = run_core(P, S, rest, dt=DTS, du=DU, out_every=50)
        assert ok, 'Ursino steady-state run became non-finite'
        _SS[key] = S
    return _SS[key]


def clark_parameters(S0, P, k_a=None):
    """Parameter overrides that switch on the uterine network after Clark 2018 (CK_MODE 1) without changing the
    resting state S0 of parameter vector P (which has CK_MODE 0). k_a = CK_KA (None: the value in P).
    CK_RUP and CK_RA are Clark's. Derived from the resting state:
      CK_RSA            so that the resting intervillous inflow and pressure stay as they are;
      M_RTISA, M_RTISV  both scaled by q_tis / (q_tis - q_anast): the bypass flow is taken out of the maternal tissue
                        flow, so maternal cardiac output and arterial and venous pressure stay as they are.
    Returns (overrides, info); info holds the resting flows (ml/s) and pressures (mmHg)."""
    g = lambda n: P[IDX[n]]
    k_a = g('CK_KA') if k_a is None else k_a
    pa = maternal_pa(S0, P); pv = (S0[SI['VM_VEN']] - g('M_VVEN0')) / g('M_CVEN')
    pt = (S0[SI['VM_TIS']] - g('M_VTIS0')) / g('M_CTIS')
    p_ivs = g('IVS_P0') + g('IVS_P1') * np.tan(np.pi * (S0[SI['VM_IVS']] / g('IVS_VMAX') - 0.5)) + g('S_PREST')   # A(3)
    q_ivs = (pa - p_ivs) / g('M_RUTA'); q_tis = (pa - pt) / g('M_RTISA')
    p_out = max(pv, k_a * g('S_PREST'))
    p_n = (pa / g('CK_RUP') + p_out / g('CK_RA') - q_ivs) / (1 / g('CK_RUP') + 1 / g('CK_RA'))
    assert p_n > p_ivs and p_n > p_out, 'CK_RUP is too high for the resting intervillous flow'
    q_an = (p_n - p_out) / g('CK_RA'); s = q_tis / (q_tis - q_an)
    over = dict(CK_MODE=1, CK_KA=k_a, CK_RSA=(p_n - p_ivs) / q_ivs, M_RTISA=g('M_RTISA') * s, M_RTISV=g('M_RTISV') * s)
    info = dict(p_node=p_n, p_out=p_out, p_ivs=p_ivs, pa=pa, pv=pv, q_ivs=q_ivs, q_anast=q_an, q_ut_art=q_ivs + q_an,
                ivs_share=q_ivs / (q_ivs + q_an), tissue_scale=s)
    return over, info


def steady_clark(placenta=0):
    """(overrides, resting state) with the Clark 2018 uterine network at the default parameters and Jongen's mother.
    placenta=1: with the Hill 1972 exchanger. Cached."""
    key = 'SC%d' % placenta
    if key not in _SS:
        P = params(H_MODE=placenta)
        over, _ = clark_parameters(steady_hill() if placenta else steady(), P)
        _SS[key] = (over, steady_state(params(H_MODE=placenta, **over)))
    return _SS[key]


def maternal_pa(S, P):
    """Maternal arterial pressure p_a,m (mmHg) of state S."""
    return (S[SI['VM_ART']] - P[IDX['M_VART0']]) / P[IDX['M_CART']]


def run(put, w_um=0.0, g=1.0, dt=DTS, du=DU, out_every=OE, mhr=0, pth=None, thorax=0, vaso=0, autotx=0, circ=0,
        placenta=0, uterus=0, **kw):
    """Model output (dict of traces every out_every*dt s), final state and ok flag for uterine pressure `put` sampled
    every `du` s. w_um is the B(2) cord-compression weight (p_um = w_um (p_ut - p_rest)); g multiplies all five
    effector gains k of A(23) (g = 1 is the published regulation). Other parameters by name; when given, the steady
    state is recomputed for them.
    mhr: maternal HR reflex (Randall 2019) mode, 0 = off (Jongen 2016, MHR fixed at 60 / M_T), 1 = open loop (MHR
    computed, maternal heart period unchanged), 2 = closed loop (maternal heart period = 60 / MHR). pth: maternal
    thoracic pressure (mmHg) sampled like put, e.g. for pushing; None = MR_PTH_REST. The SBP offset MR_SBP_OFF is
    recomputed whenever other parameters are given, so that rest stays a steady state.
    thorax: 1 = maternal thorax after Lu 2001 (MT_*; pth then also acts on the circulation), 0 = off. Rest is
    unchanged (set_thorax_rest).
    vaso: 1 = maternal vasomotor reflex (MV_*; Randall's sympathetic tone also drives tissue arterial resistance and
    venous unstressed volume; needs mhr > 0). autotx: 1 = contraction autotransfusion (MA_*). Rest is unchanged.
    circ: 1 = maternal heart, vessels and baroreflex after Ursino 1998 (MU_MODE 1; see core.c). It replaces Jongen's
    maternal hemodynamics, so mhr, thorax, vaso and autotx must be 0, and it has its own resting state
    (steady_ursino); pth is not used. Outputs mhr, pa_m (step mean of P_sa), q_co_m, pv_m, q_uta, q_utv then come
    from the Ursino circulation.
    placenta: 1 = placental O2 transfer after Hill, Power & Longo 1972 (H_MODE 1; see core.c) instead of A(16), with
    the membrane diffusing capacity H_DMC calibrated to Jongen's resting fetal arterial pO2. It changes the resting
    state (steady_hill) and can be combined with every other option and with parameter overrides.
    uterus: 1 = uterine arteries after Clark et al. 2018 (CK_MODE 1; see core.c) instead of Jongen's single R_uta:
    uterine, arcuate and radial arteries, then spiral arteries into the intervillous space in parallel with an
    arterio-venous bypass. Outputs q_ut_art (uterine artery flow), q_uta (intervillous inflow), q_anast (bypass) and
    p_node. The resting state keeps Jongen's pressures, cardiac output, intervillous flow and fetal state
    (clark_parameters). Not available with circ=1; combines with the other options and with parameter overrides
    (CK_KA sets the compression of the bypass)."""
    assert not (uterus and circ), 'the Clark 2018 uterine network is not wired into the Ursino 1998 circulation'
    if circ:
        assert not (mhr or thorax or vaso or autotx or kw), 'circ excludes the other maternal add-ons and overrides'
        return run_core(params(S_WUM=w_um, S_GAIN=g, MU_MODE=1, H_MODE=placenta), steady_ursino(placenta=placenta), put,
                        dt=dt, du=du, out_every=out_every)
    if placenta:
        kw = dict(kw, H_MODE=1)
    P = params(S_WUM=w_um, S_GAIN=g, MR_MODE=mhr, MT_MODE=thorax, MV_MODE=vaso, MA_MODE=autotx, **kw)
    assert not vaso or mhr, 'the vasomotor reflex needs the HR reflex (mhr > 0) for its sympathetic tone'
    if uterus and not {k for k in kw if not k.startswith(('MR_', 'MT_', 'MV_', 'MA_'))} - {'H_MODE'}:
        over, S0 = steady_clark(placenta)                           # default network: cached resting state
        for k, v in over.items():
            P[IDX[k]] = v
        if 'MR_SBP_OFF' not in kw:
            P[IDX['MR_SBP_OFF']] = P[IDX['MR_SW']] - maternal_pa(S0, P)
        S = S0.copy()
    elif not kw:
        S = steady().copy()
    else:
        jk = {k: v for k, v in kw.items() if not k.startswith(('MR_', 'MT_', 'MV_', 'MA_'))}
        if jk == {'H_MODE': 1}:
            S0 = steady_hill()
        else:
            S0 = steady_state(params(S_WUM=w_um, S_GAIN=g, MR_MODE=0, **jk)) if jk else steady()  # add-ons are off at rest
        if uterus:                                                  # network calibrated on the overridden rest
            over, _ = clark_parameters(S0, params(S_WUM=w_um, S_GAIN=g, MR_MODE=0, **jk))
            jk.update(over)
            for k, v in over.items():
                P[IDX[k]] = v
            S0 = steady_state(params(S_WUM=w_um, S_GAIN=g, MR_MODE=0, **jk))
        if 'MR_SBP_OFF' not in kw:
            P[IDX['MR_SBP_OFF']] = P[IDX['MR_SW']] - maternal_pa(S0, P)
        S = S0.copy()
    if mhr:
        set_reflex_rest(P, S)
        P[IDX['MV_TS0']] = S[SI['R_TS']]
    if thorax:
        set_thorax_rest(P, S)
    return run_core(P, S, put, dt=dt, du=du, pth=pth, out_every=out_every)


def simulate(put, w_um, g, nb, bin_s=2):
    """Model FHR averaged over each `bin_s`-second CTG bin (nb bins), or None if the run is non-finite.
    FHR = 60 / T, A(25)."""
    P = params(S_WUM=w_um, S_GAIN=g)
    o, _, ok = run_core(P, steady(), put, dt=DTS, du=DU, out_every=OE)
    if not ok:
        return None
    f = o['fhr'][:nb * int(bin_s / (DTS * OE))].reshape(nb, -1).mean(1)
    return f if np.isfinite(f).all() else None
