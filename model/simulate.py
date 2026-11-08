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


def load_parameter_set(path=PARAMETER_SET):
    """List of (name, value, source) in core.c order, plus the full entries. JSON floats round-trip exactly (repr),
    so the values are bit-identical to the ones written."""
    entries = json.loads(Path(path).read_text())['parameters']
    return [(e['name'], e['value'], e['source_detail']) for e in entries], entries


PARAMS, ENTRIES, NAMES, IDX = [], [], [], {}


def reload_parameters():
    """(Re)read the two parameter files into PARAMS, ENTRIES, NAMES, IDX and drop the cached steady state."""
    PARAMS.clear(); ENTRIES.clear(); NAMES.clear(); IDX.clear()
    for path in (PARAMETER_SET, MATERNAL_HR_SET):
        p, e = load_parameter_set(path); PARAMS.extend(p); ENTRIES.extend(e)
    NAMES.extend(p[0] for p in PARAMS); IDX.update({n: i for i, n in enumerate(NAMES)})
    _SS.clear()


_SS = {}
reload_parameters()
STATE = ['VM_ART', 'VM_VEN', 'VM_TIS', 'VM_IVS', 'MM_ART', 'MM_VEN', 'MM_TIS', 'MM_IVS', 'MM_H',
         'VF_ART', 'VF_VEN', 'VF_TIS', 'VF_CER', 'VF_VIL', 'MF_ART', 'MF_VEN', 'MF_TIS', 'MF_CER', 'MF_VIL', 'MF_H',
         'HM_PREV', 'HF_PREV', 'RCERA', 'RB', 'RC', 'XTV', 'XTS', 'XE', 'XR', 'XV'] + [f'PG{i}' for i in range(11)] + \
        ['R_EBC', 'R_EBA', 'R_TPB', 'R_TPR', 'R_TS', 'R_H']
SI = {n: i for i, n in enumerate(STATE)}
OUT = ['fhr', 'ptm_a', 'po2_a', 'q_uta', 'q_utv', 'v_ivs', 'po2_tis', 'po2_cer', 'q_tisa', 'q_cera', 'q_uma', 'q_umv',
       'v_villi', 'r_b', 'r_c', 'p_ut', 'e_s', 'e_v', 'z_tv', 'z_ts', 'q_co',
       'mhr', 'pa_m', 'q_co_m', 'sbp_m', 'pth_m']
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
        src = HERE / 'core.c'; tag = hashlib.sha256(src.read_bytes()).hexdigest()[:12]
        so = HERE / f'build/core_{tag}.so'
        if not so.exists():
            so.parent.mkdir(exist_ok=True)
            subprocess.run(['cc', '-O3', '-shared', '-fPIC', '-o', str(so), str(src), '-lm'], check=True)
        L = ctypes.CDLL(str(so)); d = ctypes.POINTER(ctypes.c_double)
        L.simulate.argtypes = [d, d, ctypes.c_int, ctypes.c_double, d, d, d, d, d, ctypes.c_int, ctypes.c_double, ctypes.c_int, d]
        L.simulate.restype = ctypes.c_int
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


def run_core(P, S, put, dt=0.02, du=None, pum=None, ruta=None, noise=None, pth=None, out_every=None, copy_state=True):
    """Low-level run with a full parameter vector P: len(put)*du/dt steps from state S (copied unless copy_state=False).
    pth is the maternal thoracic pressure for the HR reflex (None = MR_PTH_REST). Returns (outputs dict, final state, ok)."""
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
    return res, S, bad == 0


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


def maternal_pa(S, P):
    """Maternal arterial pressure p_a,m (mmHg) of state S."""
    return (S[SI['VM_ART']] - P[IDX['M_VART0']]) / P[IDX['M_CART']]


def run(put, w_um=0.0, g=1.0, dt=DTS, du=DU, out_every=OE, mhr=0, pth=None, **kw):
    """Model output (dict of traces every out_every*dt s), final state and ok flag for uterine pressure `put` sampled
    every `du` s. w_um is the B(2) cord-compression weight (p_um = w_um (p_ut - p_rest)); g multiplies all five
    effector gains k of A(23) (g = 1 is the published regulation). Other parameters by name; when given, the steady
    state is recomputed for them.
    mhr: maternal HR reflex (Randall 2019) mode, 0 = off (Jongen 2016, MHR fixed at 60 / M_T), 1 = open loop (MHR
    computed, maternal heart period unchanged), 2 = closed loop (maternal heart period = 60 / MHR). pth: maternal
    thoracic pressure (mmHg) sampled like put, e.g. for pushing; None = MR_PTH_REST. The SBP offset MR_SBP_OFF is
    recomputed whenever other parameters are given, so that rest stays a steady state."""
    P = params(S_WUM=w_um, S_GAIN=g, MR_MODE=mhr, **kw)
    if not kw:
        S = steady().copy()
    else:
        S0 = steady_state(params(S_WUM=w_um, S_GAIN=g, MR_MODE=0, **kw))
        if 'MR_SBP_OFF' not in kw:
            P[IDX['MR_SBP_OFF']] = P[IDX['MR_SW']] - maternal_pa(S0, P)
        S = S0.copy()
    if mhr:
        set_reflex_rest(P, S)
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
