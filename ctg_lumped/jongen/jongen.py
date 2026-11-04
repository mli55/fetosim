"""Python interface to jongen.c, the forward-Euler core of the Jongen et al. 2016 CTG model (Parts A and B).
Parameter values are those of the papers' tables (source tag per parameter in PARAMS); DERIVED values are marked.
Compiles jongen.c with the system C compiler on first use (into build/ next to this file).

    from ctg_lumped.jongen import jongen as J
    P = J.params(); S = J.steady_state(P)
    out, S_end, ok = J.run(P, S, put=J.contraction(600, 0.02), dt=0.02, du=0.02)   # dict of output arrays
"""
import ctypes
import hashlib
import math
import subprocess
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ETA1B, ETA1C = math.log(0.5), math.log(0.001 / 1.001)          # B(3) exact; tables print -0.69 / -6.91

PARAMS = [  # (name, value, source)
    ('F_TREF', 60 / 135, 'A-T1a'), ('F_VART0', 35.2, 'A-T1a'), ('F_VVEN0', 153.4, 'A-T1a'), ('F_VTIS0', 15.9, 'A-T1a'),
    ('F_VCER0', 2.6, 'A-T1a'), ('F_VVIL0', 35.9, 'A-T1a'), ('F_CART', 0.33, 'A-T1a'), ('F_CVEN', 21.91, 'A-T1a'),
    ('F_CTIS', 0.35, 'A-T1a'), ('F_CCER', 0.06, 'A-T1a'), ('F_CVIL', 2.40, 'A-T1a'), ('F_RTISA', 2.81, 'A-T1a'),
    ('F_RTISV', 0.15, 'A-T1a'), ('F_RCERA', 7.24, 'A-T1a'), ('F_RCERV', 0.38, 'A-T1a'), ('F_RUMA', 3.84, 'A-T1a'),
    ('F_RUMV', 2.56, 'A-T1a'), ('F_EMAX', 7.71, 'A-T1a'), ('F_EMIN', 0.26, 'A-T1a'), ('F_V0MAX', 0.0, 'A-T1a'),
    ('F_V0MIN', 5.8, 'A-T1a'),
    ('M_T', 60 / 80, 'A-T1b'), ('M_VART0', 721, 'A-T1b'), ('M_VVEN0', 3604, 'A-T1b'), ('M_VTIS0', 504, 'A-T1b'),
    ('M_CART', 3.9, 'A-T1b'), ('M_CVEN', 308.9, 'A-T1b'), ('M_CTIS', 6.4, 'A-T1b'), ('M_RTISA', 0.67, 'A-T1b'),
    ('M_RTISV', 0.04, 'A-T1b'), ('M_RUTA', 5.71, 'A-T1b'), ('M_RUTV', 1.43, 'A-T1b'), ('M_R0UTV', 0.358, 'A-T1b'),
    ('M_P0UTV', -2.5, 'A-T1b'), ('M_P1UTV', 25.0, 'A-T1b'), ('M_EMAX', 1.83, 'A-T1b'), ('M_EMIN', 0.06, 'A-T1b'),
    ('M_V0MAX', 0.0, 'A-T1b'), ('M_V0MIN', 44.0, 'A-T1b'),
    ('IVS_VMAX', 350.0, 'DERIVED: 2 x 175 ml (B p.2475)'), ('IVS_P0', 5.0, 'DERIVED: steady-state IVS transmural'),
    ('IVS_P1', 350.0 / (math.pi * 2.4), 'DERIVED: Vmax/(pi*Civs,max)'),
    ('U_R0UMV', 2.05, 'B-T1'), ('U_P0UMV', -17.0, 'B-T1'), ('U_P1UMV', 10.0, 'B-T1'),
    ('U_R0UMA', 2.34, 'B-T1'), ('U_P0UMA', -27.0, 'B-T1'), ('U_P1UMA', 16.0, 'B-T1'),
    ('O_DPLAC', 0.0904, 'B s2.2.2 rule (steady-state fetal arterial pO2 18.5 mmHg); table A-T2 prints 0.082, which gives 18.25 here'), ('O_DLUNG', 0.169, 'A-T2'), ('O_QTISF', 0.336, 'A-T2'), ('O_QCERF', 0.131, 'A-T2'),
    ('O_QTISM', 10.0, 'A-T2'), ('O_THCER', 5.0, 'A-T2'), ('O_THTIS', 10.0, 'A-T2'),
    ('O_THTISM', 1e-6, 'NOT STATED: maternal threshold irrelevant (B p.2475); constant uptake'),
    ('O_PAIR', 160.0, 'A-T2'), ('O_ALPHA', 1.34, 'A-T2'), ('O_BETA', 3.1e-5, 'A-T2'), ('O_HBF', 0.17, 'A-T2'),
    ('O_HBM', 0.12, 'A-T2'), ('O_C1F', 1.04e4, 'A-T2'), ('O_C1M', 2.34e4, 'A-T2'), ('O_C2F', 150.0, 'A-T2'),
    ('O_C2M', 150.0, 'A-T2'),
    ('G_CREF', 0.107, 'A-T3'), ('G_GAMMA', 0.5, 'A-T3'), ('G_TAURC', 10.0, 'A-T3'),
    ('G_ETA1B', ETA1B, 'B(3) exact (table -0.69)'), ('G_ETA2B', -0.06, 'A-T3'), ('G_RBMIN', -1.0, 'A-T3'),
    ('G_RBMAX', 1.0, 'A-T3'), ('G_TAURB', 8.0, 'A-T3 (0 = no low-pass, thesis Ch.4)'), ('G_YREFB', 45.0, 'A-T3'),
    ('G_ETA1C', ETA1C, 'B(3) exact (table -6.91)'), ('G_ETA2C', 0.52, 'A-T3'), ('G_RCMIN', -0.001, 'A-T3'),
    ('G_RCMAX', 1.0, 'A-T3'), ('G_TAURCH', 2.0, 'A-T3'), ('G_YREFC', 18.1, 'A-T3'),
    ('G_WBS', 0.25, 'A-T3'), ('G_WCS', 0.75, 'A-T3'), ('G_WBV', 0.25, 'A-T3'), ('G_WCV', 0.75, 'A-T3'),
    ('G_DTV', 0.2, 'A-T3'), ('G_TAUTV', 1.5, 'A-T3'), ('G_KTV', 2.4, 'A-T3'),
    ('G_DTS', 2.0, 'A-T3'), ('G_TAUTS', 2.0, 'A-T3'), ('G_KTS', 1.2, 'A-T3'),
    ('G_DE', 2.0, 'A-T3'), ('G_TAUE', 8.0, 'A-T3'), ('G_KE', 1.2, 'A-T3'),
    ('G_DR', 2.0, 'A-T3'), ('G_TAUR', 6.0, 'A-T3'), ('G_KR', 4.2, 'A-T3'),
    ('G_DV', 5.0, 'A-T3'), ('G_TAUV', 20.0, 'A-T3'), ('G_KV', -1.2, 'A-T3'),
    ('S_PREST', 15.0, 'A s2.4'), ('S_WUM', 0.0, 'B(2): 0 = Part A, 1 = Part B'),
    ('S_GAIN', 1.0, 'Addition (not in the papers): common multiplier of all effector gains k (1 = published)'),
    ('S_MRT', 0.0, 'Th 4.2 noise on R_tisa'), ('S_MRB', 0.0, 'Th 4.2 noise on r_b'), ('S_MEV', 0.0, 'Th 4.2 noise on e_v'),
]
NAMES = [p[0] for p in PARAMS]
IDX = {n: i for i, n in enumerate(NAMES)}
STATE = ['VM_ART', 'VM_VEN', 'VM_TIS', 'VM_IVS', 'MM_ART', 'MM_VEN', 'MM_TIS', 'MM_IVS', 'MM_H',
         'VF_ART', 'VF_VEN', 'VF_TIS', 'VF_CER', 'VF_VIL', 'MF_ART', 'MF_VEN', 'MF_TIS', 'MF_CER', 'MF_VIL', 'MF_H',
         'HM_PREV', 'HF_PREV', 'RCERA', 'RB', 'RC', 'XTV', 'XTS', 'XE', 'XR', 'XV'] + [f'PG{i}' for i in range(11)]
SI = {n: i for i, n in enumerate(STATE)}
OUT = ['fhr', 'ptm_a', 'po2_a', 'q_uta', 'q_utv', 'v_ivs', 'po2_tis', 'po2_cer', 'q_tisa', 'q_cera', 'q_uma', 'q_umv',
       'v_villi', 'r_b', 'r_c', 'p_ut', 'e_s', 'e_v', 'z_tv', 'z_ts', 'q_co']

_lib = None


def lib():
    global _lib
    if _lib is None:
        src = HERE / 'jongen.c'; tag = hashlib.sha256(src.read_bytes()).hexdigest()[:12]
        so = HERE / f'build/jongen_{tag}.so'
        if not so.exists():
            so.parent.mkdir(exist_ok=True)
            subprocess.run(['cc', '-O3', '-shared', '-fPIC', '-o', str(so), str(src), '-lm'], check=True)
        L = ctypes.CDLL(str(so)); d = ctypes.POINTER(ctypes.c_double)
        L.simulate.argtypes = [d, d, ctypes.c_int, ctypes.c_double, d, d, d, d, ctypes.c_int, ctypes.c_double, ctypes.c_int, d]
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
    """DERIVED steady-state guess (model_spec 6.6); exact steady state from steady_state()."""
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
    return S


def run(P, S, put, dt=0.02, du=None, pum=None, ruta=None, noise=None, out_every=None, copy_state=True):
    """Run len(put)*du/dt steps from state S (copied unless copy_state=False). Returns (outputs dict, final state, ok)."""
    du = dt if du is None else du
    put = np.ascontiguousarray(put, float); n_in = len(put)
    n_steps = int(round((n_in - 1) * du / dt)) + 1
    out_every = out_every or max(1, int(round(0.5 / dt)))
    nrec = (n_steps + out_every - 1) // out_every
    out = np.zeros(nrec * len(OUT)); S = np.array(S, float) if copy_state else S
    arr = [None if a is None else np.ascontiguousarray(a, float) for a in (pum, ruta, noise)]
    for a in arr:
        assert a is None or len(a) == n_in
    P = np.ascontiguousarray(P, float)
    bad = lib().simulate(_ptr(P), _ptr(S), n_steps, dt, _ptr(put), _ptr(arr[0]), _ptr(arr[1]), _ptr(arr[2]),
                         n_in, du, out_every, _ptr(out))
    o = out.reshape(nrec, len(OUT))
    res = {k: o[:, i] for i, k in enumerate(OUT)}; res['t'] = np.arange(nrec) * out_every * dt
    return res, S, bad == 0


def steady_state(P, dt=0.02, t_max=3000.0, tol=1e-4):
    """Run at rest until |d p_tm,a/dt| and |d pO2,a/dt| < tol mmHg/s (A s2.4) over the last 10 s."""
    S = initial_state(P); rest = np.full(int(100 / dt) + 1, P[IDX['S_PREST']]); t = 0.0
    while t < t_max:
        o, S, ok = run(P, S, rest, dt=dt, out_every=int(round(1 / dt)))
        assert ok, 'steady-state run became non-finite'
        t += 100
        rate = max(abs(o['ptm_a'][-1] - o['ptm_a'][-11]) / 10, abs(o['po2_a'][-1] - o['po2_a'][-11]) / 10)
        if rate < tol:
            return S
    raise RuntimeError(f'no steady state within {t_max} s (rate {rate:.2e})')


def contraction(t_total, dt, t_con=30.0, T_con=60.0, p_con=70.0, p_rest=15.0):
    """A(8) single sin^2 contraction sampled every dt."""
    t = np.arange(int(round(t_total / dt)) + 1) * dt
    p = np.full(len(t), p_rest); m = (t > t_con) & (t < t_con + T_con)
    p[m] += p_con * np.sin(np.pi * (t[m] - t_con) / T_con) ** 2
    return p
