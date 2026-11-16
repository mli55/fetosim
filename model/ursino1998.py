"""Python driver for model/ursino1998.c: Ursino 1998 closed-loop pulsatile circulation with carotid baroreflex.

    from fetosim.model import ursino1998 as U
    S = U.steady_state(U.params())                 # basal state (HR 72 beats/min)
    o, S, ok = U.run(U.params(), S, 20.0)          # 20 s of pulsatile traces

Parameters: model/parameters/maternal_circulation_ursino1998.json (all from the paper's Tables 1-3).
"""
import ctypes
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PARAMETER_SET = HERE / 'parameters/maternal_circulation_ursino1998.json'
ENTRIES = json.loads(PARAMETER_SET.read_text())['parameters']
NAMES = [e['name'] for e in ENTRIES]
IDX = {n: i for i, n in enumerate(NAMES)}
USTATE = ['PPA', 'FPA', 'PPP', 'PPV', 'PSA', 'FSA', 'PSP', 'PEV', 'PLA', 'VLV', 'PRA', 'VRV', 'XI', 'PT',
          'DELV', 'DERV', 'DRSP', 'DREP', 'DVSV', 'DVEV', 'DTS', 'DTV', 'VT']
SI = {n: i for i, n in enumerate(USTATE)}
OUT = ['psa', 'plv', 'vlv', 'fol', 'pra', 'T', 'rsp', 'rep', 'elv', 'erv', 'vusv', 'vuev', 'fcs', 'fes', 'fev',
       'psv', 'pev', 'psp', 'fsp', 'fep', 'ppa', 'pla', 'vrv', 'xi']
DT = 2e-4                       # RK4 step, s

_lib = None


def lib():
    global _lib
    if _lib is None:
        src = HERE / 'ursino1998.c'; tag = hashlib.sha256(src.read_bytes()).hexdigest()[:12]
        so = HERE / f'build/ursino1998_{tag}.so'
        if not so.exists():
            so.parent.mkdir(exist_ok=True)
            subprocess.run(['cc', '-O3', '-shared', '-fPIC', '-o', str(so), str(src), '-lm'], check=True)
        L = ctypes.CDLL(str(so)); d = ctypes.POINTER(ctypes.c_double)
        L.ursino_simulate.argtypes = [d, d, ctypes.c_int, ctypes.c_double, d, d, ctypes.c_int, ctypes.c_double,
                                      ctypes.c_int, d]
        L.ursino_simulate.restype = ctypes.c_int
        assert L.ursino_n_par() == len(NAMES) and L.ursino_n_state() == len(USTATE) and L.ursino_n_out() == len(OUT)
        _lib = L
    return _lib


def params(**kw):
    P = np.array([e['value'] for e in ENTRIES], float)
    for k, v in kw.items():
        P[IDX[k]] = v
    return P


def initial_state(P):
    """Rough basal guess; steady_state() relaxes it. Effector states start at the Table 1-2 basal values."""
    g = lambda n: P[IDX[n]]
    S = np.zeros(len(USTATE))
    for n, v in dict(PPA=15, FPA=85, PPP=12, PPV=8, PSA=92, FSA=85, PSP=87, PEV=6, PLA=6.5, VLV=120, PRA=5, VRV=120,
                     PT=90.6).items():
        S[SI[n]] = v
    x = np.log(3.24)                                    # ln(f_es - f_es,min + 1) at the basal f_es of about 4.9
    for s, gk in (('DELV', 'U_GELV'), ('DERV', 'U_GERV'), ('DRSP', 'U_GRSP'), ('DREP', 'U_GREP'), ('DVSV', 'U_GVSV'),
                  ('DVEV', 'U_GVEV'), ('DTS', 'U_GTS')):
        S[SI[s]] = g(gk) * x
    S[SI['DTV']] = 0.833 - g('U_T0') - S[SI['DTS']]
    S[SI['VT']] = g('U_VT')
    return S


def _ptr(a):
    return None if a is None else a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def run(P, S, t_total, dt=DT, pcs=None, inflow=None, du=0.01, out_dt=0.01):
    """Integrate t_total s from state S (copied). pcs (mmHg; used when U_PCS_MODE = 1) and inflow (ml/s into the
    blood volume) are arrays sampled every du s, or None. Returns (dict of traces every out_dt s, final state, ok)."""
    n_steps = int(round(t_total / dt)); oe = max(1, int(round(out_dt / dt)))
    nrec = (n_steps + oe - 1) // oe
    out = np.zeros(nrec * len(OUT)); S = np.array(S, float); P = np.ascontiguousarray(P, float)
    arr = [None if a is None else np.ascontiguousarray(a, float) for a in (pcs, inflow)]
    n_in = max([len(a) for a in arr if a is not None] or [0])
    bad = lib().ursino_simulate(_ptr(P), _ptr(S), n_steps, dt, _ptr(arr[0]), _ptr(arr[1]), n_in, du, oe, _ptr(out))
    o = out.reshape(nrec, len(OUT))
    res = {k: o[:, i] for i, k in enumerate(OUT)}; res['t'] = np.arange(nrec) * oe * dt
    return res, S, bad == 0


def beat_mean(o, keys=None):
    """Means over an integer number of heart periods at the end of a run (the paper's averaging, Results)."""
    b = np.flatnonzero(np.diff(np.floor(o['xi'])) > 0) + 1
    assert len(b) >= 2, 'fewer than two beats in the record'
    sl = slice(b[0], b[-1]); dur = o['t'][b[-1]] - o['t'][b[0]]
    m = {k: float(o[k][sl].mean()) for k in (keys or OUT) if k != 'xi'}
    m['hr'] = 60.0 * (len(b) - 1) / dur
    m['co'] = float(o['fol'][sl].mean()) * 60 / 1000            # l/min
    m['sv'] = m['co'] * 1000 / m['hr']
    m['sbp'] = float(o['psa'][sl].max()); m['dbp'] = float(o['psa'][sl].min())
    m['tpr'] = (m['psa'] - m['pra']) / float(o['fol'][sl].mean())   # total systemic resistance, mmHg s/ml
    return m


def steady_state(P, S=None, t_settle=300.0, **kw):
    """Run t_settle s (default 300 s, more than ten times the slowest effector time constant plus delay)."""
    S = initial_state(P) if S is None else S
    _, S, ok = run(P, S, t_settle, out_dt=1.0, **kw)
    assert ok, 'steady-state run became non-finite'
    return S
