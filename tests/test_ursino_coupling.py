"""Checks of the Ursino 1998 maternal circulation coupled to the fetal model (MU_MODE 1, simulate.run(circ=1)).
Run from the repo root: python tests/test_ursino_coupling.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import simulate as S  # noqa: E402

PUT = S.contraction(300, 0.5, t_con=120)
ON, S_END, OK = S.run(PUT, circ=1)
OFF = S.run(PUT)[0]
T = ON['t']; REST = T < 100; PEAK = (T >= 140) & (T < 160); LATE = T >= 280


def test_runs_and_returns_to_rest():
    assert OK
    for k, tol in (('mhr', 1.0), ('pa_m', 1.5), ('q_uta', 0.3), ('v_ivs', 1.0)):
        assert abs(ON[k][LATE].mean() - ON[k][REST].mean()) < tol, k


def test_blood_volume_conserved():
    """Ursino's total blood volume plus the intervillous volume stays constant."""
    S0 = S.steady_ursino(); tot = lambda s: s[S.SI['MU_VT']] + s[S.SI['VM_IVS']]
    assert abs(tot(S_END) - tot(S0)) < 1e-6


def test_contraction_reduces_uterine_flow_and_reflex_responds():
    """Uterine inflow falls as with Jongen's mother; the squeezed intervillous blood raises maternal arterial
    pressure and the baroreflex slows the maternal heart."""
    assert ON['q_uta'][PEAK].mean() < 0.5 * ON['q_uta'][REST].mean()
    assert OFF['q_uta'][PEAK].mean() < 0.5 * OFF['q_uta'][REST].mean()
    assert ON['pa_m'][PEAK].mean() > ON['pa_m'][REST].mean() + 1
    assert ON['mhr'][PEAK].mean() < ON['mhr'][REST].mean() - 1


def test_jongen_default_unchanged():
    """MU_MODE 0 leaves the Ursino states untouched and the maternal heart rate fixed."""
    assert np.all(OFF['mhr'] == 60.0 / S.params()[S.IDX['M_T']])


def test_circ_excludes_other_maternal_modes():
    try:
        S.run(PUT, circ=1, mhr=2)
    except AssertionError:
        return
    raise AssertionError('circ=1 with mhr should be refused')


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
