"""Checks of the Ursino 1998 circulation + carotid baroreflex (model/ursino1998.c) against the paper. The full
comparison lives in the ARPH project (research_log/experiments/ursino1998_replication_2026-11-15).
Run from the repo root: python tests/test_ursino1998.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import ursino1998 as U  # noqa: E402

P0 = U.params(); S0 = U.steady_state(P0); B = U.beat_mean(U.run(P0, S0, 20.0)[0])


def test_basal_matches_tables():
    """The closed loop settles on the basal values of Tables 1-2 (which the paper lists separately from Table 3)."""
    for k, v, tol in (('hr', 72, 0.02), ('rsp', 3.307, 0.01), ('rep', 1.407, 0.01), ('elv', 2.95, 0.01),
                      ('erv', 1.75, 0.01), ('vusv', 1121, 0.01), ('vuev', 1375, 0.01)):
        assert abs(B[k] / v - 1) < tol, (k, B[k], v)
    assert abs(B['fsp'] * 0.06 - 1.6) < 0.05                     # splanchnic flow 1.6 l/min (Methods)
    assert abs(B['sbp'] - 131) < 3 and abs(B['dbp'] - 88) < 3    # Fig. 2


def test_blood_volume_conserved():
    o, S, ok = U.run(P0, S0, 60.0)
    assert ok and S[U.SI['VT']] == P0[U.IDX['U_VT']] and abs(U.beat_mean(o)['psa'] - B['psa']) < 0.2


def test_pacing_curve():
    """Fig. 6: cardiac output 82 ml/s at 60 and 111 ml/s at 120 beats/min; stroke volume falls with rate."""
    m = {hr: U.beat_mean(U.run(P, U.steady_state(P, S0), 20.0)[0])
         for hr in (60, 120) for P in [U.params(U_PACE_T=60.0 / hr)]}
    assert abs(m[60]['co'] * 1000 / 60 - 82) < 2 and abs(m[120]['co'] * 1000 / 60 - 111) < 2
    assert m[120]['sv'] < m[60]['sv']


def test_hemorrhage():
    """Fig. 9: 10 % blood loss, 80 s later: SAP -10 %, CO -16 %, HR +8.7 %."""
    t = np.arange(0, 110, 0.01); I = np.where((t >= 10) & (t < 15), -530.0 / 5, 0.0)
    o, S, ok = U.run(P0, S0, 107.0, inflow=I)
    sel = o['t'] >= 85; m = U.beat_mean({k: v[sel] for k, v in o.items()})
    assert ok and abs(S[U.SI['VT']] - 4770) < 1e-6
    for k, pct in (('psa', -10), ('co', -16), ('hr', 8.7)):
        assert abs(100 * (m[k] / B[k] - 1) - pct) < 1.5, (k, m[k], B[k])


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
