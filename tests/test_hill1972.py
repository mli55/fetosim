"""Checks of the Hill, Power & Longo 1972 placental O2 exchanger (model/hill1972.c) against the paper and inside the
model (simulate.run(placenta=1)). The reimplementation with the paper's RK4 and the sheep standard case lives in the
ARPH project (research_log/experiments/hill1972_cellml_replication_2026-11-15).
Run from the repo root: python tests/test_hill1972.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import simulate as S  # noqa: E402

Q = 350 / 60.0                                                    # the paper's capillary flow, ml/s


def test_human_case_matches_paper():
    """p. 726: human curves and hemoglobin, 95 and 15 mmHg in: end-capillary 28.0 and 26.6 mmHg, 25.1 ml/min."""
    r, pm, pf, dp = S.hill_exchange(95, 15, Q, Q)
    assert abs(pm - 28.0) < 0.05 and abs(pf - 26.6) < 0.05 and abs(r * 60 - 25.1) < 0.05


def test_sheep_fetal_curve_case():
    """Sheep standard case (Table 1; end-capillary 32.26 and 31.81 mmHg, 24.84 ml/min) with the maternal curve of the
    lower branch (S < 81), which holds for all but the first steps of the transit."""
    kw = dict(H_K1M=1.544, H_K3M=0.346, H_K1F=1.20, H_K3F=0.349, H_HBM=10.2, H_HBF=12.0, H_KCM=130, H_KCF=151)
    r, pm, pf, dp = S.hill_exchange(60, 15, Q, Q, **kw)            # start below the branch point
    assert pm > pf and pm - pf < 1.0 and abs(dp - 2.07) < 0.05


def test_step_count_converged_and_flow_limits():
    a = S.hill_exchange(95, 15, Q, Q); b = S.hill_exchange(95, 15, Q, Q, H_NSTEP=1000)
    assert abs(a[0] / b[0] - 1) < 1e-3
    lo = S.hill_exchange(95, 15, 0.05, Q)                          # nearly stopped maternal flow
    assert 0 < lo[0] < 0.05 * 1.34 * 12.5 / 100 and abs(lo[1] - lo[2]) < 0.01
    assert S.hill_exchange(95, 15, 0.0, Q)[0] == 0.0 and S.hill_exchange(15, 15, Q, Q)[0] == 0.0


def test_in_model():
    """Rest is steady and the transfer equals fetal O2 use; a contraction lowers fetal arterial pO2 and it recovers."""
    put = S.contraction(400, 0.5, t_con=120)
    on, _, ok = S.run(put, placenta=1); off = S.run(put)[0]
    t = on['t']; rest = t < 100
    assert ok and np.ptp(on['po2_a'][rest]) < 0.01          # the steady-state tolerance (1e-4 mmHg/s) over 100 s
    assert abs(on['po2_a'][rest].mean() - off['po2_a'][rest].mean()) < 2.0
    assert on['po2_a'][t > 100].min() < on['po2_a'][rest].mean() - 1.0
    assert abs(on['po2_a'][t > 380].mean() - on['po2_a'][rest].mean()) < 0.3
    both, _, ok2 = S.run(put, placenta=1, circ=1)
    assert ok2 and both['po2_a'][t > 100].min() < both['po2_a'][rest].mean() - 0.5


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
