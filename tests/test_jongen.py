"""Quick checks against values printed in Jongen et al. 2016 (run: python tests/test_jongen.py, or pytest tests)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # repo root, no install needed

from jongen2016 import model as M  # noqa: E402


def _contraction(w_um):
    o, _, ok = M.run(M.contraction(300, 0.02), w_um=w_um, g=1.0, du=0.02)
    assert ok
    return o


def test_steady_state():
    o = _contraction(0.0)
    assert abs(o['fhr'][0] - 135) < 1                    # A-T1a: T_ref = 60/135 s
    assert abs(o['po2_a'][0] - 18.5) < 0.1               # B s2.2.2: steady-state fetal arterial pO2


def test_part_a_late_deceleration():
    o = _contraction(0.0)
    assert abs((o['fhr'][0] - o['fhr'].min()) - 3) < 1   # A: reference contraction gives about 3 bpm


def test_part_b_variable_deceleration():
    o = _contraction(1.0)
    assert abs(o['fhr'].min() - 85) < 1.5                # B Fig. 1: nadir about 85 bpm


def test_fit_runs_on_synthetic_window():
    rng = np.random.default_rng(0)
    nb = 400; t = (np.arange(nb) + 1) * M.DT
    toco = 10 + sum(40 * np.exp(-0.5 * ((t - c) / 20) ** 2) for c in range(120, 800, 180)) + rng.normal(0, 1, nb)
    fhr = 140 + rng.normal(0, 2, nb)
    r = M.fit(fhr, toco, train=300)
    assert len(r['peaks_s']) >= 3 and np.isfinite(r['fhr_model']).all()
    assert M.A_BOUNDS[0] <= r['a'] <= M.A_BOUNDS[1] and M.W_BOUNDS[0] <= r['w_um'] <= M.W_BOUNDS[1]


def test_parameter_classes():
    rows = M.parameter_table()
    fitted = sorted(r['symbol'] for r in rows if r['cls'] == 3)
    assert fitted == ['a', 'b', 'g', 'w_um']


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
