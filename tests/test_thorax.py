"""Checks of the maternal thorax graft (MT_MODE 1, after Lu et al. 2001). Calibration and the comparison with data
live in the ARPH project (research_log/experiments/maternal_thorax_graft_2026-11-07).
Run from the repo root: python tests/test_thorax.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import simulate as S  # noqa: E402


def test_rest_unchanged():
    put = np.full(1201, S.TONE)                                     # 120 s at 0.1 s
    on, _, ok1 = S.run(put, du=0.1, mhr=2, thorax=1)
    off, _, ok0 = S.run(put, du=0.1, mhr=2, thorax=0)
    assert ok1 and ok0
    for k in ('mhr', 'pa_m', 'q_co_m', 'fhr', 'po2_a'):
        assert np.abs(on[k] - off[k]).max() < 1e-6, k


def test_valsalva_direction():
    t = np.arange(1201) * 0.1; dp = np.where((t >= 60) & (t < 75), 40.0, 0.0)
    o, _, ok = S.run(S.TONE + dp, du=0.1, mhr=2, thorax=1, pth=4.5 + dp)
    assert ok
    c = (o['t'] >= 55) & (o['t'] < 60); late = (o['t'] >= 74) & (o['t'] < 75)
    assert o['q_co_m'][late].mean() < 0.6 * o['q_co_m'][c].mean()     # venous return falls
    assert o['mhr'][late].mean() > o['mhr'][c].mean() + 10            # heart rate rises during the strain
    assert o['sbp_m'][late].mean() < o['sbp_m'][c].mean()             # phase-2 pressure below control


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
