"""Checks of the maternal vasomotor reflex (MV_MODE 1) and contraction autotransfusion (MA_MODE 1). Calibration and the
comparison with data live in the ARPH project (research_log/experiments/maternal_vasomotor_autotransfusion_2026-11-07).
Run from the repo root: python tests/test_vasomotor.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import simulate as S  # noqa: E402


def test_rest_unchanged():
    put = np.full(1201, S.TONE)
    on, _, ok1 = S.run(put, du=0.1, mhr=2, thorax=1, vaso=1, autotx=1)
    off, _, ok0 = S.run(put, du=0.1)
    assert ok1 and ok0
    for k in ('mhr', 'pa_m', 'q_co_m', 'fhr', 'po2_a'):
        assert np.abs(on[k] - off[k]).max() < 1e-6, k


def test_valsalva_overshoot():
    """With the vasomotor reflex, release is followed by a pressure overshoot and a heart-rate fall (phase 4)."""
    t = np.arange(1201) * 0.1; dp = np.where((t >= 60) & (t < 75), 40.0, 0.0)
    kw = dict(du=0.1, mhr=2, thorax=1, pth=4.5 + dp)
    on, _, ok1 = S.run(S.TONE + dp, vaso=1, **kw)
    off, _, ok0 = S.run(S.TONE + dp, vaso=0, **kw)
    assert ok1 and ok0
    c = (on['t'] >= 55) & (on['t'] < 60); ph4 = (on['t'] >= 79) & (on['t'] < 95)
    assert on['pa_m'][ph4].max() > on['pa_m'][c].mean() + 10
    assert on['pa_m'][ph4].max() > off['pa_m'][ph4].max() + 10
    assert on['mhr'][ph4].min() < on['mhr'][c].mean() - 5
    assert on['rtisa_m'][(on['t'] >= 70) & (on['t'] < 75)].mean() > 1.1 * on['rtisa_m'][c].mean()


def test_autotransfusion_raises_stroke_volume():
    put = S.contraction(300, 0.5, t_con=120)
    on, _, ok1 = S.run(put, mhr=2, thorax=1, autotx=1)
    off, _, ok0 = S.run(put, mhr=2, thorax=1, autotx=0)
    assert ok1 and ok0
    sv = lambda o: o['q_co_m'] * 60 / o['mhr']
    pk = (on['t'] >= 140) & (on['t'] < 160); b = on['t'] < 100
    assert sv(on)[pk].mean() > 1.1 * sv(on)[b].mean()
    assert sv(on)[pk].mean() > sv(off)[pk].mean() + 5
    assert np.isclose(on['v_sh_m'][pk].max(), S.params()[S.IDX['MA_VSQ']], rtol=0.05)


def test_vaso_needs_hr_reflex():
    try:
        S.run(np.full(11, S.TONE), du=0.1, mhr=0, vaso=1)
    except AssertionError:
        return
    raise AssertionError('vaso=1 without mhr should be refused')


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
