"""Checks of the uterine arterial network after Clark et al. 2018 (CK_MODE 1, S.run(uterus=1)). The port, its
calibration and the comparison of coupling choices live in the ARPH project
(research_log/experiments/clark_port_2026-11-20). Run from the repo root: python tests/test_clark2018.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import simulate as S  # noqa: E402


def _contraction(uterus, w_um=0.0, **kw):
    c = {k: kw.pop(k) for k in ('p_con', 'T_con', 't_con') if k in kw}
    o, _, ok = S.run(S.contraction(300, 0.02, **c), w_um=w_um, du=0.02, uterus=uterus, **kw)
    assert ok
    return o


def test_network_resistances():
    """Clark's baseline geometry: 0.9027 and 1.6841 Pa s/mm3 for one side; x 7.5006 / 2 in the parameter file."""
    P = S.params()
    assert np.isclose(P[S.IDX['CK_RUP']], 0.90265 * 7.50064 / 2, rtol=1e-4)
    assert np.isclose(P[S.IDX['CK_RA']], 1.68408 * 7.50064 / 2, rtol=1e-4)


def test_rest_unchanged():
    put = np.full(241, S.TONE)
    for kw in (dict(), dict(mhr=2, thorax=1, vaso=1, autotx=1), dict(placenta=1)):
        on, _, ok1 = S.run(put, uterus=1, **kw); off, _, ok0 = S.run(put, **kw)
        assert ok1 and ok0
        for k in ('fhr', 'mhr', 'pa_m', 'pv_m', 'q_co_m', 'po2_a', 'q_uta', 'v_ivs'):
            assert np.abs(on[k] - off[k]).max() < 2e-3 * max(1.0, np.abs(off[k]).max()), (kw, k)
        assert np.allclose(on['q_ut_art'], on['q_uta'] + on['q_anast'])
        assert np.allclose(off['q_ut_art'], off['q_uta']) and not off['q_anast'].any()
    _, info = S.clark_parameters(S.steady(), S.params())
    assert 800 < 60 * info['q_ut_art'] < 900 and 0.70 < info['ivs_share'] < 0.80


def test_rest_unchanged_with_other_tone():
    """The network is calibrated on the resting state of the overridden parameters (here another resting pressure)."""
    put = np.full(241, 10.0)
    on, _, _ = S.run(put, uterus=1, S_PREST=10.0); off, _, _ = S.run(put, S_PREST=10.0)
    for k in ('fhr', 'pa_m', 'po2_a', 'q_uta'):
        assert np.abs(on[k] - off[k]).max() < 2e-3 * np.abs(off[k]).max(), k


def test_jongen_scenarios_kept():
    """Published FHR values of Jongen 2016 (deceleration Part A: 3, 14, 21 bpm; minimum Part B: 85 bpm)."""
    dec = lambda o: o['fhr'][0] - o['fhr'].min()
    assert abs(dec(_contraction(1)) - 3) < 0.5
    assert abs(dec(_contraction(1, p_con=110)) - 14) < 1.0
    assert abs(dec(_contraction(1, T_con=120, t_con=0)) - 21) < 1.0
    b = _contraction(1, w_um=1.0)
    assert abs(b['fhr'].min() - 85) < 1.0
    assert abs(b['fhr'].min() - _contraction(0, w_um=1.0)['fhr'].min()) < 0.2


def test_flow_falls_in_contraction():
    o = _contraction(1)
    drop = lambda k: 1 - o[k].min() / o[k][0]
    assert 0.55 < drop('q_ut_art') < 0.80 and 0.60 < drop('q_uta') < 0.85
    assert o['q_anast'].min() >= 0 and o['q_anast'].max() <= o['q_anast'][0] * 1.01


def test_uncompressed_bypass_steals_flow():
    """CK_KA 0: the bypass is not compressed, takes more flow during the contraction and deepens the deceleration."""
    o = _contraction(1, CK_KA=0.0)
    assert o['q_anast'].max() > 1.5 * o['q_anast'][0]
    assert o['fhr'][0] - o['fhr'].min() > 2 * (_contraction(1)['fhr'][0] - _contraction(1)['fhr'].min())


def test_maternal_blood_volume_conserved():
    over, S0 = S.steady_clark()
    _, S1, ok = S.run_core(S.params(**over), S0, S.contraction(300, 0.5), du=0.5)
    vol = lambda s: sum(s[S.SI[k]] for k in ('VM_ART', 'VM_VEN', 'VM_TIS', 'VM_IVS', 'HM_PREV'))
    assert ok and abs(vol(S1) - vol(S0)) < 1e-6


def test_not_with_ursino_circulation():
    try:
        S.run(np.full(11, S.TONE), circ=1, uterus=1)
    except AssertionError:
        return
    raise AssertionError('uterus=1 with circ=1 should be refused')


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
