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


def test_with_ursino_circulation_rest():
    """With the Ursino 1998 mother the network keeps the resting means of Ursino mode (arterial pressure, intervillous
    inflow, fetal arterial pO2 and FHR) and adds the bypass."""
    put = np.full(241, S.TONE)
    m = {}
    for u in (0, 1):
        o, _, ok = S.run(put, du=0.5, circ=1, uterus=u, out_every=1); assert ok
        w = o['t'] >= 60; m[u] = {k: float(np.mean(o[k][w])) for k in ('pa_m', 'q_uta', 'q_anast', 'po2_a', 'fhr')}
    assert abs(m[1]['pa_m'] - m[0]['pa_m']) < 0.1
    assert abs(m[1]['q_uta'] / m[0]['q_uta'] - 1) < 0.005
    assert abs(m[1]['po2_a'] - m[0]['po2_a']) < 0.1 and abs(m[1]['fhr'] - m[0]['fhr']) < 0.05
    assert m[0]['q_anast'] == 0 and m[1]['q_anast'] > 2


def test_with_ursino_circulation_contraction():
    """With the Ursino mother a contraction lowers uterine artery flow and intervillous inflow, and the bypass is
    squeezed but stays open (k_A p_ut below the node pressure at 85 mmHg is not required; flow only has to fall)."""
    put = S.contraction(300, 0.5, t_con=120)
    o, _, ok = S.run(put, du=0.5, circ=1, uterus=1, out_every=50); assert ok
    rest, con = o['t'] < 110, (o['t'] > 120) & (o['t'] < 200)
    for k in ('q_ut_art', 'q_uta', 'q_anast'):
        assert o[k][con].min() < 0.7 * o[k][rest].mean(), k


def test_rlc_parameters_add_up():
    """CK_MODE 2 uses the same resistances as CK_MODE 1, split by vessel class."""
    P = S.params(); g = lambda n: P[S.IDX[n]]
    assert abs(g('CK_RU') + g('CK_RAR') + g('CK_RR') - g('CK_RUP')) < 1e-9
    assert abs(g('CK_RAV') + g('CK_RD') - g('CK_RA')) < 1e-9


def test_rlc_rest_matches_resistive_network():
    """With inertance and compliance the mean resting flows, pressure and fetal state equal those of the resistive
    network (the network is linear while the bypass is open, so L and C change only the waveform)."""
    put = np.full(241, S.TONE); m = {}
    for u in (1, 2):
        o, _, ok = S.run(put, du=0.5, circ=1, uterus=u, out_every=1); assert ok
        w = o['t'] >= 60; m[u] = {k: float(np.mean(o[k][w])) for k in ('pa_m', 'q_ut_art', 'q_uta', 'q_anast', 'fhr')}
    for k in ('q_ut_art', 'q_uta', 'q_anast'):
        assert abs(m[2][k] / m[1][k] - 1) < 0.01, k
    assert abs(m[2]['pa_m'] - m[1]['pa_m']) < 0.3 and abs(m[2]['fhr'] - m[1]['fhr']) < 0.05


def test_rlc_doppler_matches_clark_baseline():
    """Resting Doppler indices at the insonation site are close to the authors' code at baseline geometry
    (network-wave-transmission RUNME.py: PI 0.527, RI 0.392, S/D 1.645, no notch), although the inlet here is Ursino's
    pressure waveform instead of their prescribed flow waveform."""
    o, _, ok = S.run(np.full(241, S.TONE), du=0.5, circ=1, uterus=2, out_every=1); assert ok
    d = S.doppler_indices(o, 60, 120)
    assert len(d['PI']) > 30 and not d['notch'].any()
    assert abs(np.median(d['PI']) - 0.527) < 0.08 and abs(np.median(d['RI']) - 0.392) < 0.05
    assert abs(np.median(d['SD']) - 1.645) < 0.15


def test_rlc_blood_volume_conserved():
    """Ursino's total volume + intervillous volume + volume stored in the network compliances stays constant
    through a contraction."""
    over, S0 = S.steady_clark_ursino(rlc=1); P = S.params(MU_MODE=1, **over); g = lambda n: P[S.IDX[n]]
    C = [g('CK_CU') / 5] * 5 + [g('CK_CAR'), g('CK_CR'), g('CK_CAV')]
    vol = lambda X: X[S.SI['MU_VT']] + X[S.SI['VM_IVS']] + sum(c * X[S.SI['CK_P%d' % i]] for i, c in enumerate(C))
    _, S1, ok = S.run_core(P, S0, S.contraction(200, 0.5), du=0.5)
    assert ok and abs(vol(S1) - vol(S0)) < 1e-6


def test_rlc_needs_ursino_mother():
    try:
        S.run(np.full(11, S.TONE), uterus=2)
    except AssertionError:
        return
    raise AssertionError('uterus=2 without circ=1 should be refused')


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
