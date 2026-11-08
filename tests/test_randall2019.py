"""Self-contained checks of the Randall 2019 heart-rate port (model/randall2019.py). The full reproduction against the
authors' Fortran and data lives in the ARPH project (research_log/experiments/randall2019_reproduction_2026-11-07).
Run from the repo root: python tests/test_randall2019.py, or pytest tests."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.model import randall2019 as R  # noqa: E402

# resting values of the authors' representative control subject (paper "control subject 2")
REST = dict(Pbar=148.85822826, Pthbar=4.49329347, Hbar=96.96345644, HminR=80.9879, HmaxR=112.4263, age=21.0)


def test_rest_stays_at_resting_heart_rate():
    p = R.nominal_parameters(**REST)
    t = np.arange(0, 60.01, 0.1)
    P = np.full(len(t), REST['Pbar']); Pth = np.full(len(t), REST['Pthbar'])
    H = R.heart_rate(t, P, Pth, p, R.initial_conditions(p, P[0], Pth[0], REST['Hbar']), h=1e-3)
    assert np.abs(H - REST['Hbar']).max() < 1e-6     # nominal parameters are built so that rest is a steady state


def test_spline_matches_clamped_cubic_spline():
    from scipy.interpolate import CubicSpline
    x = np.linspace(0, 10, 41); y = np.sin(x) + 0.1 * x ** 2
    xs = np.linspace(-1, 11, 500)
    ref = CubicSpline(x, y, bc_type=((1, 0.0), (1, 0.0)))(np.clip(xs, x[0], x[-1]))
    assert np.abs(R.nr_splint(x, y, R.nr_spline(x, y), xs) - ref).max() < 1e-10


def test_baroreflex_direction():
    """A step rise in SBP slows the heart (vagal) and a step fall speeds it up."""
    p = R.nominal_parameters(**REST)
    t = np.arange(0, 40.01, 0.1)
    out = {}
    for d in (+20, -20):
        P = np.where(t < 5, REST['Pbar'], REST['Pbar'] + d); Pth = np.full(len(t), REST['Pthbar'])
        H = R.heart_rate(t, P, Pth, p, R.initial_conditions(p, P[0], Pth[0], REST['Hbar']), h=1e-3)
        out[d] = H[-1] - REST['Hbar']
    assert out[+20] < -1 and out[-20] > 1


def test_fetosim_maternal_hr_coupling():
    """Rest stays at MHR = 60 / M_T in closed loop; a reference contraction lowers MHR (baroreflex); mode 0 is Jongen."""
    from fetosim.model import simulate as S
    rest = np.full(int(120 / 0.02) + 1, 15.0)
    o, _, ok = S.run(rest, mhr=2, du=0.02)
    assert ok and np.abs(o['mhr'] - 80).max() < 1e-9
    put = S.contraction(300, 0.02)
    o0, _, _ = S.run(put, mhr=0, du=0.02); o1, _, _ = S.run(put, mhr=1, du=0.02)
    assert np.all(o0['mhr'] == 80) and np.array_equal(o0['fhr'], o1['fhr'])   # open loop leaves the circulation alone
    assert -8 < o1['mhr'].min() - 80 < -2


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f(); print('ok', name)
