"""Reproduce the published Jongen et al. 2016 scenarios at the published regulation (g = 1).

Part A (w_um = 0): uterine flow reduction only -> late deceleration (reference contraction, 110 mmHg, 120 s,
uterine-artery block of 20 / 30 s). Part B (w_um = 1): + umbilical cord compression -> variable deceleration
(reference, 110 / 30 mmHg, 120 / 30 s, cord occlusion at 25 / 33 / 41 / 90 mmHg).
Published values: text values where the paper states them, otherwise figure readings (about +-2 bpm).
Writes outputs/replication.json (relative to the current folder) and prints simulated vs published.

    python scripts/replicate_jongen.py

Known differences from the papers:
  - D_plac = 0.0904 (paper's rule: steady-state fetal arterial pO2 18.5 mmHg) instead of the table's 0.082.
  - Part A reference: minimum uterine flow 163 ml/min here vs about 180 in A Fig. 3 (figure reading).
  - Contraction onset t_con and the block / occlusion onsets are not stated; 30 s is read from the figures.
  - IVS V_max, p_0, p_1 are derived from the stated rules (not printed in the tables).
"""
import json
from pathlib import Path

import numpy as np

from ctg_lumped.jongen import model as M

DT, TT = 0.02, 300.0


def scenario(w_um=0.0, p_con=70.0, T_con=60.0, t_con=30.0, put=None, pum=None, ruta=None):
    P = M.J.params(S_WUM=w_um)
    S = M.steady()
    put = M.contraction(TT, DT, t_con=t_con, T_con=T_con, p_con=p_con) if put is None else put
    o, _, ok = M.J.run(P, S, put, dt=DT, pum=pum, ruta=ruta, out_every=5)
    assert ok
    return o


def summary(o):
    i = int(np.argmin(o['fhr'])); f0 = o['fhr'][0]
    return dict(fhr0=float(f0), fhr_min=float(o['fhr'].min()), t_min_min=float(o['t'][i] / 60),
                fhr_max=float(o['fhr'].max()), decel=float(f0 - o['fhr'].min()),
                fhr_pct=float(100 * (o['fhr'].min() / f0 - 1)), po2a_min=float(o['po2_a'].min()),
                po2a_pct=float(100 * (o['po2_a'].min() / o['po2_a'][0] - 1)), ptm_max=float(o['ptm_a'].max()),
                ptm_pct=float(100 * (o['ptm_a'].max() / o['ptm_a'][0] - 1)), q_uta_min=float(60 * o['q_uta'].min()),
                v_ivs_min=float(o['v_ivs'].min()), v_villi_max=float(o['v_villi'].max()))


def ramp(t_on, rise, hold, fall, peak, base=0.0):
    t = np.arange(int(round(TT / DT)) + 1) * DT
    return base + peak * np.interp(t, [0, t_on, t_on + rise, t_on + rise + hold, t_on + rise + hold + fall, TT],
                                   [0, 0, 1, 1, 0, 0])


def main():
    rest = np.full(int(round(TT / DT)) + 1, M.TONE)
    R = {}
    R['A_ref'] = (scenario(), dict(decel=3, po2a_min=17.1, q_uta_min=180, v_ivs_min=121, ptm_max=45.4))
    R['A_110mmHg'] = (scenario(p_con=110), dict(decel=14, fhr_min=121, po2a_min=16.0))
    R['A_120s'] = (scenario(T_con=120, t_con=0), dict(decel=21, fhr_min=114, po2a_min=15.8))
    for hold, pub in [(20, dict(fhr_pct=-11, po2a_pct=-13, ptm_pct=4)), (30, dict(fhr_pct=-24, po2a_pct=-18, ptm_pct=10))]:
        R[f'A_block{hold}s'] = (scenario(put=rest, ruta=ramp(30, 5, hold, 5, 99, 1.0)), pub)
    R['B_ref'] = (scenario(w_um=1), dict(fhr_min=85, fhr_max=145.5, t_min_min=1.33, po2a_min=12.3, ptm_max=61.5, v_villi_max=125))
    for name, kw, pub in [('B_110mmHg', dict(p_con=110), 73), ('B_30mmHg', dict(p_con=30), 132),
                          ('B_120s', dict(T_con=120, t_con=0), 64), ('B_30s', dict(T_con=30, t_con=45), 120)]:
        R[name] = (scenario(w_um=1, **kw), dict(fhr_min=pub))
    for pum, pub in [(25, dict(fhr_pct=-2, po2a_pct=-7, ptm_pct=2)), (33, dict(fhr_pct=-12, po2a_pct=-16, ptm_pct=12)),
                     (41, dict(fhr_pct=-32, po2a_pct=-27, ptm_pct=34)), (90, dict(fhr_pct=-51, po2a_pct=-43, ptm_pct=42))]:
        R[f'B_occl{pum}'] = (scenario(put=rest, pum=ramp(30, 5, 35, 5, pum)), pub)
    out = {}
    for k, (o, pub) in R.items():
        s = summary(o); out[k] = dict(simulated=s, published=pub)
        print(f'{k:11s}', '  '.join(f'{m} {s[m]:.1f} (paper {p})' for m, p in pub.items()))
    Path('outputs').mkdir(exist_ok=True)
    Path('outputs/replication.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
