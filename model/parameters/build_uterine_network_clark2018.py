"""Write uterine_network_clark2018.json: the uterine arterial network of Clark et al. 2018 (Placenta 66:74-81,
doi:10.1016/j.placenta.2018.05.001) as lumped resistances, in the index order of model/core.c (CK_*).

Geometry and viscosity are the paper's baseline as coded by the authors (github.com/VirtualPregnancy/
network-wave-transmission, VesselDefinition.py). Resistance of n parallel vessels of radius r and length L:
R = 81 mu L / (8 pi r^4) / n (the authors' Poiseuille form). The venous resistance of the anastomotic path is half the
upstream arterial resistance (Clark). Clark models one uterine artery; N_SIDES = 2 puts both sides in parallel.
Not from the paper: CK_RSA (derived at run time, see simulate.clark_parameters) and CK_KA (the paper has no
contraction; calibrated on the Part A scenarios of Jongen 2016, tests/test_clark2018.py).

Run from the repo root: python model/parameters/build_uterine_network_clark2018.py
"""
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent / 'uterine_network_clark2018.json'
MU = 3.4e-3                       # blood viscosity, Pa s
VESSELS = [('uterine', 1, 2.0, 100.0), ('arcuate', 1, 2.0, 18.0), ('radial', 50, 0.2, 6.0),
           ('anastomosis', 50, 0.2, 9.0)]                     # (class, number, radius mm, length mm)
R_L = 1.6                         # Pa s/mm3: Clark's spiral artery + intervillous terminal load (replaced by CK_RSA)
N_SIDES = 2
TO_MMHG_S_ML = 1000.0 / 133.322   # Pa s/mm3 -> mmHg s/ml


def resistances():
    """Resistance of each vessel class of one side, Pa s/mm3, plus up (arteries in series), venous and A (bypass)."""
    R = {n: 81.0 * MU * L / (8.0 * math.pi * r ** 4) / k for n, k, r, L in VESSELS}
    R['up'] = R['uterine'] + R['arcuate'] + R['radial']
    R['venous'] = R['up'] / 2.0
    R['A'] = R['anastomosis'] + R['venous']
    return R


def E(name, symbol, value, unit, meaning, detail, cls='paper', note=''):
    return dict(name=name, symbol=symbol, value=float(value), unit=unit, meaning=meaning, part='uterus', equations='',
                **{'class': cls}, source='Clark 2018' if cls == 'paper' else 'ours', source_detail=detail, note=note)


def entries():
    R = resistances(); k = TO_MMHG_S_ML / N_SIDES
    return [
        E('CK_MODE', 'switch', 0, '–', 'uterine arteries: 0 Jongen 2016 (single R_uta), 1 Clark 2018 network',
          'off by default', cls='our choice'),
        E('CK_RUP', 'R_up', R['up'] * k, 'mmHg s/ml', 'uterine + arcuate + radial arteries in series, both sides',
          'baseline geometry, 81 mu L/(8 pi r^4)/n; x 7.5006 (unit) / 2 (sides)', note='radial arteries are 91% of it'),
        E('CK_RSA', 'R_sa', 1.0, 'mmHg s/ml', 'spiral arteries, from the network node to the intervillous space',
          'placeholder; simulate.clark_parameters sets it so that the resting intervillous inflow and pressure are Jongen\'s',
          cls='our choice', note='replaces Clark\'s lumped spiral artery + intervillous load R_L = 1.6 Pa s/mm3'),
        E('CK_RA', 'R_A', R['A'] * k, 'mmHg s/ml', 'arterio-venous anastomoses and their veins (bypass), both sides',
          'anastomoses 50 x (r 0.2 mm, L 9 mm) + venous resistance = R_up / 2'),
        E('CK_KA', 'k_A', 0.7, '–', 'share of uterine pressure that acts as outflow pressure of the bypass',
          'Starling resistor against max(p_v,m, k_A p_ut); 0 = not compressed. Calibrated on the three Part A '
          'contraction scenarios of Jongen 2016 (grid step 0.05)', cls='our choice',
          note='Clark 2018 has no contraction'),
    ]


def main():
    R = resistances()
    doc = dict(
        model='Uterine arterial network after Clark et al. 2018 (Placenta 66:74-81); appended after placental_o2_hill1972.json',
        order='Entries continue the index order of model/core.c after H_DMC; do not reorder.',
        clark_one_side_Pa_s_mm3={n: R[n] for n in ('uterine', 'arcuate', 'radial', 'anastomosis', 'venous', 'up', 'A')},
        clark_terminal_load_Pa_s_mm3=R_L, n_sides=N_SIDES, parameters=entries())
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')
    print({e['name']: e['value'] for e in doc['parameters']})


if __name__ == '__main__':
    main()
