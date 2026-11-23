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
RHO = 1.05e-3                     # blood density, g/mm3 (with g, mm, s: Pa = g/(mm s2))
E_WALL, H_WALL = 1.5e6, 0.1       # wall elastance, Pa, and wall thickness / radius (Clark 2018 code)
VESSELS = [('uterine', 1, 2.0, 100.0), ('arcuate', 1, 2.0, 18.0), ('radial', 50, 0.2, 6.0),
           ('anastomosis', 50, 0.2, 9.0)]                     # (class, number, radius mm, length mm)
R_L = 1.6                         # Pa s/mm3: Clark's spiral artery + intervillous terminal load (replaced by CK_RSA)
N_SIDES = 2
TO_MMHG_S_ML = 1000.0 / 133.322   # Pa s/mm3 -> mmHg s/ml (also Pa s2/mm3 -> mmHg s2/ml)
TO_ML_MMHG = 133.322 / 1000.0     # mm3/Pa -> ml/mmHg


def resistances():
    """Resistance of each vessel class of one side, Pa s/mm3, plus up (arteries in series), venous and A (bypass)."""
    R = {n: 81.0 * MU * L / (8.0 * math.pi * r ** 4) / k for n, k, r, L in VESSELS}
    R['up'] = R['uterine'] + R['arcuate'] + R['radial']
    R['venous'] = R['up'] / 2.0
    R['A'] = R['anastomosis'] + R['venous']
    return R


def inertances():
    """Inertance L (Pa s2/mm3) and compliance C (mm3/Pa) of each vessel class of one side, from the per-length terms of
    the authors' code (after Mo et al. 1988): L' = 9 rho / (4 pi r^2) per vessel, C' = 3 pi r^3 / (2 h r E) per vessel;
    n vessels in parallel divide L and multiply C."""
    L = {n: 9.0 * RHO * l / (4.0 * math.pi * r ** 2) / k for n, k, r, l in VESSELS}
    C = {n: 3.0 * math.pi * r ** 3 * l / (2.0 * H_WALL * r * E_WALL) * k for n, k, r, l in VESSELS}
    return L, C


def E(name, symbol, value, unit, meaning, detail, cls='paper', note=''):
    return dict(name=name, symbol=symbol, value=float(value), unit=unit, meaning=meaning, part='uterus', equations='',
                **{'class': cls}, source='Clark 2018' if cls == 'paper' else 'ours', source_detail=detail, note=note)


def entries():
    R = resistances(); k = TO_MMHG_S_ML / N_SIDES
    return [
        E('CK_MODE', 'switch', 0, '–', 'uterine arteries: 0 Jongen 2016 (single R_uta), 1 Clark 2018 network (resistances), 2 with inertance and compliance (Ursino mother only)',
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
    ] + rlc_entries()


def rlc_entries():
    """Per-class R, L, C for CK_MODE 2 (network with inertance and compliance, Ursino mother only), both sides."""
    R = resistances(); L, C = inertances(); kr = TO_MMHG_S_ML / N_SIDES; kc = TO_ML_MMHG * N_SIDES
    out = []
    for key, cls in (('U', 'uterine'), ('AR', 'arcuate'), ('R', 'radial'), ('AV', 'anastomosis')):
        out.append(E('CK_R' + key, 'R_' + key, R[cls] * kr, 'mmHg s/ml', cls + ' resistance, both sides',
                     '81 mu L/(8 pi r^4)/n; x 7.5006 / 2'))
        out.append(E('CK_L' + key, 'L_' + key, L[cls] * kr, 'mmHg s2/ml', cls + ' inertance, both sides',
                     '9 rho L/(4 pi r^2)/n (Clark code, after Mo 1988); x 7.5006 / 2'))
        out.append(E('CK_C' + key, 'C_' + key, C[cls] * kc, 'ml/mmHg', cls + ' compliance, both sides',
                     '3 pi r^3 L/(2 h r E) n, E 1.5 MPa, h 0.1 (Clark code); x 0.133322 x 2'))
    out.append(E('CK_RD', 'R_D', R['venous'] * kr, 'mmHg s/ml', 'venous resistance of the bypass, both sides',
                 'R_up / 2 (Clark); CK_RAV + CK_RD = CK_RA'))
    return out


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
