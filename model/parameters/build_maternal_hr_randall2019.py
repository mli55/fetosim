"""Write maternal_hr_randall2019.json: the maternal heart-rate reflex parameters appended after jongen2016.json.

Randall et al. 2019 (J Appl Physiol 127:1386; code Zenodo 10.5281/zenodo.3234934, load_global.m) build their nominal
parameters from a subject's resting values. We apply the same rules to the resting values below. Rest is then a steady
state with MHR = 60 / M_T, so the coupled model starts from Jongen's steady state unchanged.

SBP input: Jongen's maternal heart is not pulsatile, so p_a,m is a mean pressure. We use SBP = p_a,m + MR_SBP_OFF
(a constant pulse offset, so changes in SBP equal changes in p_a,m). MR_SBP_OFF = SBP_REST - p_a,m at Jongen's rest.
Because Randall's wall-strain sigmoid is centred on the resting SBP (s_w = SBP_REST), the value of SBP_REST does not
change any dynamics; only Delta p_a,m does.

Run from the repo root: python model/parameters/build_maternal_hr_randall2019.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'maternal_hr_randall2019.json'
sys.path.insert(0, str(HERE.parents[2]))                       # folder holding fetosim/

# resting values of a labouring woman (our choices; see notes)
REST = dict(SBP_REST=115.0, PTH_REST=4.5, AGE=30.0, RSA_RANGE=10.0)


def entries(sbp_off):
    from fetosim.model import randall2019 as R
    hbar = 60.0 / 0.75                                          # Jongen 2016 A-T1b, M_T = 0.75 s
    p = R.nominal_parameters(Pbar=REST['SBP_REST'], Pthbar=REST['PTH_REST'], Hbar=hbar,
                             HminR=hbar - REST['RSA_RANGE'] / 2, HmaxR=hbar + REST['RSA_RANGE'] / 2, age=REST['AGE'])
    meta = {
        'A': ('A', '–', 'cross-sectional area ratio (wall strain)', 'paper'),
        'B': ('B', '–', 'carotid share of the neural integration (1 - B is aortic)', 'paper'),
        'Kb': ('K_b', '–', 'afferent baroreceptor gain', 'paper'),
        'Kpb': ('K_pb', '–', 'baroreflex parasympathetic gain', 'paper'),
        'Kpr': ('K_pr', '–', 'respiratory parasympathetic gain', 'paper'),
        'Ks': ('K_s', '–', 'baroreflex sympathetic gain', 'paper'),
        'taub': ('tau_b', 's', 'baroreceptor strain time constant', 'paper'),
        'taupb': ('tau_pb', 's', 'baroreflex parasympathetic time constant', 'paper'),
        'taupr': ('tau_pr', 's', 'respiratory parasympathetic time constant', 'paper'),
        'taus': ('tau_s', 's', 'baroreflex sympathetic time constant', 'paper'),
        'tauH': ('tau_H', 's', 'heart-rate time constant', 'paper'),
        'qw': ('q_w', '1/mmHg', 'wall-strain sigmoid steepness', 'paper'),
        'qpb': ('q_pb', 's', 'parasympathetic sigmoid steepness', 'paper'),
        'qpr': ('q_pr', '1/mmHg', 'respiratory sigmoid steepness', 'paper'),
        'qs': ('q_s', 's', 'sympathetic sigmoid steepness', 'paper'),
        'sw': ('s_w', 'mmHg', 'wall-strain half saturation = resting SBP', 'our choice'),
        'spb': ('s_pb', '1/s', 'parasympathetic half saturation (80 % parasympathetic at rest)', 'our choice'),
        'spr': ('s_pr', 'mmHg', 'respiratory half saturation = resting Pth', 'our choice'),
        'ss': ('s_s', '1/s', 'sympathetic half saturation (20 % sympathetic at rest)', 'our choice'),
        'HI': ('H_I', 'bpm', 'intrinsic heart rate, 118 - 0.57 age', 'our choice'),
        'Hpb': ('H_pb', '–', 'parasympathetic HR gain (rest MHR = 60 / M_T)', 'our choice'),
        'Hpr': ('H_pr', '–', 'respiratory HR gain (resting RSA range)', 'our choice'),
        'Hs': ('H_s', '–', 'sympathetic HR gain, from maximal HR 208 - 0.7 age', 'our choice'),
        'Ds': ('D_s', 's', 'sympathetic delay', 'paper'),
    }
    rows = [
        dict(name='MR_MODE', symbol='mode', value=0.0, unit='–',
             meaning='maternal HR reflex: 0 off (Jongen 2016), 1 open loop, 2 closed loop (maternal period = 60/MHR)',
             part='mother', equations='', **{'class': 'our choice'}, source='ours', source_detail='off by default',
             note=''),
        dict(name='MR_SBP_OFF', symbol='Delta_SBP', value=sbp_off, unit='mmHg',
             meaning='SBP input = p_a,m + offset (Jongen maternal heart is not pulsatile)', part='mother',
             equations='', **{'class': 'our choice'}, source='ours',
             source_detail=f"SBP_REST {REST['SBP_REST']} - p_a,m at Jongen rest",
             note='recomputed by simulate.run when other parameters are changed'),
        dict(name='MR_PTH_REST', symbol='P_th,rest', value=REST['PTH_REST'], unit='mmHg',
             meaning='resting maternal thoracic pressure (no breathing input)', part='mother', equations='',
             **{'class': 'our choice'}, source='ours',
             source_detail='Randall subject rest 4.49 mmHg; their breathing signal spans 3.5-6 mmHg', note=''),
    ]
    for name, v in zip(R.NAMES, p):
        sym, unit, meaning, cls = meta[name]
        rows.append(dict(name='MR_' + name.upper(), symbol=sym, value=float(v), unit=unit, meaning=meaning,
                         part='mother', equations='Randall 2019 Eqs. 6-16', **{'class': cls},
                         source='Randall 2019 Table 3 / load_global.m',
                         source_detail='nominal value' if cls == 'paper' else
                         f"load_global.m rule with SBP {REST['SBP_REST']}, Pth {REST['PTH_REST']}, MHR {hbar:g}, "
                         f"RSA range {REST['RSA_RANGE']} bpm, age {REST['AGE']:g}",
                         note=''))
    return rows


def write(sbp_off):
    OUT.write_text(json.dumps({
        'model': 'Maternal heart-rate reflex: Randall et al. 2019 (J Appl Physiol 127:1386), nominal parameters for '
                 'a labouring woman at rest; appended after jongen2016.json',
        'order': 'Entries continue the index order of model/core.c after S_MEV; do not reorder.',
        'rest': REST,
        'parameters': entries(sbp_off)}, indent=1) + '\n')


if __name__ == '__main__':
    if not OUT.exists():
        write(0.0)                                              # bootstrap; the offset does not affect rest at mode 0
    from fetosim.model import simulate as S
    s0 = S.steady_state(S.params(MR_MODE=0))
    write(REST['SBP_REST'] - S.maternal_pa(s0, S.params()))
    print(OUT, 'p_a,m rest', S.maternal_pa(s0, S.params()))
