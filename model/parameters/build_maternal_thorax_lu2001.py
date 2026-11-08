"""Write maternal_thorax_lu2001.json: maternal thorax parameters appended after maternal_hr_randall2019.json.

Structure after Lu et al. 2001 (Am J Physiol Heart 281:H2661; runnable form Physiome JSim model 0154, VS001): pleural
pressure surrounds the heart, the intrathoracic aorta, the intrathoracic vena cava and the pulmonary blood. Jongen's
maternal circulation is lumped (one heart, arteries, tissue, veins), so the graft is lumped too:
  - a fraction MT_FA of the arterial compliance is intrathoracic (absolute p_a,m rises by MT_FA dpth);
  - venous return crosses MT_RVR (no backflow) from the venous pool into an intrathoracic reservoir (vena cava +
    pulmonary blood; compliance MT_CTH, surrounded by dpth);
  - the heart fills from the reservoir's transmural pressure and ejects against the transmural aortic pressure;
  - a fraction MT_FVAB of the venous compliance is abdominal, and abdominal pressure rises by MT_KAB dpth.
MT_V0TH and MT_DPF are derived from Jongen's rest (simulate.set_thorax_rest) so rest is unchanged; the values written
here are for the default parameter set.

MT_FA and MT_FVAB are calibrated on a non-pregnant Valsalva manoeuvre (Lu 2001 Table 13 measured values, Fox et al.);
see the ARPH experiment maternal_thorax_graft_2026-11-07. Run from the repo root:
python model/parameters/build_maternal_thorax_lu2001.py [MT_FA MT_FVAB]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'maternal_thorax_lu2001.json'
sys.path.insert(0, str(HERE.parents[2]))                       # folder holding fetosim/

# Lu 2001 / VS001 values (HighlyIntHuman_wIntervention_0.mod)
K1_VC = 0.15                                                   # vena cava P-V slope, mmHg/ml
C_PULM = 0.8 + 1.5 + 3.75 + 8 + 15                             # Ctpap + Ctpad + Cpa + Cpc + Cpv, ml/mmHg
V_TH = 200.65 + 18.5 + 30.85 + 68.74 + 114.7 + 201.55          # initial Vvc + pulmonary volumes, ml
R_VC = 0.001 * (350 / 200.65) ** 2 + 0.025                     # Rvc = KR (Vmax_vc/Vvc)^2 + Ro at the initial Vvc


def entries(fa, fvab, v0th, dpf):
    E = lambda n, s, v, u, m, c, src, det, note='': dict(
        name=n, symbol=s, value=float(v), unit=u, meaning=m, part='mother', equations='', **{'class': c},
        source=src, source_detail=det, note=note)
    lu = 'Lu 2001 / Physiome VS001 (.mod)'
    return [
        E('MT_MODE', 'mode', 0, '–', 'maternal thorax: 0 off (Jongen 2016), 1 on', 'our choice', 'ours', 'off by default'),
        E('MT_FA', 'f_a', fa, '–', 'intrathoracic share of maternal arterial compliance', 'fitted', 'ours',
          'fitted to the phase-1 SBP rise of a 40 mmHg Valsalva (Fox, measured, in Lu 2001 Table 13)',
          'Lu/VS001 itself: Caop/(Caop+Caod+Csap) = 0.11'),
        E('MT_FVAB', 'f_v,ab', fvab, '–', 'abdominal share of maternal venous compliance', 'fitted', 'ours',
          'fitted to the phase-2 cardiac output of a 40 mmHg Valsalva (Fox, measured, in Lu 2001 Table 13)',
          'Lu 2001 has no abdomen; its stiff venous P-V curve plays this role'),
        E('MT_KAB', 'k_ab', 1.0, '–', 'abdominal pressure rise per thoracic pressure rise', 'our choice', 'ours',
          'strain and pushing are held with the abdominal wall; intra-abdominal and intrathoracic pressures rise together'),
        E('MT_CTH', 'C_th', 1 / K1_VC + C_PULM, 'ml/mmHg', 'intrathoracic venous + pulmonary compliance', 'paper', lu,
          '1/K1 (vena cava) + Ctpap + Ctpad + Cpa + Cpc + Cpv'),
        E('MT_VTH', 'V_th', V_TH, 'ml', 'intrathoracic venous + pulmonary blood volume at rest (moved from the venous pool)',
          'paper', lu, 'initial Vvc + Vpap + Vpad + Vpa + Vpc + Vpv', 'bookkeeping only; does not change dynamics'),
        E('MT_V0TH', 'V_0,th', v0th, 'ml', 'reservoir unstressed volume', 'our choice', 'ours',
          'MT_VTH - MT_CTH (p_v,m - MT_DPF) at rest', 'recomputed by simulate.run'),
        E('MT_DPF', 'Delta p_f', dpf, 'mmHg', 'resting pressure drop across MT_RVR, added back to the filling pressure',
          'our choice', 'ours', 'MT_RVR x resting maternal CO', 'recomputed by simulate.run'),
        E('MT_RVR', 'R_vr', R_VC, 'mmHg·s/ml', 'venous return resistance into the thorax', 'paper', lu,
          'Rvc = KR (Vmax_vc/Vvc)^2 + Ro at the initial Vvc'),
    ]


def write(fa, fvab):
    base = entries(fa, fvab, 0.0, 0.0)
    doc = {'model': 'Maternal thorax after Lu et al. 2001 (lumped graft onto Jongen 2016); appended after '
                    'maternal_hr_randall2019.json',
           'order': 'Entries continue the index order of model/core.c after MR_DS; do not reorder.',
           'parameters': base}
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')
    from fetosim.model import simulate as S
    S.reload_parameters()
    P = S.params(MT_MODE=1); St = S.steady().copy(); S.set_thorax_rest(P, St)
    doc['parameters'] = entries(fa, fvab, P[S.IDX['MT_V0TH']], P[S.IDX['MT_DPF']])
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    a = [float(x) for x in sys.argv[1:3]] if len(sys.argv) > 2 else [0.11, 0.5]
    write(*a)
    print('wrote', OUT)
