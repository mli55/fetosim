"""Write placental_o2_hill1972.json: parameters of the Hill, Power & Longo 1972 placental O2 exchanger (Am J Physiol
222:721-729, doi:10.1152/ajplegacy.1972.222.3.721), in the index order of model/hill1972.c.

Values are the paper's: human dissociation constants and hemoglobin concentrations (section "Results using human O2
dissociation curves and hemoglobin concentration"), and the Table 1 capillary volumes and membrane diffusing capacity,
which the paper measured in sheep and also used for its human calculation. The paper does not print the fetal human
k'c; the adult human value is used (the paper's human results change by less than 0.1 mmHg for 151-250).
H_DMC is not from the paper: it is the membrane diffusing capacity used when the exchanger runs inside the model,
calibrated to Jongen's resting fetal arterial pO2 (tests/test_hill1972.py checks this).

Run from the repo root: python model/parameters/build_placental_o2_hill1972.py
"""
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent / 'placental_o2_hill1972.json'
HUM, T1 = 'Hill 1972, human dissociation curves (p. 726)', 'Hill 1972 Table 1'


def E(name, symbol, value, unit, meaning, detail, eq='', cls='published', note=''):
    return dict(name=name, symbol=symbol, value=float(value), unit=unit, meaning=meaning, part='placenta', equations=eq,
                **{'class': cls}, source='Hill 1972' if cls == 'published' else 'ours', source_detail=detail, note=note)


def entries():
    return [
        E('H_MODE', 'switch', 0, '–', 'placental O2 transfer: 0 Jongen 2016 A(16), 1 Hill 1972 exchanger', 'off by default',
          cls='our choice'),
        E('H_K1M', 'k1,m', 1.421, '–', 'maternal dissociation curve constant', HUM, 'H(3)'),
        E('H_K2M', 'k2,m', 0.456, '–', 'maternal Bohr constant', HUM, 'H(3)'),
        E('H_K3M', 'k3,m', 0.373, '–', 'maternal dissociation curve slope constant', HUM, 'H(3)'),
        E('H_K1F', 'k1,f', 1.302, '–', 'fetal dissociation curve constant', HUM, 'H(3)'),
        E('H_K2F', 'k2,f', 0.464, '–', 'fetal Bohr constant', HUM, 'H(3)'),
        E('H_K3F', 'k3,f', 0.395, '–', 'fetal dissociation curve slope constant', HUM, 'H(3)'),
        E('H_HBM', 'Hb_m', 12.5, 'g/100 ml', 'maternal hemoglobin concentration', HUM, 'H(4)'),
        E('H_HBF', 'Hb_f', 15.5, 'g/100 ml', 'fetal hemoglobin concentration', HUM, 'H(4)'),
        E('H_PHM', 'pH_m', 7.4, '–', 'maternal blood pH', T1, 'H(3)'),
        E('H_PHF', 'pH_f', 7.35, '–', 'fetal blood pH', T1, 'H(3)'),
        E('H_VM', 'V_m', 10.0, 'ml', 'maternal placental capillary blood volume', T1, 'H(1), H(5)',
          note='sheep value, also used by the paper for humans'),
        E('H_VF', 'V_f', 10.0, 'ml', 'fetal placental capillary blood volume', T1, 'H(2), H(5)',
          note='sheep value, also used by the paper for humans'),
        E('H_DM', 'D_M', 2.73, 'ml/(min mmHg)', 'diffusing capacity of the placental membrane', T1, 'H(5)',
          note='sheep value, also used by the paper for humans'),
        E('H_KCM', "k'c,m", 164.0, '1/(mM s)', 'O2 - red cell reaction rate constant, maternal blood', 'Hill 1972 p. 722 (human blood)', 'H(6)'),
        E('H_KCF', "k'c,f", 164.0, '1/(mM s)', 'O2 - red cell reaction rate constant, fetal blood',
          'adult human value; the fetal human value is not printed in the paper', 'H(6)', cls='our choice'),
        E('H_KCREF', "k'c,human", 164.0, '1/(mM s)', 'reference reaction rate constant of the theta fit', 'Hill 1972 p. 722', 'H(6)'),
        E('H_C1', 'c1', 3.287, 'ml/(ml min mmHg)', 'theta fit constant', 'Hill 1972 Eq. 6', 'H(6)'),
        E('H_C2', 'c2', 0.1117, '–', 'theta fit constant', 'Hill 1972 Eq. 6', 'H(6)'),
        E('H_C3', 'c3', 7.05e-3, 'ml/(ml min mmHg)', 'theta fit constant', 'Hill 1972 Eq. 6', 'H(6)'),
        E('H_C4', 'c4', 0.8142, 'ml/(ml min mmHg)', 'theta fit constant', 'Hill 1972 Eq. 6', 'H(6)'),
        E('H_NSTEP', 'n', 100, '–', 'integration steps along the capillary', 'as the paper (100 steps)', cls='our choice'),
        E('H_QMIN', 'q_min', 0.01, 'ml/s', 'capillary flow below which no exchange is computed', 'numerical guard', cls='our choice'),
        E('H_DMC', 'D_M,model', 2.0, 'ml/(min mmHg)', 'membrane diffusing capacity used inside the model (H_MODE 1)',
          'calibrated: resting fetal arterial pO2 18.5 mmHg with the Jongen 2016 mother', 'H(5)', cls='our choice',
          note='replaces H_DM = 2.73 (19.3 mmHg at rest), because the fetal reflex parameters are tuned to 18.5 mmHg'),
    ]


if __name__ == '__main__':
    doc = {'model': 'Hill, Power & Longo 1972 concurrent placental O2 exchanger (human constants)',
           'order': 'Entries follow the index order of model/hill1972.c; appended after maternal_circulation_ursino1998.json.',
           'parameters': entries()}
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')
    print(len(doc['parameters']), 'parameters ->', OUT.name)
