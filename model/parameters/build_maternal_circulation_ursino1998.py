"""Write maternal_circulation_ursino1998.json: every parameter of Ursino 1998 (Am J Physiol 275:H1733-H1747,
doi:10.1152/ajpheart.1998.275.5.H1733), in the index order of model/ursino1998.c.

All values are copied from the paper's Tables 1-3 and their footnotes. They describe a non-pregnant 70-kg adult
(total blood volume 5,300 ml, heart rate 72 beats/min). The seven switches at the end are not in the paper; they
reproduce its open-loop experiments (isolated carotid sinus, pacing, vagotomy, single-effector knock-out) and are 0
for the published closed-loop model.

Run from the repo root: python model/parameters/build_maternal_circulation_ursino1998.py
"""
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent / 'maternal_circulation_ursino1998.json'
T1, T2, T3 = 'Ursino 1998 Table 1', 'Ursino 1998 Table 2', 'Ursino 1998 Table 3'


def E(name, symbol, value, unit, meaning, detail, eq='', cls='published', note=''):
    return dict(name=name, symbol=symbol, value=float(value), unit=unit, meaning=meaning, part='mother', equations=eq,
                **{'class': cls}, source='Ursino 1998' if cls == 'published' else 'ours', source_detail=detail, note=note)


def entries():
    C, V, R = 'ml/mmHg', 'ml', 'mmHg s/ml'
    e = [
        E('U_CSA', 'C_sa', 0.28, C, 'systemic arterial compliance', T1, 'U(5)'),
        E('U_CSP', 'C_sp', 2.05, C, 'splanchnic peripheral compliance', T1, 'U(7)'),
        E('U_CEP', 'C_ep', 1.67, C, 'extrasplanchnic peripheral compliance', T1, 'U(7)'),
        E('U_CSV', 'C_sv', 61.11, C, 'splanchnic venous compliance', T1, 'U(9)'),
        E('U_CEV', 'C_ev', 50.0, C, 'extrasplanchnic venous compliance', T1, 'U(8)'),
        E('U_CPA', 'C_pa', 0.76, C, 'pulmonary arterial compliance', T1, 'U(1)'),
        E('U_CPP', 'C_pp', 5.80, C, 'pulmonary peripheral compliance', T1, 'U(3)'),
        E('U_CPV', 'C_pv', 25.37, C, 'pulmonary venous compliance', T1, 'U(4)'),
        E('U_VUSA', 'V_u,sa', 0.0, V, 'systemic arterial unstressed volume', T1, 'U(10)'),
        E('U_VUSP', 'V_u,sp', 274.4, V, 'splanchnic peripheral unstressed volume', T1, 'U(10)'),
        E('U_VUEP', 'V_u,ep', 336.6, V, 'extrasplanchnic peripheral unstressed volume', T1, 'U(10)'),
        E('U_VUPA', 'V_u,pa', 0.0, V, 'pulmonary arterial unstressed volume', T1, 'U(10)'),
        E('U_VUPP', 'V_u,pp', 123.0, V, 'pulmonary peripheral unstressed volume', T1, 'U(10)'),
        E('U_VUPV', 'V_u,pv', 120.0, V, 'pulmonary venous unstressed volume', T1, 'U(10)'),
        E('U_RSA', 'R_sa', 0.06, R, 'systemic arterial resistance', T1, 'U(6)'),
        E('U_RSV', 'R_sv', 0.038, R, 'splanchnic venous resistance', T1, 'U(23)'),
        E('U_REV', 'R_ev', 0.016, R, 'extrasplanchnic venous resistance', T1, 'U(8)'),
        E('U_RPA', 'R_pa', 0.023, R, 'pulmonary arterial resistance', T1, 'U(2)'),
        E('U_RPP', 'R_pp', 0.0894, R, 'pulmonary peripheral resistance', T1, 'U(3)'),
        E('U_RPV', 'R_pv', 0.0056, R, 'pulmonary venous resistance', T1, 'U(4)'),
        E('U_LSA', 'L_sa', 0.22e-3, 'mmHg s2/ml', 'systemic arterial inertance', T1, 'U(6)',
          note='Table 1 prints the unit as mmHg ml s-2; dimensionally it is mmHg s2/ml'),
        E('U_LPA', 'L_pa', 0.18e-3, 'mmHg s2/ml', 'pulmonary arterial inertance', T1, 'U(2)'),
        E('U_VT', 'V_t', 5300.0, V, 'total blood volume', T1 + ' footnote', 'U(9), U(11)'),
        E('U_CLA', 'C_la', 19.23, C, 'left atrial compliance', T2, 'U(12)'),
        E('U_VULA', 'V_u,la', 25.0, V, 'left atrial unstressed volume', T2, 'U(10)'),
        E('U_RLA', 'R_la', 2.5e-3, R, 'mitral valve resistance', T2, 'U(13)'),
        E('U_P0LV', 'P_0,lv', 1.5, 'mmHg', 'left ventricle end-diastolic pressure-volume scale', T2, 'U(18)'),
        E('U_KELV', 'k_E,lv', 0.014, '1/ml', 'left ventricle end-diastolic pressure-volume exponent', T2, 'U(18)'),
        E('U_VULV', 'V_u,lv', 16.77, V, 'left ventricle unstressed volume', T2, 'U(18)'),
        E('U_KRLV', 'k_R,lv', 3.75e-4, 's/ml', 'left ventricle viscous resistance per unit isometric pressure', T2, 'U(16)'),
        E('U_CRA', 'C_ra', 31.25, C, 'right atrial compliance', T2, 'U(23)'),
        E('U_VURA', 'V_u,ra', 25.0, V, 'right atrial unstressed volume', T2, 'U(10)'),
        E('U_RRA', 'R_ra', 2.5e-3, R, 'tricuspid valve resistance', T2, 'U(24)'),
        E('U_P0RV', 'P_0,rv', 1.5, 'mmHg', 'right ventricle end-diastolic pressure-volume scale', T2, 'U(29)'),
        E('U_KERV', 'k_E,rv', 0.011, '1/ml', 'right ventricle end-diastolic pressure-volume exponent', T2, 'U(29)'),
        E('U_VURV', 'V_u,rv', 40.8, V, 'right ventricle unstressed volume', T2, 'U(29)'),
        E('U_KRRV', 'k_R,rv', 1.4e-3, 's/ml', 'right ventricle viscous resistance per unit isometric pressure', T2, 'U(27)'),
        E('U_TSYS0', 'T_sys,0', 0.5, 's', 'systole duration at infinite heart period', T2 + ' footnote', 'U(22)'),
        E('U_KSYS', 'k_sys', 0.075, 's2', 'shortening of systole with heart rate', T2 + ' footnote', 'U(22)'),
        E('U_PN', 'P_n', 92.0, 'mmHg', 'carotid sinus pressure at the centre of the baroreceptor sigmoid', T3, 'U(31)'),
        E('U_FMIN', 'f_min', 2.52, 'spikes/s', 'lower saturation of sinus nerve firing', T3, 'U(31)'),
        E('U_FMAX', 'f_max', 47.78, 'spikes/s', 'upper saturation of sinus nerve firing', T3, 'U(31)'),
        E('U_KA', 'k_a', 11.758, 'mmHg', 'slope parameter of the baroreceptor sigmoid', T3, 'U(31), U(32)'),
        E('U_TAUZ', 'tau_z', 6.37, 's', 'baroreceptor zero (rate sensitivity) time constant', T3, 'U(30)'),
        E('U_TAUP', 'tau_p', 2.076, 's', 'baroreceptor pole time constant', T3, 'U(30)'),
        E('U_FESINF', 'f_es,inf', 2.10, 'spikes/s', 'sympathetic firing at high sinus nerve firing', T3, 'U(33)'),
        E('U_FES0', 'f_es,0', 16.11, 'spikes/s', 'sympathetic firing at zero sinus nerve firing', T3, 'U(33)'),
        E('U_KES', 'k_es', 0.0675, 's', 'sympathetic firing decay with sinus nerve firing', T3, 'U(33)'),
        E('U_FESMIN', 'f_es,min', 2.66, 'spikes/s', 'sympathetic firing threshold of the effectors', T3, 'U(35)'),
        E('U_FEV0', 'f_ev,0', 3.2, 'spikes/s', 'vagal firing at low sinus nerve firing', T3, 'U(34)'),
        E('U_FEVINF', 'f_ev,inf', 6.3, 'spikes/s', 'vagal firing at high sinus nerve firing', T3, 'U(34)'),
        E('U_KEV', 'k_ev', 7.06, 'spikes/s', 'slope parameter of the vagal sigmoid', T3, 'U(34)'),
        E('U_FCS0', 'f_cs,0', 25.0, 'spikes/s', 'sinus nerve firing at the centre of the vagal sigmoid', T3, 'U(34)'),
    ]
    eff = [('ELV', 'E_max,lv', 0.475, 8, 2, 2.392, 'mmHg/ml', 'left ventricle end-systolic elastance'),
           ('ERV', 'E_max,rv', 0.282, 8, 2, 1.412, 'mmHg/ml', 'right ventricle end-systolic elastance'),
           ('RSP', 'R_sp', 0.695, 6, 2, 2.49, R, 'splanchnic peripheral resistance'),
           ('REP', 'R_ep', 0.53, 6, 2, 0.78, R, 'extrasplanchnic peripheral resistance'),
           ('VSV', 'V_u,sv', -265.4, 20, 5, 1435.4, V, 'splanchnic venous unstressed volume'),
           ('VEV', 'V_u,ev', -132.5, 20, 5, 1537.0, V, 'extrasplanchnic venous unstressed volume')]
    for k, sym, g, tau, d, th0, unit, what in eff:
        e += [E(f'U_G{k}', f'G_{sym}', g, f'{unit} per nu', f'sympathetic gain on {what}', T3, 'U(35)'),
              E(f'U_TAU{k}', f'tau_{sym}', tau, 's', f'time constant of {what} control', T3, 'U(36)'),
              E(f'U_D{k}', f'D_{sym}', d, 's', f'pure delay of {what} control', T3, 'U(35)'),
              E(f'U_{"VU" if k[0] == "V" else ""}{k[1:] if k[0] == "V" else k}0', f'{sym},0', th0, unit,
                f'{what} without sympathetic drive', T3, 'U(37)')]
    e += [
        E('U_GTS', 'G_T,s', -0.13, 's per nu', 'sympathetic gain on heart period', T3, 'U(38)'),
        E('U_TAUTS', 'tau_T,s', 2.0, 's', 'time constant of the sympathetic heart period control', T3, 'U(39)'),
        E('U_DTS', 'D_T,s', 2.0, 's', 'pure delay of the sympathetic heart period control', T3, 'U(38)'),
        E('U_GTV', 'G_T,v', 0.09, 's per nu', 'vagal gain on heart period', T3, 'U(40)'),
        E('U_TAUTV', 'tau_T,v', 1.5, 's', 'time constant of the vagal heart period control', T3, 'U(41)'),
        E('U_DTV', 'D_T,v', 0.2, 's', 'pure delay of the vagal heart period control', T3, 'U(40)'),
        E('U_T0', 'T_0', 0.58, 's', 'heart period without nervous drive', T3, 'U(42)'),
    ]
    sw = [('U_PCS_MODE', '0: carotid sinus pressure = systemic arterial pressure; 1: given as input (isolated sinus)'),
          ('U_PACE_T', '> 0: heart period fixed at this value in s (artificial pacing)'),
          ('U_VAGOTOMY', '1: vagal firing f_ev = 0'),
          ('U_HOLD_E', '1: end-systolic elastances frozen'), ('U_HOLD_R', '1: peripheral resistances frozen'),
          ('U_HOLD_V', '1: venous unstressed volumes frozen'), ('U_HOLD_T', '1: heart period frozen')]
    e += [E(n, 'switch', 0, '–', m, 'experiment switch; 0 = published closed-loop model', cls='our choice') for n, m in sw]
    return e


if __name__ == '__main__':
    doc = {'model': 'Ursino 1998 closed-loop pulsatile circulation with carotid baroreflex (non-pregnant adult)',
           'order': 'Entries follow the index order of model/ursino1998.c; do not reorder.',
           'parameters': entries()}
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')
    print(len(doc['parameters']), 'parameters ->', OUT.name)
