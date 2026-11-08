"""Write maternal_vasomotor.json: maternal vasomotor reflex (MV_*) and contraction autotransfusion (MA_*) parameters,
appended after maternal_thorax_lu2001.json.

Vasomotor reflex (MV_MODE 1, needs MR_MODE > 0). Jongen's mother has fixed resistances and the Randall 2019 model
controls only heart rate, so a strain could not build up arterial pressure (no Valsalva phase-4 overshoot) and stage-2
contractions lowered rather than raised maternal pressure (Kuhn 2017: SAP +36 %, SVR +88 %). The Randall sympathetic
tone Ts, which already drives HR, now also drives
  - the tissue arterial resistance: R_tis,a = M_RTISA max(MV_RMIN, 1 + MV_KR x_R);
  - the venous unstressed volume: MV_KV x_V ml move to stressed volume (venoconstriction);
x_R and x_V are Ts - MV_TS0 low-pass filtered with the effector time constants MV_TAUR and MV_TAUV (Ursino 1998).
This follows the efferent structure of Lu 2001 and Ursino 1998 (sympathetic tone -> systemic resistance and venous
unstressed volume), with Randall's Ts (its delay Ds and time constant tau_s) as the shared efferent. MV_KR and MV_KV are
fitted to a non-pregnant Valsalva (Fox, measured, Lu 2001 Table 13); the uterine bed is not constricted.

Contraction autotransfusion (MA_MODE 1). The contracting myometrium squeezes blood into the maternal veins:
v_sq = MA_VSQ x_c / MA_PREF ml, with x_c the uterine pressure above rest not due to abdominal (push) pressure. MA_VSQ
is fitted to the cardiac-output rise during contractions in Lee et al. 1989 (+11 %, epidural, left lateral).

Calibration: ARPH experiment maternal_vasomotor_autotransfusion_2026-11-07. Run from the repo root:
python model/parameters/build_maternal_vasomotor.py [MV_KR MV_KV MA_VSQ]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'maternal_vasomotor.json'


def entries(kr, kv, vsq):
    E = lambda n, s, v, u, m, c, src, det, note='': dict(
        name=n, symbol=s, value=float(v), unit=u, meaning=m, part='mother', equations='', **{'class': c},
        source=src, source_detail=det, note=note)
    return [
        E('MV_MODE', 'mode', 0, '–', 'maternal vasomotor reflex: 0 off, 1 on (needs MR_MODE > 0)', 'our choice', 'ours',
          'off by default'),
        E('MV_TS0', 'T_s,0', 0.2, '–', 'resting Randall sympathetic tone', 'our choice', 'ours',
          'Randall 2019 resting Ts (20 % sympathetic)', 'recomputed by simulate.run from the resting state'),
        E('MV_KR', 'k_R', kr, '1/Ts', 'tissue arterial resistance gain per unit Ts above rest', 'fitted', 'ours',
          'fitted to a 40 mmHg Valsalva (Fox, measured, in Lu 2001 Table 13)',
          'structure as Lu 2001 / Ursino 1998 (sympathetic -> systemic resistance)'),
        E('MV_KV', 'k_V', kv, 'ml/Ts', 'venous unstressed volume moved to stressed volume per unit Ts above rest',
          'fitted', 'ours', 'fitted to a 40 mmHg Valsalva (Fox, measured, in Lu 2001 Table 13)',
          'structure as Lu 2001 / Ursino 1998 (sympathetic -> venous unstressed volume)'),
        E('MV_RMIN', 'r_min', 0.3, '–', 'lower bound of the resistance factor', 'our choice', 'ours',
          'keeps R_tis,a positive when Ts falls'),
        E('MV_TAUR', 'tau_R', 6.0, 's', 'time constant of the resistance effector', 'paper', 'Ursino 1998 (Am J Physiol 275:H1733)',
          'tau_R,sp (systemic peripheral resistance)'),
        E('MV_TAUV', 'tau_V', 20.0, 's', 'time constant of the venous unstressed-volume effector', 'paper',
          'Ursino 1998 (Am J Physiol 275:H1733)', 'tau_V,u,sv (systemic venous unstressed volume)'),
        E('MA_MODE', 'mode', 0, '–', 'contraction autotransfusion: 0 off, 1 on', 'our choice', 'ours', 'off by default'),
        E('MA_VSQ', 'V_sq', vsq, 'ml', 'myometrial blood squeezed into the maternal veins at MA_PREF', 'fitted', 'ours',
          'fitted to the CO rise during contractions in Lee et al. 1989 (+11 %)',
          'compare the commonly quoted 300-500 ml per contraction'),
        E('MA_PREF', 'p_ref', 70.0, 'mmHg', 'contraction amplitude at which MA_VSQ is squeezed', 'our choice', 'ours',
          'Jongen reference contraction amplitude above rest (A 2.4)'),
    ]


def write(kr, kv, vsq):
    doc = {'model': 'Maternal vasomotor reflex (Randall Ts -> resistance, venous unstressed volume) and contraction '
                    'autotransfusion; appended after maternal_thorax_lu2001.json',
           'order': 'Entries continue the index order of model/core.c after MT_RVR; do not reorder.',
           'parameters': entries(kr, kv, vsq)}
    OUT.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + '\n')


if __name__ == '__main__':
    a = [float(x) for x in sys.argv[1:4]] if len(sys.argv) > 3 else [0.0, 0.0, 0.0]
    write(*a)
    print('wrote', OUT)
