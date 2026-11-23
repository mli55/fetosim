"""Run the three model versions of the block-by-block section of docs/index.html on one synthetic TOCO and save what
it plots.

Input (input.json): a synthetic 20-minute TOCO at 1 s, its 2-s bins, the four detected peak times and the uterine
pressure made from them (a = 1). Runs, all with w_um = 1, g = 1, du = 1 s, recorded every 0.02 s step:
  cm0  Jongen 2016 (maternal heart period fixed);
  cm1  + maternal layers: Randall 2019 heart-rate reflex closed loop, vasomotor reflex, autotransfusion
       (mhr=2, vaso=1, autotx=1);
  cm2  Ursino 1998 mother and the Clark 2018 uterine network with inertance and compliance (circ=1, uterus=2).
The model is run from a private copy of model/ whose core.c records a few more internal signals (Randall wall strain,
afferent signal and vagal tones; Ursino afferent, efferent and effector signals). The copy adds output lines only;
the dynamics are those of model/core.c.
Saved per run: 1-s samples (every 50th step); for the pulsatile maternal signals of cm2, the mean, max and min over one
heart period centred on each sample (beat mean, systolic, diastolic); and 2-s bin means of FHR and MHR (of the 0.1-s samples).
Run from anywhere: python docs/blocks/run.py, then python docs/blocks/build.py"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

EXTRA = ['x_ewc', 'x_nb', 'x_tpb', 'x_tpr', 'x_vsq', 'x_vvaso',
         'u_fcs', 'u_fes', 'u_fev', 'u_rsp', 'u_rep', 'u_elv', 'u_erv', 'u_vusv', 'u_vuev', 'u_ts', 'u_tv']


def instrumented_model(dst):
    """Copy model/ to dst/fetosim/model and add the EXTRA outputs to core.c and simulate.OUT."""
    m = dst / 'fetosim' / 'model'
    shutil.copytree(REPO / 'model', m, ignore=shutil.ignore_patterns('build', '__pycache__'))
    (dst / 'fetosim' / '__init__.py').touch()
    c = (m / 'core.c').read_text()
    n0 = int(c.split('#define NOUT ')[1].split()[0])
    sub = [('#define NOUT %d' % n0, '#define NOUT %d' % (n0 + len(EXTRA))),
           ('    double sbp = pa_m + P[MR_SBP_OFF], d_r[6] = {0, 0, 0, 0, 0, 0};',
            '    double sbp = pa_m + P[MR_SBP_OFF], d_r[6] = {0, 0, 0, 0, 0, 0}, x_ewc = 0, x_nb = 0;'),
           ('      double nb_r = P[MR_B]*(ewc - S[R_EBC]) + (1 - P[MR_B])*(ewa - S[R_EBA]);',
            '      double nb_r = P[MR_B]*(ewc - S[R_EBC]) + (1 - P[MR_B])*(ewa - S[R_EBA]); x_ewc = ewc; x_nb = nb_r;'),
           ('      o[38] = q_ins; o[39] = p_ins;',
            '      o[38] = q_ins; o[39] = p_ins;\n'
            '      { double *x = o + %d; for (int j = 0; j < %d; j++) x[j] = 0.0;\n' % (n0, len(EXTRA)) +
            '        x[0] = x_ewc; x[1] = x_nb; x[2] = mr ? S[R_TPB] : 0.0; x[3] = mr ? S[R_TPR] : 0.0;\n'
            '        x[4] = ma ? P[MA_VSQ]*x_c/P[MA_PREF] : 0.0; x[5] = mv ? P[MV_KV]*S[MV_XV] : 0.0;\n'
            '        if (mu) { x[6] = ua.fcs; x[7] = ua.fes; x[8] = ua.fev; x[9] = ua.rsp; x[10] = ua.rep; x[11] = ua.elv;\n'
            '          x[12] = ua.erv; x[13] = ua.vusv; x[14] = ua.vuev; x[15] = SU[X_DTS]; x[16] = SU[X_DTV]; } }')]
    for a, b in sub:
        assert c.count(a) == 1, a
        c = c.replace(a, b)
    (m / 'core.c').write_text(c)
    s = (m / 'simulate.py').read_text()
    a = "'q_ins', 'p_ins']"
    assert s.count(a) == 1
    (m / 'simulate.py').write_text(s.replace(a, "'q_ins', 'p_ins', " + ', '.join(repr(k) for k in EXTRA) + ']'))


def main():
    d = json.loads((HERE / 'input.json').read_text())
    put = np.array(d['p_ut']); n = len(put)
    tmp = Path(tempfile.mkdtemp()); instrumented_model(tmp)
    sys.path.insert(0, str(tmp)); from fetosim.model import simulate as S   # noqa: E402
    r = lambda a, nd=2: [round(float(v), nd) for v in a]
    inst = lambda x: x[np.arange(n) * 50]
    def win(x, f=np.mean, per=None):
        """f over a window around each 1-s sample: one heart period per[i] (s) long, so a pulsatile signal gives its
        beat mean, max and min; 1 s when per is None."""
        out = np.empty(n)
        for i in range(n):
            w = int(round((1.0 if per is None else per[i]) / 0.02)); a = min(max(0, i * 50 - w // 2), len(x) - w)
            out[i] = f(x[a:a + w])
        return out

    def common(o, P, pulsatile):
        g = lambda k: P[S.IDX[k]]
        per = 60.0 / inst(o['mhr'])
        pick = lambda k: win(o[k], per=per) if pulsatile else inst(o[k])
        pvl = (o['v_villi'] - g('F_VVIL0')) / g('F_CVIL'); pa = o['ptm_a']                       # A(1), villi
        pum = g('S_WUM') * (o['p_ut'] - g('S_PREST'))                                           # B(2)
        f = 0.5 + np.arctan((0.5 * (pvl + pa) - pum - g('U_P0UMA')) / g('U_P1UMA')) / np.pi     # B(1), A(6)
        x = {k: r(inst(o[k]), nd) for k, nd in [
            ('fhr', 2), ('mhr', 2), ('po2_a', 2), ('ptm_a', 2), ('r_c', 4), ('r_b', 4), ('e_v', 4), ('e_s', 4),
            ('z_tv', 4), ('z_ts', 4), ('spo2_f', 2), ('sto2_f_tis', 2), ('sto2_f_cer', 2), ('spo2_m', 2), ('sto2_m', 2)]}
        for k in ('q_uta', 'q_utv'): x[k] = r(pick(k) * 60, 1)                                 # ml/min
        for k in ('q_uma', 'q_cera'): x[k] = r(inst(o[k]) * 60, 1)
        x['pa_m'] = r(pick('pa_m')); x['pv_m'] = r(pick('pv_m'))
        x['p_vil'] = r(inst(pvl)); x['r_uma'] = r(inst(np.maximum(g('U_R0UMA') / f ** 2, g('F_RUMA'))), 3)
        x['T'] = r(60 / inst(o['fhr']), 4); x['T_m'] = r(60 / inst(o['mhr']), 4)
        x['fhr_bins'] = r(o['fhr'][:(n - 1) * 50:5].reshape(-1, 20).mean(1))    # 0.1-s samples, as model/simulate.simulate
        x['mhr_bins'] = r(o['mhr'][:(n - 1) * 50:5].reshape(-1, 20).mean(1))
        return x

    res = {}
    o, _, ok = S.run(put, du=1.0, out_every=1, w_um=1.0, g=1.0); assert ok
    res['cm0'] = common(o, S.params(S_WUM=1.0, S_GAIN=1.0), False)

    o, _, ok = S.run(put, du=1.0, out_every=1, w_um=1.0, g=1.0, mhr=2, vaso=1, autotx=1); assert ok
    x = common(o, S.params(S_WUM=1.0, S_GAIN=1.0), False)
    x['rtisa_m'] = r(inst(o['rtisa_m']), 4); x['ts_m'] = r(inst(o['ts_m']), 4)
    for k, kk in (('ewc', 'x_ewc'), ('nb', 'x_nb'), ('tpb', 'x_tpb'), ('tpr', 'x_tpr')): x[k] = r(inst(o[kk]), 4)
    x['v_sq'] = r(inst(o['x_vsq']), 1); x['v_vaso'] = r(inst(o['x_vvaso']), 1)
    res['cm1'] = x

    over, _ = S.steady_clark_ursino(0, rlc=1)
    o, _, ok = S.run(put, du=1.0, out_every=1, w_um=1.0, g=1.0, circ=1, uterus=2); assert ok
    x = common(o, S.params(S_WUM=1.0, S_GAIN=1.0, MU_MODE=1, **over), True)
    per = 60.0 / inst(o['mhr'])
    for k in ('q_ut_art', 'q_anast'): x[k] = r(win(o[k], per=per) * 60, 1)
    for k in ('p_node', 'u_fcs', 'u_fes', 'u_fev'): x[k] = r(win(o[k], per=per))
    x['sbp'] = r(win(o['pa_m'], np.max, per)); x['dbp'] = r(win(o['pa_m'], np.min, per))
    for k in ('u_rsp', 'u_rep', 'u_elv', 'u_erv', 'u_ts', 'u_tv'): x[k] = r(inst(o[k]), 4)
    for k in ('u_vusv', 'u_vuev'): x[k] = r(inst(o[k]), 1)
    res['cm2'] = x

    out = {k: d[k] for k in ('t', 'toco', 'tb', 'bins', 'peaks', 'p_ut')}; out.update(res)
    (HERE / 'traces.json').write_text(json.dumps(out, separators=(',', ':')))
    for m in ('cm0', 'cm1', 'cm2'):
        x = res[m]
        print(m, 'FHR rest %.1f min %.1f | MHR rest %.1f range %.1f-%.1f' % (x['fhr'][0], min(x['fhr']), x['mhr'][0], min(x['mhr']), max(x['mhr'])))


if __name__ == '__main__':
    main()
