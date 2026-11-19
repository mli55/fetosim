"""Run the pushing configuration on one measured Valsalva strain and save what docs/pushing_run.html shows.

Input: expiratory pressure of PhysioNet rvmh1 subject 27, strain at 622.4 s for 16.9 s (mean 40 mmHg; the effort
closest to 40 mmHg for 15 s, chosen before looking at the response), 1-s mean, negatives set to 0. It is the thoracic
pressure above rest and is also added to the uterine tone (the strain presses on the uterus; no contraction).
Model: Jongen 2016 mother + Lu 2001 thorax + Randall 2019 heart-rate reflex (closed loop) + Hill 1972 placenta,
default parameters, nothing fitted to the subject.
Writes traces.json (every 10th 0.02 s step of the plotted window, per block) and outputs.png (the four outputs).
Run from the repo root: RVMH1_DIR=... python docs/pushing_run/run.py, then python docs/pushing_run/build.py"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); import rvmh1 as R                      # noqa: E402
sys.path.insert(0, str(HERE.parents[2])); from fetosim.model import simulate as S  # noqa: E402

SUBJECT, START, DUR = 27, 622.4, 16.9                                 # the strain (s in the record)
PRE = POST = 60.0                                                     # simulated before / after it
a, c = START, START + DUR; t0, t1 = a - 20, c + 45                    # plotted window

p, fp = R.read('norm%03d_P' % SUBJECT)
PTH0 = float(S.params()[S.IDX['MR_PTH_REST']])
i0, i1 = int(round((a - PRE) * fp)), int(round((c + POST) * fp))
n = int(round(fp)); sm = np.convolve(np.pad(p, (n // 2, n - 1 - n // 2), mode='edge'), np.ones(n) / n, mode='valid')
dp = np.clip(sm[i0:i1], 0, None)
o, _, ok = S.run(S.TONE + dp, du=1 / fp, mhr=2, thorax=1, placenta=1, pth=PTH0 + dp, out_every=1); assert ok

dt = S.DTS; t = np.arange(len(o['mhr'])) * dt + i0 / fp
P = S.params(MR_MODE=2, MT_MODE=1, H_MODE=1); g = lambda nm: float(P[S.IDX[nm]])

# Randall 2019 internals (not outputs of core.c), replayed with core.c's Euler scheme from the recorded SBP, P_th, T_s
sbp, pth = o['sbp_m'], o['pth_m']
ws = lambda x: 1 - np.sqrt((1 + np.exp(-g('MR_QW') * (x - g('MR_SW')))) / (g('MR_A') + np.exp(-g('MR_QW') * (x - g('MR_SW')))))
ewc, ewa = ws(sbp), ws(sbp - pth)
N = len(t); ebc, eba, tpb, tpr, H = (np.zeros(N) for _ in range(5)); ts = o['ts_m']
ebc[0], eba[0] = g('MR_KB') * ewc[0], g('MR_KB') * ewa[0]
nb0 = g('MR_B') * (ewc[0] - ebc[0]) + (1 - g('MR_B')) * (ewa[0] - eba[0])
tpb[0] = g('MR_KPB') / (1 + np.exp(-g('MR_QPB') * (nb0 - g('MR_SPB'))))
tpr[0] = g('MR_KPR') / (1 + np.exp(g('MR_QPR') * (pth[0] - g('MR_SPR'))))
H[0] = o['mhr'][0]; nb = np.zeros(N)
for i in range(N):
    nb[i] = g('MR_B') * (ewc[i] - ebc[i]) + (1 - g('MR_B')) * (ewa[i] - eba[i])
    if i == N - 1:
        break
    gpb = 1 / (1 + np.exp(-g('MR_QPB') * (nb[i] - g('MR_SPB')))); gpr = 1 / (1 + np.exp(g('MR_QPR') * (pth[i] - g('MR_SPR'))))
    htil = g('MR_HI') * (1 - g('MR_HPB') * tpb[i] + g('MR_HPR') * tpr[i] + g('MR_HS') * ts[i])
    ebc[i + 1] = ebc[i] + dt * (-ebc[i] + g('MR_KB') * ewc[i]) / g('MR_TAUB')
    eba[i + 1] = eba[i] + dt * (-eba[i] + g('MR_KB') * ewa[i]) / g('MR_TAUB')
    tpb[i + 1] = tpb[i] + dt * (-tpb[i] + g('MR_KPB') * gpb) / g('MR_TAUPB')
    tpr[i + 1] = tpr[i] + dt * (-tpr[i] + g('MR_KPR') * gpr) / g('MR_TAUPR')
    H[i + 1] = H[i] + dt * (-H[i] + htil) / g('MR_TAUH')
assert np.abs(H - o['mhr']).max() < 1e-9, 'Randall replay does not match core.c'

b, fb = R.read('norm%03d_BP' % SUBJECT); tb, msbp, mhr = R.beats(b, fb)
sel = (t >= t0) & (t <= t1); idx = np.flatnonzero(sel)[::10]
base = lambda tt, y: y - np.nanmedian(y[(tt > a - 10) & (tt < a)])


def smo(y, w=4.0):
    m = int(w / dt); return np.convolve(np.pad(y, (m // 2, m - 1 - m // 2), mode='edge'), np.ones(m) / m, mode='valid')


qd = np.full(len(t), np.nan)                                          # placental transfer: same C function as core.c
for i_ in idx:
    qd[i_] = S.hill_exchange(o['po2_m_a'][i_], o['po2_a'][i_], o['q_uta'][i_], o['q_uma'][i_], H_DM=g('H_DMC'))[0] * 60
tr = dict(dpth=o['pth_m'] - PTH0, p_ut=o['p_ut'], pa_m=o['pa_m'], pv_m=o['pv_m'], pth_tm=o['pth_tm_m'], q_vr=smo(o['q_vr_m']),
          q_co=smo(o['q_co_m']), sbp=o['sbp_m'], dsbp=base(t, o['sbp_m']), ewc=ewc, ewa=ewa, nb=nb, tpb=tpb, tpr=tpr, ts=ts,
          mhr=o['mhr'], dmhr=base(t, o['mhr']), q_uta=smo(o['q_uta']), q_uma=smo(o['q_uma']), po2_a=o['po2_a'], r_c=o['r_c'],
          z_ts=o['z_ts'], z_tv=o['z_tv'], fhr=o['fhr'],
          spo2_m=o['spo2_m'], sto2_m=o['sto2_m'], spo2_f=o['spo2_f'], sto2_f_tis=o['sto2_f_tis'], sto2_f_cer=o['sto2_f_cer'],
          po2_m_a=o['po2_m_a'], po2_m_tis=o['po2_m_tis'], po2_tis=o['po2_tis'], po2_cer=o['po2_cer'], pa_f=o['ptm_a'],
          q_co_f=smo(o['q_co']), r_b=o['r_b'], e_v=o['e_v'], e_s=o['e_s'], qd_pl=qd)
out = {'t': [round(float(x), 2) for x in t[idx]], 'win': [t0, t1], 'strain': [a, c],
       **{kk: [round(float(x), 5) for x in v[idx]] for kk, v in tr.items()}}
ms = (tb >= t0) & (tb <= t1)
out['meas'] = {'t': [round(float(x), 2) for x in tb[ms]], 'dsbp': [round(float(x), 2) for x in base(tb, msbp)[ms]],
               'dmhr': [round(float(x), 2) for x in base(tb, mhr)[ms]]}
pi = np.arange(int(t0 * fp), int(t1 * fp)); out['meas']['tp'] = [round(float(x), 2) for x in (pi / fp)[::5]]
out['meas']['p'] = [round(float(x), 2) for x in p[pi][::5]]
out['par'] = {nm: g(nm) for nm in ('MT_FA', 'MT_FVAB', 'MT_KAB', 'MT_RVR', 'MT_CTH', 'MR_PTH_REST', 'MR_SBP_OFF', 'MR_B', 'MR_DS',
                                     'MR_HI', 'M_RUTA', 'MR_TAUB', 'MR_TAUPB', 'MR_TAUPR', 'MR_TAUS', 'MR_TAUH', 'MR_SW', 'MR_KB',
                                     'H_DM', 'H_DMC', 'F_TREF')}
(HERE / 'traces.json').write_text(json.dumps(out))

# the four outputs: thoracic pressure in; maternal heart rate (with the measured one), FHR, SpO2, StO2 out
plt.rcParams.update({'font.family': ['Avenir Next LT Pro', 'Avenir', 'Arial'], 'font.size': 13, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': .8})
INK, MEAS = '#1f2a44', '#23838a'; tp = np.arange(i0, i1) / fp
panels = [('Thoracic pressure\n(mmHg above rest)', [(np.arange(len(p)) / fp, p, MEAS, 1.0), (tp, dp, INK, 1.6)]),
          ('Maternal heart rate,\nchange (bpm)', [(tb, base(tb, mhr), MEAS, 1.6), (t, base(t, o['mhr']), INK, 1.6)]),
          ('FHR (bpm)', [(t, o['fhr'], INK, 1.6)]),
          ('Maternal SpO$_2$ (%)', [(t, o['spo2_m'], INK, 1.6)]),
          ('Fetal SpO$_2$ (%)', [(t, o['spo2_f'], INK, 1.6)]),
          ('Maternal StO$_2$ (%)', [(t, o['sto2_m'], INK, 1.6)]),
          ('Fetal cerebral\nStO$_2$ (%)', [(t, o['sto2_f_cer'], INK, 1.6)])]
fig, ax = plt.subplots(len(panels), 1, figsize=(7.6, 1.15 * len(panels) + .7), sharex=True)
for a_, (lab, lines) in zip(ax, panels):
    a_.axvspan(a, c, color='#f1f2f5', lw=0); a_.yaxis.set_major_locator(plt.MaxNLocator(3)); a_.tick_params(length=3)
    for x, y, col, lw in lines:
        a_.plot(x, y, color=col, lw=lw)
    ys = np.concatenate([y[(x > t0) & (x < t1)] for x, y, _, _ in lines]); ys = ys[np.isfinite(ys)]
    pad = max(.08 * (ys.max() - ys.min()), .2); a_.set_ylim(ys.min() - pad, ys.max() + pad)
    a_.set_ylabel(lab, rotation=0, ha='right', va='center', labelpad=8)
ax[-1].set_xlim(t0, t1); ax[-1].set_xlabel('Time (s)'); fig.align_ylabels(ax); fig.tight_layout(h_pad=.5)
fig.savefig(HERE / 'outputs.png', dpi=220); plt.close(fig)
print('wrote traces.json, outputs.png')
