"""Build docs/pushing_run.html from traces.json and outputs.png (written by run.py) and the current model code:
the call, then every block of S.run() with its code lines, inputs, outputs and traces.
Run from the repo root: python docs/pushing_run/build.py"""
import base64
import html as H
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
tr = json.loads((HERE / 'traces.json').read_text())
png = base64.b64encode((HERE / 'outputs.png').read_bytes()).decode()
CORE = (REPO / 'model/core.c').read_text().splitlines()
SIMPY = (REPO / 'model/simulate.py').read_text().splitlines()
DRV = (HERE / 'run.py').read_text().splitlines()
GH = 'https://github.com/mli55/fetosim/blob/main/model/core.c'
T0, T1 = tr['win']; A, C = tr['strain']
PAR = tr['par']


# ---------------------------------------------------------------- sparklines
def xs(t):
    return 3 + 294 * (t - T0) / (T1 - T0)


def spark(series, lo=None, hi=None, label=''):
    """series: list of (t list, y list, css class)."""
    ys = [v for _, y, _ in series for v in y]
    lo = min(ys) if lo is None else lo
    hi = max(ys) if hi is None else hi
    if hi - lo < 1e-9:
        hi = lo + 1
    pad = 4
    Y = lambda v: pad + (46 - 2 * pad) * (1 - (v - lo) / (hi - lo))
    lines = ''.join('<polyline class="ln %s" points="%s"/>' % (c, ' '.join('%.1f,%.1f' % (xs(a), Y(b)) for a, b in zip(t, y)))
                    for t, y, c in series)
    band = '<rect class="band" x="%.1f" y="0" width="%.1f" height="46"/>' % (xs(A), xs(C) - xs(A))
    span = hi - lo
    f = (lambda v: '%.0f' % v) if span >= 8 else (lambda v: '%.1f' % v) if span >= 0.8 else (lambda v: '%.2f' % v)
    return ('<div class="sp"><svg viewBox="0 0 300 46" preserveAspectRatio="none" role="img" aria-label="%s">%s%s</svg>'
            '<span class="rg"><span>%s</span><span>%s</span></span></div>' % (label, band, lines, f(hi), f(lo)))


def S1(key, cls='out'):
    return spark([(tr['t'], tr[key], cls)], label=key)


M = tr['meas']


def row(kind, name, desc, link='', sp=''):
    tag = {'in': 'Input', 'loop': 'Loop', 'out': 'Output', 'given': 'Given'}[kind]
    return ('<div class="r %s"><span class="tag">%s</span><code class="nm">%s</code><span class="ds">%s</span>'
            '<span class="lk">%s</span>%s</div>' % (kind, tag, name, desc, link, sp or '<span></span>'))


def code(src, ranges):
    out = []
    for a, b in ranges:
        for n in range(a, b + 1):
            out.append('<span class="l"><i>%d</i>%s</span>' % (n, H.escape(src[n - 1].rstrip())))
    return '<pre class="code">' + ''.join(out) + '</pre>'


def L(n):
    return '<a href="#b%d">← %d</a>' % (n, n)


def F(n):
    return '<a href="#b%d">→ %d</a>' % (n, n)


def card(n, title, fn, ranges, rows):
    first = ranges[0][0]
    key = {8: 'R_TPB', 10: 'VM_ART', 15: 'rcs'}.get(n) or re.match(r'[A-Za-z_][A-Za-z_0-9]*', fn).group(0)
    # the code shown must contain the block's variable, so a shift in core.c line numbers stops the build
    assert any(key in CORE[i - 1] for a, b in ranges for i in range(a, b + 1)), (n, key)
    return ('<section class="step" id="b%d"><header><span class="n">%d</span><span class="tt">%s</span><code class="fn">%s</code>'
            '<a class="at" href="%s#L%d">core.c:%d</a></header>%s<div class="rows">%s</div></section>'
            % (n, n, title, fn, GH, first, first, code(CORE, ranges), ''.join(rows)))


# ---------------------------------------------------------------- the call
i_run = next(i for i, l in enumerate(DRV, 1) if 'S.run(S.TONE + dp' in l)
i_dp = next(i for i, l in enumerate(DRV, 1) if 'dp = np.clip' in l)
run_card = ('<section class="step" id="run"><header><span class="n">▶</span><span class="tt">The call</span>'
            '<code class="fn"><b>S.run</b>(S.TONE + dp, du=1/fp, mhr=2, thorax=1, placenta=1, pth=PTH0 + dp, out_every=1)</code>'
            '<a class="at" href="https://github.com/mli55/fetosim/blob/main/docs/pushing_run/run.py#L%d">run.py:%d</a></header>%s<div class="rows">%s</div></section>'
            % (i_run, i_run, code(DRV, [(i_dp - 2, i_dp), (i_run, i_run)]), ''.join([
                row('in', 'dp', 'thoracic pressure above rest, mmHg: measured expiratory pressure (PhysioNet rvmh1, subject 27), 1-s mean, negatives set to 0',
                    '→ 1', spark([(M['tp'], M['p'], 'meas'), (tr['t'], tr['dpth'], 'in')], lo=-2, label='dp')),
                row('given', 'pth', '= PTH0 + dp, PTH0 = %.1f mmHg resting thoracic pressure' % PAR['MR_PTH_REST']),
                row('given', 'p_ut', '= S.TONE + dp = 15 mmHg uterine tone + dp; the strain also presses on the uterus', '→ 6'),
                row('given', 'mhr=2', 'Randall 2019 heart-rate reflex, closed loop (blocks 7–9)'),
                row('given', 'thorax=1', 'Lu 2001 thorax mechanics lumped onto the Jongen 2016 mother (blocks 1–5)'),
                row('given', 'placenta=1', 'Hill 1972 placental O2 exchanger in place of Jongen\'s A(16) diffusion (block 12)'),
                row('out', 'mhr', 'maternal heart rate, change from the 10 s before the strain, bpm', '← 9',
                    spark([(M['t'], M['dmhr'], 'meas'), (tr['t'], tr['dmhr'], 'out')], label='MHR change')),
                row('out', 'fhr', 'fetal heart rate, bpm', '← 16', S1('fhr')),
                row('out', 'spo2_m', 'maternal SpO2 (arterial O2 saturation), %', '← 17', S1('spo2_m')),
                row('out', 'spo2_f', 'fetal SpO2, %', '← 17', S1('spo2_f')),
                row('out', 'sto2_m', 'maternal StO2 (tissue O2 saturation, NIRS), %', '← 17', S1('sto2_m')),
                row('out', 'sto2_f_cer', 'fetal cerebral StO2, %', '← 17', S1('sto2_f_cer')),
            ])))

# ---------------------------------------------------------------- block cards
cards = [
    card(1, 'Thoracic excess', 'dpth', [(180, 182)], [
        row('in', 'pth', 'thoracic pressure, mmHg', '← call'),
        row('given', 'MR_PTH_REST', '= %.1f mmHg, resting thoracic pressure' % PAR['MR_PTH_REST']),
        row('out', 'dpth', 'thoracic pressure above rest, mmHg', '→ 2 3 4 5', S1('dpth'))]),
    card(2, 'Arterial pressure', 'pa_m', [(183, 183)], [
        row('loop', 'V_art', 'arterial volume, ml, from the previous step', L(10)),
        row('in', 'dpth', 'thoracic excess, mmHg', L(1)),
        row('given', 'MT_FA', '= %.3f, share of the arterial compliance inside the chest (Lu 2001; fitted to Fox phase-1 SBP +18 mmHg)' % PAR['MT_FA']),
        row('out', 'pa_m', 'maternal arterial pressure, mmHg', '→ 5 6 7', S1('pa_m'))]),
    card(3, 'Venous pressure', 'pv_m', [(188, 188)], [
        row('loop', 'V_ven', 'venous volume, ml, from the previous step', L(10)),
        row('in', 'dpth', 'thoracic excess, mmHg', L(1)),
        row('given', 'MT_FVAB', '= %.3f, share of the venous compliance in the abdomen, pressed by the strain (our addition; fitted to Fox phase-2 CO 39 %%)' % PAR['MT_FVAB']),
        row('out', 'pv_m', 'maternal venous pressure, mmHg', '→ 4 6', S1('pv_m'))]),
    card(4, 'Venous return', 'q_vr', [(197, 198)], [
        row('loop', 'V_TH', 'intrathoracic venous and pulmonary blood, ml', L(10)),
        row('in', 'pv_m', 'venous pressure, mmHg', L(3)),
        row('in', 'dpth', 'thoracic excess, mmHg', L(1)),
        row('given', 'MT_CTH', '= %.1f ml/mmHg, compliance of the intrathoracic reservoir (Lu 2001: vena cava + pulmonary)' % PAR['MT_CTH']),
        row('given', 'MT_RVR', '= %.3f mmHg·s/ml, vena cava resistance (Lu 2001); no backflow' % PAR['MT_RVR']),
        row('out', 'pth_tm', 'transmural pressure of the intrathoracic reservoir, mmHg', '→ 5', S1('pth_tm')),
        row('out', 'q_vr', 'venous return into the chest, ml/s (4-s mean)', '→ 10', S1('q_vr'))]),
    card(5, 'Heart', 'qco_m', [(199, 204)], [
        row('in', 'pth_tm', 'filling pressure, mmHg', L(4)),
        row('in', 'pa_m', 'arterial pressure, mmHg', L(2)),
        row('in', 'dpth', 'the heart ejects against pa_m − dpth', L(1)),
        row('loop', 'H', 'heart rate, bpm; heart period t_m = 60 / H', '<a href="#b9">← 9</a>'),
        row('given', 'M_EMIN, M_EMAX', 'Jongen 2016 heart: end-diastolic and end-systolic volume from the pressures'),
        row('out', 'qco_m', 'cardiac output, ml/s (4-s mean)', '→ 10', S1('q_co'))]),
    card(6, 'Uterine inflow', 'quta', [(190, 193), (208, 208)], [
        row('in', 'pa_m', 'arterial pressure, mmHg', L(2)),
        row('in', 'pv_m', 'venous pressure, mmHg', L(3)),
        row('in', 'p_ut', 'uterine pressure = 15 + dp, mmHg; compresses the uterine veins', '← call', S1('p_ut', 'in')),
        row('loop', 'V_IVS', 'intervillous volume, ml', L(10)),
        row('given', 'M_RUTA', '= %.2f mmHg·s/ml, uterine arterial resistance (Jongen 2016)' % PAR['M_RUTA']),
        row('out', 'quta', 'uterine blood flow, ml/s (4-s mean)', '→ 10 11 12', S1('q_uta'))]),
    card(7, 'Baroreceptors', 'nb_r', [(319, 323), (330, 331)], [
        row('in', 'pa_m', 'SBP = pa_m + MR_SBP_OFF (%.1f mmHg); Jongen\'s heart is not pulsatile' % PAR['MR_SBP_OFF'], L(2)),
        row('in', 'pthv', 'thoracic pressure; the aortic receptor sees SBP − P_th', '← call'),
        row('out', 'sbp', 'systolic pressure, mmHg', '', S1('sbp')),
        row('out', 'ewc', 'carotid wall strain, from SBP', '', S1('ewc')),
        row('out', 'ewa', 'aortic wall strain, from SBP − P_th', '', S1('ewa')),
        row('out', 'nb_r', 'afferent signal: strain minus its slow adapting part, carotid and aortic weighted by MR_B = %.1f' % PAR['MR_B'], '→ 8', S1('nb'))]),
    card(8, 'Autonomic outflow', 'T_pb, T_pr, T_s', [(324, 325), (328, 328), (332, 334)], [
        row('in', 'nb_r', 'afferent signal', L(7)),
        row('in', 'pthv', 'thoracic pressure, drives the respiratory parasympathetic part', '← call'),
        row('given', 'MR_TAU*', 'time constants %.1f / %.1f / %.1f s; sympathetic delay MR_DS = %.0f s' % (PAR['MR_TAUPB'], PAR['MR_TAUPR'], PAR['MR_TAUS'], PAR['MR_DS'])),
        row('out', 'T_pb', 'parasympathetic tone from the baroreflex', '→ 9', S1('tpb')),
        row('out', 'T_pr', 'parasympathetic tone from thoracic pressure', '→ 9', S1('tpr')),
        row('out', 'T_s', 'sympathetic tone', '→ 9', S1('ts'))]),
    card(9, 'Heart rate', 'H', [(326, 326), (335, 335)], [
        row('in', 'T_pb, T_pr, T_s', 'autonomic tones', L(8)),
        row('given', 'MR_HI', '= %.1f bpm intrinsic heart rate; MR_TAUH = %.1f s' % (PAR['MR_HI'], PAR['MR_TAUH'])),
        row('out', 'H', 'maternal heart rate, bpm', '→ 5, out', S1('mhr'))]),
    card(10, 'Volume update', 'Euler, dt = 0.02 s', [(351, 354), (359, 359)], [
        row('in', 'qco_m', 'cardiac output', L(5)),
        row('in', 'quta', 'uterine inflow', L(6)),
        row('in', 'q_vr', 'venous return', L(4)),
        row('out', 'V_art, V_ven, V_TH, V_IVS', 'volumes for the next step', '→ 2 3 4 6 11')]),
    card(11, 'Maternal O2 pools', 'pg[0], pg[2]', [(207, 207), (263, 265), (268, 272), (282, 283), (286, 291), (355, 360)], [
        row('in', 'qco_m', 'cardiac output', L(5)),
        row('in', 'quta, qutv', 'uterine inflow and outflow', L(6)),
        row('in', 'q_vr', 'venous return into the chest', L(4)),
        row('in', 'qtisa_m, qtisv_m', 'maternal tissue flows, line 207 (from pa_m, pv_m and tissue pressure)'),
        row('loop', 'qd_pl', 'O2 leaving the intervillous pool into the fetus', '<a href="#b12">← 12</a>'),
        row('loop', 'V_*, MM_*', 'blood volumes and O2 contents of the maternal pools, from the previous step', L(10)),
        row('given', 'O_DLUNG, O_PAIR', 'lung uptake qd_lu = O_DLUNG (O_PAIR − pg[1]) into the venous pool, A(16)'),
        row('given', 'O_QTISM, O_THTISM', 'maternal tissue O2 use, reduced below a tissue pO2 threshold, A(17)'),
        row('out', 'pg[0]', 'maternal arterial pO2, mmHg', '→ 12 17', S1('po2_m_a')),
        row('out', 'pg[2]', 'maternal tissue pO2, mmHg', '→ 17', S1('po2_m_tis'))]),
    card(12, 'Placental O2 transfer', 'qd_pl', [(281, 281)], [
        row('in', 'pg[0]', 'maternal arterial pO2 = pO2 of the blood entering the exchanger, mmHg', L(11)),
        row('in', 'pg[5]', 'fetal arterial pO2 = pO2 of the fetal blood entering, mmHg', '<a href="#b14">← 14</a>'),
        row('in', 'quta', 'maternal capillary flow, ml/s (4-s mean)', L(6), S1('q_uta', 'in')),
        row('in', 'quma', 'fetal capillary flow (umbilical arterial), ml/s (4-s mean)', '<a href="#b13">← 13</a>', S1('q_uma', 'in')),
        row('given', 'H_DMC', '= %.2f ml/(min·mmHg), membrane diffusing capacity in the model; calibrated so that fetal arterial pO2 is 18.5 mmHg at rest (Hill 1972 Table 1: %.2f)' % (PAR['H_DMC'], PAR['H_DM'])),
        row('given', 'H_K*, H_HB*, H_VM, H_VF', 'Hill 1972 human dissociation curves, hemoglobin, capillary volumes 10 ml'),
        row('out', 'qd_pl', 'O2 transfer, ml/min: out of the intervillous pool, into the villous pool (lines 357, 369)', '→ 11 14', S1('qd_pl'))]),
    card(13, 'Fetal circulation', 'quma, pa', [(246, 260)], [
        row('in', 'T', 'fetal heart period, s', '<a href="#b16">← 16</a>'),
        row('in', 'zE, zR, zV', 'effector signals: heart contractility, tissue resistance, venous unstressed volume', '<a href="#b16">← 16</a>'),
        row('loop', 'R_cera', 'cerebral arterial resistance (autoregulation)', '<a href="#b15">← 15</a>'),
        row('loop', 'V_F*', 'fetal blood volumes, from the previous step (update lines 361–365)'),
        row('given', 'p_um', '= 0: S_WUM = 0 in this run, so the cord is not compressed; p_ut cancels in the fetal transmural pressures'),
        row('out', 'quma', 'umbilical arterial flow, ml/s (4-s mean)', '→ 12 14', S1('q_uma')),
        row('out', 'qco, qtisa, qcera', 'fetal cardiac output (shown, 4-s mean) and tissue and cerebral flows, ml/s', '→ 14', S1('q_co_f')),
        row('out', 'pa', 'fetal arterial pressure, mmHg', '→ 15', S1('pa_f'))]),
    card(14, 'Fetal O2 pools', 'pg[5], pg[7], pg[8]', [(266, 267), (273, 278), (284, 285), (298, 301), (366, 371)], [
        row('in', 'qd_pl', 'O2 entering the villous pool, ml/min', '<a href="#b12">← 12</a>'),
        row('in', 'qco, quma, qtisa, qcera', 'fetal flows that carry O2 between the pools', '<a href="#b13">← 13</a>'),
        row('loop', 'MF_*', 'O2 contents of the fetal pools, from the previous step'),
        row('given', 'O_QTISF, O_QCERF', 'fetal tissue and cerebral O2 use, A(17)'),
        row('out', 'pg[5]', 'fetal arterial pO2, mmHg', '→ 12 15 17', S1('po2_a')),
        row('out', 'pg[7], pg[8]', 'fetal tissue and cerebral pO2, mmHg (cerebral shown)', '→ 17', S1('po2_cer'))]),
    card(15, 'Fetal reflexes', 'r_b, r_c, e_v, e_s', [(304, 316), (373, 380)], [
        row('in', 'pg[5]', 'fetal arterial pO2, mmHg', '<a href="#b14">← 14</a>'),
        row('in', 'pa', 'fetal arterial pressure, mmHg', '<a href="#b13">← 13</a>'),
        row('given', 'G_*', 'Jongen 2016 baro- and chemoreceptor curves (Gompertz), weights, delays and time constants; noise N = 0'),
        row('out', 'r_c', 'chemoreceptor activity', '', S1('r_c')),
        row('out', 'r_b', 'baroreceptor activity', '', S1('r_b')),
        row('out', 'e_v', 'vagal efferent', '', S1('e_v')),
        row('out', 'e_s', 'sympathetic efferent', '', S1('e_s')),
        row('loop', 'x_Tv, x_Ts, x_E, x_R, x_V', 'delayed, filtered efferents (state), read in the next step', '→ 16'),
        row('loop', 'R_cera', 'cerebral resistance from fetal arterial O2 content, A(18), A(19)', '→ 13')]),
    card(16, 'Fetal heart rate', 'T', [(243, 245)], [
        row('loop', 'x_Tv, x_Ts, x_E, x_R, x_V', 'filtered efferents from the previous step', '<a href="#b15">← 15</a>'),
        row('given', 'F_TREF', '= %.3f s, resting fetal heart period; S_GAIN = 1' % PAR['F_TREF']),
        row('out', 'T', 'fetal heart period; FHR = 60 / T, bpm', '→ 13, out', S1('fhr')),
        row('out', 'zE, zR, zV', 'effector signals', '→ 13')]),
    ('<section class="step" id="b17"><header><span class="n">17</span><span class="tt">O2 saturation</span>'
     '<code class="fn">SpO2, StO2</code><a class="at" href="%s#L203">simulate.py:203</a></header>'
     '<p class="note">Not part of the 0.02 s loop: run_core computes it from the recorded pO2 after core.c returns.</p>%s<div class="rows">%s</div></section>'
     % (GH.replace('core.c', 'simulate.py'), code(SIMPY, [(75, 75), (203, 216)]), ''.join([
        row('in', 'pg[0], pg[2]', 'maternal arterial and tissue pO2, mmHg', L(11)),
        row('in', 'pg[5], pg[7], pg[8]', 'fetal arterial, tissue and cerebral pO2, mmHg', '<a href="#b14">← 14</a>'),
        row('given', 'O_C1M, O_C1F, O_C2*', 'Jongen 2016 dissociation curves, A(13): saturation = g / (g + c1), g = p³ + c2·p'),
        row('given', 'STO2_ART', '= 0.25, arterial share of the NIRS volume (our choice: the usual 25 : 75 arterial : venous)'),
        row('out', 'spo2_m', 'maternal arterial saturation, %', 'out', S1('spo2_m')),
        row('out', 'sto2_m', 'maternal tissue saturation, %', 'out', S1('sto2_m')),
        row('out', 'spo2_f', 'fetal arterial saturation, %', 'out', S1('spo2_f')),
        row('out', 'sto2_f_tis', 'fetal tissue saturation, %', 'out', S1('sto2_f_tis')),
        row('out', 'sto2_f_cer', 'fetal cerebral saturation, %', 'out', S1('sto2_f_cer'))]))),
]

# ---------------------------------------------------------------- block diagram
BW, BH = 176, 62
B = {
    1: (190, 160, 'Thoracic excess', 'dpth', ''),
    2: (420, 44, 'Arterial pressure', 'pa_m', ''),
    3: (420, 160, 'Venous pressure', 'pv_m', ''),
    4: (420, 276, 'Venous return', 'q_vr', ''),
    7: (660, 44, 'Baroreceptors', 'nb_r', ''),
    5: (660, 160, 'Heart', 'qco_m', ''),
    6: (660, 276, 'Uterine inflow', 'quta', ''),
    8: (900, 44, 'Autonomic outflow', 'T_pb, T_pr, T_s', ''),
    10: (900, 218, 'Volume update', 'V_art, V_ven, V_TH', ''),
    9: (1180, 44, 'Heart rate', 'H', 'out: MHR'),
    11: (660, 440, 'Maternal O2 pools', 'pg[0], pg[2]', ''),
    12: (900, 440, 'Placental O2', 'qd_pl (Hill 1972)', ''),
    14: (1180, 440, 'Fetal O2 pools', 'pg[5], pg[7], pg[8]', ''),
    17: (420, 590, 'O2 saturation', 'after the run', 'out: SpO2, StO2'),
    16: (660, 590, 'Fetal heart rate', 'T', 'out: FHR'),
    13: (900, 590, 'Fetal circulation', 'quma, pa', ''),
    15: (1180, 590, 'Fetal reflexes', 'r_b, r_c, e_v, e_s', ''),
}
g = ['<svg class="dia" viewBox="0 0 1390 730" role="img" aria-label="Blocks inside S.run for the pushing run">',
     '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
     '<path d="M0,0 L10,5 L0,10 z" class="ahf"/></marker>'
     '<marker id="ahl" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
     '<path d="M0,0 L10,5 L0,10 z" class="ahl"/></marker></defs>']


def path(pts, loop=False):
    d = 'M' + ' L'.join('%d,%d' % p for p in pts)
    g.append('<path class="%s" d="%s" marker-end="url(#%s)"/>' % ('ed lp' if loop else 'ed', d, 'ahl' if loop else 'ah'))


def lab(x, y, s, loop=False, anchor='start'):
    g.append('<text class="el%s" x="%d" y="%d" text-anchor="%s">%s</text>' % (' lp' if loop else '', x, y, anchor, s))


R_ = lambda n: B[n][0] + BW
CY = lambda n: B[n][1] + BH // 2
X_ = lambda n: B[n][0]
Y_ = lambda n: B[n][1]
# input boxes
g.append('<rect class="inb" x="20" y="160" width="130" height="62" rx="8"/><text class="it" x="85" y="187">p_th(t)</text>'
         '<text class="is" x="85" y="206">thoracic pressure</text>')
g.append('<rect class="inb" x="20" y="300" width="130" height="50" rx="8"/><text class="it" x="85" y="322">p_ut(t)</text>'
         '<text class="is" x="85" y="339">15 + dp</text>')
# edges
path([(150, CY(1)), (X_(1), CY(1))])
path([(85, 160), (85, 18), (700, 18), (700, Y_(7))]); path([(700, 18), (940, 18), (940, Y_(8))]); lab(300, 13, 'p_th')
sx = 392
path([(R_(1), CY(1)), (sx, CY(1)), (sx, CY(2) - 8), (X_(2), CY(2) - 8)])
path([(sx, CY(3) - 8), (X_(3), CY(3) - 8)]); path([(sx, CY(1)), (sx, CY(4) - 8), (X_(4), CY(4) - 8)])
path([(sx, 133), (720, 133), (720, Y_(5))]); lab(372, 128, 'dpth')
path([(R_(2), CY(2)), (X_(7), CY(2))]); lab(R_(2) + 8, CY(2) - 6, 'pa_m')
path([(616, CY(2)), (616, CY(5) - 10), (X_(5), CY(5) - 10)]); path([(616, CY(5) - 10), (616, CY(6) - 12), (X_(6), CY(6) - 12)])
path([(R_(3), CY(3) + 12), (632, CY(3) + 12), (632, CY(6) + 4), (X_(6), CY(6) + 4)]); lab(R_(3) + 6, CY(3) + 28, 'pv_m')
path([(508, Y_(3) + BH), (508, Y_(4))]); lab(514, Y_(4) - 22, 'pv_m')
path([(R_(4), CY(4) - 18), (604, CY(4) - 18), (604, CY(5) + 12), (X_(5), CY(5) + 12)]); lab(R_(4) + 6, CY(4) - 24, 'pth_tm')
path([(150, 325), (180, 325), (180, 398), (740, 398), (740, Y_(6) + BH)]); lab(190, 393, 'p_ut')
path([(R_(7), CY(7)), (X_(8), CY(7))]); lab(R_(7) + 8, CY(7) - 6, 'nb_r')
path([(R_(8), CY(8)), (X_(9), CY(8))]); lab(R_(8) + 8, CY(8) - 6, 'tones')
path([(R_(5), CY(5)), (872, CY(5)), (872, CY(10) - 14), (X_(10), CY(10) - 14)]); lab(R_(5) + 6, CY(5) - 6, 'qco_m')
path([(R_(6), CY(6)), (872, CY(6)), (872, CY(10) + 14), (X_(10), CY(10) + 14)]); lab(R_(6) + 6, CY(6) - 6, 'quta')
path([(508, Y_(4) + BH), (508, 360), (988, 360), (988, Y_(10) + BH)]); lab(516, 355, 'q_vr')
# loops (dashed)
path([(1268, Y_(9) + BH), (1268, 140), (780, 140), (780, Y_(5))], loop=True); lab(1000, 135, 'H', loop=True)
path([(R_(10), CY(10)), (1110, CY(10)), (1110, 378), (404, 378), (404, CY(4) + 10), (X_(4), CY(4) + 10)], loop=True)
path([(404, CY(4) + 10), (404, CY(3) + 10), (X_(3), CY(3) + 10)], loop=True)
path([(404, CY(3) + 10), (404, CY(2) + 10), (X_(2), CY(2) + 10)], loop=True)
lab(780, 373, 'V_art, V_ven, V_TH, V_IVS', loop=True)
path([(720, 378), (720, Y_(6) + BH)], loop=True)
path([(930, Y_(10) + BH), (930, 412), (800, 412), (800, Y_(11))]); lab(808, 430, 'flows, V_*')
path([(R_(6), CY(6) + 18), (856, CY(6) + 18), (856, 452), (X_(12), 452)]); lab(862, 446, 'quta')
path([(R_(11), 482), (X_(12), 482)]); lab(R_(11) + 6, 497, 'pg[0]')
path([(940, Y_(12) + BH), (940, 530), (760, 530), (760, Y_(11) + BH)], loop=True); lab(850, 545, 'qd_pl', loop=True, anchor='middle')
path([(R_(12), 461), (X_(14), 461)]); lab(R_(12) + 8, 455, 'qd_pl')
path([(X_(14), 487), (R_(12), 487)]); lab(R_(12) + 8, 500, 'pg[5]')
path([(990, Y_(13)), (990, Y_(12) + BH)]); lab(996, 572, 'quma')
path([(1060, Y_(13)), (1060, 562), (1220, 562), (1220, Y_(14) + BH)]); lab(1120, 556, 'fetal flows')
path([(1300, Y_(14) + BH), (1300, Y_(15))]); lab(1306, 572, 'pg[5]')
path([(R_(13), 640), (X_(15), 640)]); lab(R_(13) + 8, 634, 'pa')
path([(X_(15), 606), (R_(13), 606)], loop=True); lab(R_(13) + 8, 600, 'R_cera', loop=True)
path([(1268, Y_(15) + BH), (1268, 684), (748, 684), (748, Y_(16) + BH)], loop=True); lab(1000, 678, 'x_Tv, x_Ts, x_E, x_R, x_V', loop=True, anchor='middle')
path([(R_(16), 621), (X_(13), 621)]); lab(R_(16) + 10, 615, 'T, z')
path([(700, Y_(11) + BH), (700, 556), (508, 556), (508, Y_(17))]); lab(514, 574, 'pg[0], pg[2]')
path([(R_(14), 471), (1372, 471), (1372, 712), (470, 712), (470, Y_(17) + BH)]); lab(900, 706, 'pg[5], pg[7], pg[8]', anchor='middle')
# boxes
for n, (x, y, t, c, o) in B.items():
    g.append('<a href="#b%d"><rect class="bx" x="%d" y="%d" width="%d" height="%d" rx="8"/>'
             '<text class="bt" x="%d" y="%d">%d · %s</text><text class="bc" x="%d" y="%d">%s</text>'
             % (n, x, y, BW, BH, x + 12, y + 22, n, t, x + 12, y + 40, c))
    if o:
        g.append('<text class="bo" x="%d" y="%d">%s</text>' % (x + 12, y + 55, o))
    g.append('</a>')
g.append('</svg>')
dia = ''.join(g)

css = (HERE / 'style.css').read_text()
page = f'''<title>Pushing: Thoracic Pressure to MHR</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>{css}</style>
<main>
<h1>Pushing: thoracic pressure to MHR</h1>
<p class="sub">Built by <code>docs/pushing_run/run.py</code> and <code>build.py</code> from the current fetosim code. Input: expiratory pressure of PhysioNet rvmh1 subject 27 as thoracic pressure. Mother: Jongen 2016 circulation, Lu 2001 thorax, Randall 2019 heart-rate reflex. Placenta: Hill 1972 O2 exchanger. Fetus: Jongen 2016. Nothing fitted to this subject. Traces span {T0:.0f} to {T1:.0f} s of the record; the strain ({A:.0f}–{C:.0f} s) is shaded. Teal is measured.</p>
<div class="legend"><span><span class="tag" style="color:var(--in);background:var(--in-bg);font:600 11px var(--mono);padding:1px 6px;border-radius:4px">INPUT</span>read from an earlier block or the call</span><span><span class="tag" style="color:var(--loop);border:1px dashed var(--loop);font:600 11px var(--mono);padding:1px 6px;border-radius:4px">LOOP</span>returns from a later block, read from the previous step</span><span><span class="tag" style="color:var(--out);background:var(--out-bg);font:600 11px var(--mono);padding:1px 6px;border-radius:4px">OUTPUT</span>returned</span><span><span class="tag" style="color:var(--given);background:var(--given-bg);font:600 11px var(--mono);padding:1px 6px;border-radius:4px">GIVEN</span>fixed parameter</span></div>

{run_card}

<h2>Inside S.run(): blocks 1–16 run every 0.02 s in core.c; block 17 runs after</h2>
<div class="diabox">{dia}</div>
{''.join(cards)}

<h2>Result</h2>
<div class="vfig"><img alt="Thoracic pressure; maternal heart rate change, measured and simulated; simulated FHR, maternal and fetal SpO2, maternal StO2 and fetal cerebral StO2" src="data:image/png;base64,{png}"></div>
</main>
'''
(REPO / 'docs/pushing_run.html').write_text(page)
print(len(page))
