"""Jongen 2016 maternal-fetal CTG model: physiology, TOCO -> pressure input, signal cleaning, per-window fitting and
the parameter table.

Sources (abbreviations used in every docstring below):
  A   Jongen GJLM, van der Hout-van der Jagt MB, van de Vosse FN, Oei SG, Bovendeerd PHM. Simulation of fetal heart
      rate variations during labor, Part A: uterine flow reduction. J Biomech 49(12):2466-2473 (2016).
      doi:10.1016/j.jbiomech.2016.01.036
  B   Same authors. Part B: combined uterine and umbilical flow reduction. J Biomech 49(12):2474-2480 (2016).
      doi:10.1016/j.jbiomech.2016.01.046
  Th  Jongen GJLM, PhD thesis, Eindhoven University of Technology (2016), Ch. 4 (FHR variability noise).
  A(8) = Part A equation 8, A-T1a = Part A Table 1a, B(2) = Part B equation 2, and so on.

Layers (each function names the equations it implements):
  1. physiology   jongen.c + jongen.py, the forward-Euler core with all published parameters.
  2. signals      4 Hz CTG -> 2-s bins, optional "TOCO = 0 is a valid reading" rule, row masks, FHR artefact
                  cleaning, moving TOCO tone.
  3. input        TOCO -> uterine pressure p_ut(t). Default: clean pulses (15 mmHg tone + one Gaussian bell per
                  detected TOCO peak, height a x 70 mmHg). Two earlier mappings are kept for comparison.
  4. fit          a, w_um, g (grid, then bounded Nelder-Mead on the training MSE) and b (closed form) for one window.
  5. parameters   parameter_table(): every parameter with symbol, value/range, unit, meaning and class
                  (1 fixed at the paper's value, 2 fixed by our choice, 3 fitted per window).

Entry point for any recording given as 2-s bins (NaN = missing):
    from jongen2016 import model as M        # with the repo root on sys.path
    res = M.fit(fhr_bins, toco_bins)        # dict: a, w_um, g, b, fhr_model, p_ut, peaks_s, ...
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.signal import find_peaks, medfilt

try:
    from . import jongen as J
except ImportError:                     # run as a script: python jongen2016/model.py
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from jongen2016 import jongen as J

# ---- numerics and window layout -------------------------------------------------------------------------------
DT = 2                          # 2-s bins (8 samples at 4 Hz)
TRAIN, GAP, TEST = 600, 60, 300  # bins: 20 min fit, 2 min unused gap, 10 min comparison
NB = TRAIN + GAP + TEST
MEMORY = 60                     # bins: a row needs no TOCO gap in the previous 120 s
DTS = 0.02                      # model time step, s (A s2.4: forward Euler, dt 0.02 s)
DU = 0.5                        # input grid for p_ut, s (linearly interpolated inside the core)
OE = 5                          # core output every 5 steps = 0.1 s; averaged to 2-s bins
TONE = 15.0                     # resting uterine pressure p_rest, mmHg (A s2.4); also the tone of our pulse input
P_CON = 70.0                    # reference contraction amplitude above rest, mmHg (A s2.4; peak 85 mmHg)


# =================================================================================================================
# 1. physiology
# =================================================================================================================
_SS = {}


def steady():
    """Steady state of the published model at rest (A s2.4: run until |dp_tm,a/dt| and |dpO2,a/dt| < 1e-4 mmHg/s).
    Cached. w_um and g do not change the resting state (p_um = w_um (p_ut - p_rest) = 0 and all effector inputs
    are 0 at rest), so every fit starts from this one state."""
    if 'S' not in _SS:
        _SS['S'] = J.steady_state(J.params())
    return _SS['S']


def run(put, w_um=0.0, g=1.0, dt=DTS, du=DU, out_every=OE, **params):
    """Full core output (dict of traces at out_every*dt) for uterine pressure `put` sampled every `du` s.
    w_um is B(2) cord-compression weight (p_um = w_um (p_ut - p_rest)); g multiplies all five effector gains k of
    A(23) (our addition; g = 1 is the published regulation). Other core parameters by name (see jongen.PARAMS)."""
    P = J.params(S_WUM=w_um, S_GAIN=g, **params)
    return J.run(P, steady() if not params else J.steady_state(P), put, dt=dt, du=du, out_every=out_every)


def simulate(put, w_um, g, nb=NB):
    """Model FHR averaged over each 2-s bin (the observation of fit.py / audit.sim), or None if non-finite.
    FHR = 60 / T (A(25); the paper gives no explicit formula)."""
    P = J.params(S_WUM=w_um, S_GAIN=g)
    o, _, ok = J.run(P, steady(), put, dt=DTS, du=DU, out_every=OE)
    if not ok:
        return None
    f = o['fhr'][:nb * int(DT / (DTS * OE))].reshape(nb, -1).mean(1)
    return f if np.isfinite(f).all() else None


def contraction(t_total, dt, t_con=30.0, T_con=60.0, p_con=P_CON, p_rest=TONE):
    """A(8): p_ut = p_rest + p_con sin^2(pi (t - t_con) / T_con) on (t_con, t_con + T_con), else p_rest."""
    return J.contraction(t_total, dt, t_con=t_con, T_con=T_con, p_con=p_con, p_rest=p_rest)


# =================================================================================================================
# 2. signals
# =================================================================================================================
FLAT_S, FS = 60, 4              # TOCO flat-line rule (s), raw sampling rate (Hz)


def uc_valid(uc):
    """TOCO = 0 counts as a reading. Valid = finite and >= 0, except runs of one exactly constant value (0 or not)
    lasting >= FLAT_S s, which are treated as flat-line dropout / sensor off."""
    ok = np.isfinite(uc) & (uc >= 0)
    x = np.where(np.isfinite(uc), uc, np.nan)
    same = np.r_[False, x[1:] == x[:-1]]                       # sample equals the previous one
    start = np.flatnonzero(~same); ends = np.r_[start[1:], len(x)]
    for s, e in zip(start, ends):
        if e - s >= FLAT_S * FS and np.isfinite(x[s]):
            ok[s:e] = False
    return ok


def bin_signals(raw, uc_zero_valid=False):
    """raw: (n, 2) array of 4 Hz (FHR bpm, TOCO). 2-s bins of 8 samples; a bin needs >= 4 valid samples.
    FHR valid in 50-210 bpm. TOCO valid if > 0 (every 0 missing), or with uc_zero_valid by uc_valid().
    Returns (fhr_bins, toco_bins, fraction valid per bin and channel)."""
    count = raw.shape[0] // 8
    if not uc_zero_valid:
        x = raw[:count * 8].reshape(count, 8, 2)
        ok = np.isfinite(x)
        ok[:, :, 0] &= (x[:, :, 0] >= 50) & (x[:, :, 0] <= 210)
        ok[:, :, 1] &= x[:, :, 1] > 0
        n = ok.sum(axis=1)
        sums = np.where(ok, x, 0).sum(axis=1)
        out = np.full((count, 2), np.nan)
        np.divide(sums, n, out=out, where=n >= 4)
        return out[:, 0], out[:, 1], n / 8
    x = raw[:count * 8].copy()
    ok = np.zeros((len(x), 2), bool)
    ok[:, 0] = np.isfinite(x[:, 0]) & (x[:, 0] >= 50) & (x[:, 0] <= 210)
    ok[:, 1] = uc_valid(x[:, 1])
    x = x.reshape(count, 8, 2); ok = ok.reshape(count, 8, 2)
    n = ok.sum(axis=1); sums = np.where(ok, x, 0).sum(axis=1)
    out = np.full((count, 2), np.nan)
    np.divide(sums, n, out=out, where=n >= 4)
    return out[:, 0], out[:, 1], n / 8


def carry_uc(u):
    """TOCO gaps of up to 5 bins (10 s) carry the last value forward. Returns (filled, imputed mask)."""
    filled = u.copy()
    imputed = np.zeros(len(u), bool)
    last = np.nan
    age = 6
    for i, val in enumerate(u):
        if np.isfinite(val):
            last, age = val, 0
        else:
            age += 1
            if age <= 5 and np.isfinite(last):
                filled[i] = last
                imputed[i] = True
    return filled, imputed


def valid_rows(y, u):
    """Rows with FHR and no TOCO gap (after the 10-s carry) in the previous MEMORY bins; the first MEMORY bins are
    never used. Returns (rows, filled TOCO, imputed mask, history mask)."""
    filled, imputed = carry_uc(u)
    bad = (~np.isfinite(filled)).astype(int)
    invalid_history = np.convolve(bad, np.ones(MEMORY + 1, int), mode='full')[:len(u)] > 0
    history = ~invalid_history
    history[:MEMORY] = False
    return np.isfinite(y) & history, filled, imputed, history


# ---- cleaning ---------------------------------------------------------------------------------------------------
JUMP = 25.0; NEWLEVEL_S = 30; RESET_S = 10; INTERP_S = 15; BASE_S = 240


def clean_fhr(f):
    """Raw 4 Hz FHR artefact removal, before binning:
    1 out of range (< 50 or > 210 bpm) -> missing; 2 jumps > JUMP bpm from the last accepted value are removed
    unless they last > NEWLEVEL_S s and are stable (then a new level); a gap > RESET_S s resets the anchor;
    3 halving / doubling (ratio 0.4-0.6 or 1.7-2.3 to the centred 5-min median); 4 gaps <= INTERP_S s interpolated."""
    x = np.where((f >= 50) & (f <= 210), f, np.nan)
    out = np.full_like(x, np.nan); anchor = np.nan; gap = 0; pend = []

    def flush_as_level():
        nonlocal anchor
        for j in pend:
            out[j] = x[j]
        anchor = x[pend[-1]]; pend.clear()
    for i, v in enumerate(x):
        if not np.isfinite(v):
            gap += 1
            if gap > RESET_S * FS:
                anchor = np.nan; pend.clear()
            continue
        gap = 0
        if not np.isfinite(anchor) or abs(v - anchor) <= JUMP:
            pend.clear(); out[i] = v; anchor = v                       # back near the anchor: pending run was an artefact
        else:
            pend.append(i)
            seg = x[pend]
            if len(pend) > NEWLEVEL_S * FS and np.nanmax(seg[-NEWLEVEL_S * FS:]) - np.nanmin(seg[-NEWLEVEL_S * FS:]) <= 2 * JUMP:
                flush_as_level()                                      # sustained, stable: a real level change
    med = pd.Series(out).rolling(300 * FS, center=True, min_periods=30 * FS).median().to_numpy()
    r = out / med
    out[((r > .4) & (r < .6)) | ((r > 1.7) & (r < 2.3))] = np.nan
    s = pd.Series(out)
    filled = s.interpolate(limit=INTERP_S * FS, limit_area='inside')
    isn = s.isna().to_numpy()                                      # interpolate() also fills the start of long gaps; undo
    gl = np.zeros(len(s), int)
    i = 0
    while i < len(s):
        if isn[i]:
            j = i
            while j < len(s) and isn[j]:
                j += 1
            gl[i:j] = j - i; i = j
        else:
            i += 1
    out2 = filled.to_numpy(); out2[isn & (gl > INTERP_S * FS)] = np.nan
    return out2


def toco_tone(uc_bins, dt=DT):
    """Moving resting TOCO tone: centred rolling 10th percentile over BASE_S s (follows belt shifts)."""
    w = BASE_S // dt
    return pd.Series(uc_bins).rolling(w, center=True, min_periods=w // 4).quantile(.1).bfill().ffill().to_numpy()


def bin_fhr(f, raw, uc_zero_valid=False):
    """2-s bins as bin_signals, with the cleaned FHR `f` in place of the raw one. Returns (fhr_bins, uc_bins)."""
    r2 = raw.copy(); r2[:, 0] = np.where(np.isfinite(f), f, np.nan)
    y, u, _ = bin_signals(r2, uc_zero_valid)
    return y, u


# =================================================================================================================
# 3. input: TOCO -> uterine pressure p_ut(t)
# =================================================================================================================
# Default ("clean pulses"): TOCO has no absolute pressure scale, so we use only WHEN contractions happen.
PEAK = dict(medfilt=5, prominence=15, distance=40, width=8)     # in 2-s bins: 10 s, 15 TOCO units, 80 s, 16 s
SIGMA = 15.0                                                    # Gaussian bell sd, s (about 60 s wide at the base)


def time_grid(nb=NB):
    """Input grid every DU s; a 2-s bin's value is known at the bin's end (causal)."""
    return np.arange(int(nb * DT / DU) + 1) * DU


def detect_peaks(uc_bins):
    """Contraction peaks: find_peaks on medfilt(carried TOCO, 5 bins) with PEAK thresholds; missing bins after the
    10-s carry are set to the window's median TOCO. Returns peak times in s (bin ends)."""
    filled, _ = carry_uc(uc_bins)
    sm = medfilt(np.nan_to_num(filled, nan=np.nanmedian(uc_bins)), PEAK['medfilt'])
    pk, _ = find_peaks(sm, prominence=PEAK['prominence'], distance=PEAK['distance'], width=PEAK['width'])
    return (pk + 1) * DT


def pulses(peaks_s, t_grid, height, sigma=SIGMA, tone=TONE, rel=None):
    """Clean-pulse input: p_ut(t) = tone + height * sum_k r_k exp(-(t - t_k)^2 / (2 sigma^2)).
    Stands in for A(8) (one sin^2 contraction, p_rest 15 + p_con 70, 60 s) repeated at every TOCO peak t_k.
    r_k = 1 unless `rel` (relative strengths) is given."""
    p = np.full_like(t_grid, tone, dtype=float)
    rel = np.ones(len(peaks_s)) if rel is None else rel
    for tp, r in zip(peaks_s, rel):
        p += height * r * np.exp(-0.5 * ((t_grid - tp) / sigma) ** 2)
    return p


def linear_input(u, train=TRAIN):
    """Earlier mapping (first CTU fit): p_ut = 15 + a max(0, TOCO - rest), rest = 10th percentile of the training
    TOCO. Returns (excess on the input grid, rest); p_ut = TONE + a * excess."""
    nb = len(u)
    filled, _ = carry_uc(u)
    obs = u[:train][np.isfinite(u[:train])]; rest = float(np.percentile(obs, 10))
    x = np.where(np.isfinite(filled), filled, rest)
    excess = np.maximum(0.0, x - rest)
    return np.interp(time_grid(nb), (np.arange(nb) + 1) * DT, excess), rest


def sin2_contractions(uc_bins, dt=DT):
    """Earlier mapping (audit): contractions on the excess over the moving tone, with onset / offset
    at 20% of peak, relative strength r_k (0.5-1.5) and saturation flag. Returns (list of dicts, tone)."""
    filled, _ = carry_uc(uc_bins)
    tone = toco_tone(filled, dt)
    x = np.where(np.isfinite(filled), filled, tone)
    sm = medfilt(np.maximum(0, x - tone), 5)
    pk, _ = find_peaks(sm, prominence=max(8.0, .25 * np.percentile(sm, 95)), distance=30, width=5)
    out = []
    xmax = np.nanmax(x)
    for k, p in enumerate(pk):
        lo = pk[k - 1] + int(np.argmin(sm[pk[k - 1]:p])) if k else max(0, p - 60)
        hi = p + int(np.argmin(sm[p:pk[k + 1]])) if k + 1 < len(pk) else min(len(sm) - 1, p + 60)
        lo, hi = max(lo, p - 60), min(hi, p + 60)
        thr = .2 * sm[p]
        on = p - int(np.argmax(sm[lo:p + 1][::-1] < thr)) if np.any(sm[lo:p + 1] < thr) else lo
        off = p + int(np.argmax(sm[p:hi + 1] < thr)) if np.any(sm[p:hi + 1] < thr) else hi
        sat = xmax >= 95 and np.sum(x[max(0, p - 5):p + 6] >= xmax - .5) >= 2
        out.append(dict(peak=int(p), on=int(on), off=int(off), excess=float(sm[p]), saturated=bool(sat)))
    if out:
        med = np.median([c['excess'] for c in out])
        for c in out:
            r = float(np.clip(c['excess'] / med, .5, 1.5))
            c['r'] = max(r, 1.0) if c['saturated'] else r
    return out, tone


def sin2_pressure(cons, t_grid, A, dt=DT):
    """A(8) per detected contraction: p_ut = 15 + A r_k sin^2(pi (t - on_k) / D_k), D_k >= 30 s."""
    p = np.full(len(t_grid), TONE)
    for c in cons:
        t0 = (c['on'] + 1) * dt; D = max(30.0, (c['off'] - c['on']) * dt)
        m = (t_grid >= t0) & (t_grid <= t0 + D)
        p[m] += A * c['r'] * np.sin(np.pi * (t_grid[m] - t0) / D) ** 2
    return p


def shift(put, lag_s):
    """Shift an input-grid series by lag_s seconds (positive = later), linear interpolation, edges held."""
    n = np.arange(len(put)); return np.interp(n - lag_s / DU, n, put)


def to_grid(bins):
    """2-s bin series -> input grid (value known at each bin's end)."""
    nb = len(bins)
    return np.interp(time_grid(nb), (np.arange(nb) + 1) * DT, bins)


# =================================================================================================================
# 4. fit
# =================================================================================================================
A_BOUNDS = (0.25, 1.6)          # a: pulse height / 70 mmHg (17.5-112 mmHg above tone; covers B's 30 and A/B's 110)
W_BOUNDS = (0.0, 2.0)           # w_um, B(2) (B uses 1)
G_BOUNDS = (0.25, 4.0)          # g, common effector-gain multiplier (published = 1)
GRID = dict(a=(0.5, 1.0, 1.4), w_um=(0.0, 0.5, 1.0, 1.5, 2.0), g=(0.5, 1.0, 2.0))
NM = dict(maxfev=150, xatol=1e-3, fatol=1e-3)


def observe(f, y, rows):
    """Observation model: FHR_hat = b + (model FHR - model FHR at window start); b closed form = mean residual on
    the fit rows (least squares). Returns (training MSE, prediction, b)."""
    dev = f - f[0]; b = float(np.mean(y[rows] - dev[rows])); p = b + dev
    return float(np.mean((y[rows] - p[rows]) ** 2)), p, b


def minimise(loss, grid, bounds, nm=NM):
    """Grid search, then bounded Nelder-Mead from the best grid point; keep the grid point if NM does not improve.
    Returns (theta, n_evaluations)."""
    vals = [loss(np.array(t, float)) for t in grid]
    x0 = np.array(grid[int(np.argmin(vals))], float)
    r = minimize(loss, x0, method='Nelder-Mead', bounds=bounds, options=dict(nm))
    return (r.x if r.fun <= min(vals) else x0), len(grid) + r.nfev


def fit_pulses(y, rows, peaks_s, rel=None):
    """Fit a, w_um, g (theta = log a, w_um, log g) and b to FHR bins y on `rows` with the clean-pulse input."""
    tg = time_grid(len(y))
    put = lambda th: pulses(peaks_s, tg, np.exp(th[0]) * P_CON, rel=rel)

    def loss(th):
        f = simulate(put(th), th[1], np.exp(th[2]), len(y))
        return np.inf if f is None else observe(f, y, rows)[0]
    grid = [(np.log(a), w, np.log(g)) for a in GRID['a'] for w in GRID['w_um'] for g in GRID['g']]
    bounds = [np.log(A_BOUNDS), W_BOUNDS, np.log(G_BOUNDS)]
    th, nev = minimise(loss, grid, bounds)
    p_ut = put(th); f = simulate(p_ut, th[1], np.exp(th[2]), len(y))
    mse, pred, b = observe(f, y, rows)
    return dict(a=float(np.exp(th[0])), w_um=float(th[1]), g=float(np.exp(th[2])), b=b, train_mse=mse,
                evaluations=nev, fhr_model=pred, p_ut=p_ut, t_input=tg)


def fit(fhr_bins, toco_bins, train=TRAIN):
    """ENTRY POINT. Fit the model to one recording window given as 2-s bins (NaN = missing), using the first
    `train` bins: TOCO peaks -> clean pulses -> a / w_um / g / b. Rows = FHR present and no TOCO gap in the
    previous 120 s. Returns fitted parameters, the model FHR over the whole window and the pressure input."""
    y = np.asarray(fhr_bins, float); u = np.asarray(toco_bins, float)
    rows, _, _, _ = valid_rows(y, u)
    tr = np.zeros(len(y), bool); tr[np.flatnonzero(rows[:train])] = True
    pk = detect_peaks(u)
    res = fit_pulses(y, tr, pk)
    res.update(peaks_s=pk, rows=rows, train_rows=tr)
    return res


# =================================================================================================================
# 5. parameter table
# =================================================================================================================
CLASS = {1: 'fixed at the paper\'s value', 2: 'fixed by our choice', 3: 'fitted per 20-min window'}
_TAB = {'A-T1a': 'A Table 1a', 'A-T1b': 'A Table 1b', 'B-T1': 'B Table 1', 'A-T2': 'A Table 2', 'A-T3': 'A Table 3'}
H = 'mmHg·s/ml'; CO = 'ml/mmHg'; EL = 'mmHg/ml'
# name: (symbol, unit, meaning, part, equations)
META = {
    'F_TREF': ('T_ref', 's', 'fetal reference heart period (60/135 s = 135 bpm)', 'heart', 'A(25)'),
    'F_VART0': ('V_art,0', 'ml', 'fetal arterial unstressed volume', 'heart', 'A(1)'),
    'F_VVEN0': ('V_ven,0,ref', 'ml', 'fetal venous unstressed volume (reference; sympathetically regulated)', 'heart', 'A(1), A(24)'),
    'F_VTIS0': ('V_tis,0', 'ml', 'fetal body-tissue unstressed volume', 'heart', 'A(1)'),
    'F_VCER0': ('V_cer,0', 'ml', 'fetal brain unstressed volume', 'heart', 'A(1)'),
    'F_VVIL0': ('V_villi,0', 'ml', 'placental villi unstressed volume (fetal side of the placenta)', 'placenta', 'A(1)'),
    'F_CART': ('C_art', CO, 'fetal arterial compliance', 'heart', 'A(1)'),
    'F_CVEN': ('C_ven', CO, 'fetal venous compliance', 'heart', 'A(1)'),
    'F_CTIS': ('C_tis', CO, 'fetal body-tissue compliance', 'heart', 'A(1)'),
    'F_CCER': ('C_cer', CO, 'fetal brain compliance', 'heart', 'A(1)'),
    'F_CVIL': ('C_villi', CO, 'villi compliance', 'placenta', 'A(1)'),
    'F_RTISA': ('R_tisa,ref', H, 'fetal body arterial resistance (reference; sympathetically regulated)', 'heart', 'A(5), A(24)'),
    'F_RTISV': ('R_tisv', H, 'fetal body venous resistance', 'heart', 'A(5)'),
    'F_RCERA': ('R_cera,ref', H, 'fetal brain arterial resistance (reference; autoregulated)', 'heart', 'A(18)-(19)'),
    'F_RCERV': ('R_cerv', H, 'fetal brain venous resistance', 'heart', 'A(5)'),
    'F_RUMA': ('R_uma,ref', H, 'umbilical artery resistance at rest', 'cord', 'A(6), B(1)'),
    'F_RUMV': ('R_umv,ref', H, 'umbilical vein resistance at rest', 'cord', 'A(6), B(1)'),
    'F_EMAX': ('E_max,ref', EL, 'fetal ventricle end-ejection elastance (reference; regulated)', 'heart', 'A(9), A(24)'),
    'F_EMIN': ('E_min', EL, 'fetal ventricle end-diastolic elastance', 'heart', 'A(9)'),
    'F_V0MAX': ('V_max,0', 'ml', 'fetal ventricle unstressed volume at end-ejection', 'heart', 'A(9)'),
    'F_V0MIN': ('V_min,0', 'ml', 'fetal ventricle unstressed volume at end-diastole', 'heart', 'A(9)'),
    'M_T': ('T_m', 's', 'maternal heart period (60/80 s = 80 bpm, not regulated)', 'mother', 'A(10)'),
    'M_VART0': ('V_art,0', 'ml', 'maternal arterial unstressed volume', 'mother', 'A(1)'),
    'M_VVEN0': ('V_ven,0', 'ml', 'maternal venous unstressed volume', 'mother', 'A(1)'),
    'M_VTIS0': ('V_tis,0', 'ml', 'maternal tissue unstressed volume', 'mother', 'A(1)'),
    'M_CART': ('C_art', CO, 'maternal arterial compliance', 'mother', 'A(1)'),
    'M_CVEN': ('C_ven', CO, 'maternal venous compliance', 'mother', 'A(1)'),
    'M_CTIS': ('C_tis', CO, 'maternal tissue compliance', 'mother', 'A(1)'),
    'M_RTISA': ('R_tisa', H, 'maternal tissue arterial resistance', 'mother', 'A(5)'),
    'M_RTISV': ('R_tisv', H, 'maternal tissue venous resistance', 'mother', 'A(5)'),
    'M_RUTA': ('R_uta', H, 'uterine (spiral) artery resistance; NOT squeezed by contractions', 'spiral', 'A(5)'),
    'M_RUTV': ('R_utv,ref', H, 'uterine vein resistance at rest', 'uterus', 'A(6)-(7)'),
    'M_R0UTV': ('R_0,utv', H, 'uterine vein collapse law: base resistance', 'uterus', 'A(6)'),
    'M_P0UTV': ('p_0,utv', 'mmHg', 'uterine vein collapse law: midpoint transmural pressure', 'uterus', 'A(6)'),
    'M_P1UTV': ('p_1,utv', 'mmHg', 'uterine vein collapse law: width', 'uterus', 'A(6)'),
    'M_EMAX': ('E_max', EL, 'maternal ventricle end-ejection elastance', 'mother', 'A(9)'),
    'M_EMIN': ('E_min', EL, 'maternal ventricle end-diastolic elastance', 'mother', 'A(9)'),
    'M_V0MAX': ('V_max,0', 'ml', 'maternal ventricle unstressed volume at end-ejection', 'mother', 'A(9)'),
    'M_V0MIN': ('V_min,0', 'ml', 'maternal ventricle unstressed volume at end-diastole', 'mother', 'A(9)'),
    'IVS_VMAX': ('V_max,ivs', 'ml', 'intervillous space: maximal blood volume (2 x 175 ml, B p.2475 rule)', 'placenta', 'A(3)'),
    'IVS_P0': ('p_0,ivs', 'mmHg', 'intervillous space: transmural pressure at rest (B p.2475 rule)', 'placenta', 'A(3)'),
    'IVS_P1': ('p_1,ivs', 'mmHg', 'intervillous space: stiffness width, V_max / (pi C_ivs,max) (B p.2475 rule)', 'placenta', 'A(3)'),
    'U_R0UMV': ('R_0,umv', H, 'umbilical vein collapse law: base resistance', 'cord', 'A(6), B(1)'),
    'U_P0UMV': ('p_0,umv', 'mmHg', 'umbilical vein collapse law: midpoint', 'cord', 'A(6), B(1)'),
    'U_P1UMV': ('p_1,umv', 'mmHg', 'umbilical vein collapse law: width', 'cord', 'A(6), B(1)'),
    'U_R0UMA': ('R_0,uma', H, 'umbilical artery collapse law: base resistance', 'cord', 'A(6), B(1)'),
    'U_P0UMA': ('p_0,uma', 'mmHg', 'umbilical artery collapse law: midpoint', 'cord', 'A(6), B(1)'),
    'U_P1UMA': ('p_1,uma', 'mmHg', 'umbilical artery collapse law: width', 'cord', 'A(6), B(1)'),
    'O_DPLAC': ('D_plac', 'ml O2/(s·mmHg)', 'placental O2 diffusion capacity, intervillous blood -> villi', 'placenta', 'A(16)'),
    'O_DLUNG': ('D_lung', 'ml O2/(s·mmHg)', 'maternal lung O2 uptake capacity', 'mother', 'A(16)'),
    'O_QTISF': ('Q_M,0,tis,f', 'ml O2/s', 'fetal body O2 consumption', 'heart', 'A(17)'),
    'O_QCERF': ('Q_M,0,cer,f', 'ml O2/s', 'fetal brain O2 consumption', 'heart', 'A(17)'),
    'O_QTISM': ('Q_M,0,tis,m', 'ml O2/s', 'maternal O2 consumption', 'mother', 'A(17)'),
    'O_THCER': ('pO2,th,cer,f', 'mmHg', 'fetal brain: pO2 below which consumption falls', 'heart', 'A(17)'),
    'O_THTIS': ('pO2,th,tis,f', 'mmHg', 'fetal body: pO2 below which consumption falls', 'heart', 'A(17)'),
    'O_THTISM': ('pO2,th,tis,m', 'mmHg', 'maternal threshold: not stated, "irrelevant" (B p.2475); we use constant uptake', 'mother', 'A(17)'),
    'O_PAIR': ('pO2,air', 'mmHg', 'inspired-air pO2', 'mother', 'A(16)'),
    'O_ALPHA': ('alpha', 'ml O2/g Hb', 'O2 binding capacity of haemoglobin', 'placenta', 'A(12)'),
    'O_BETA': ('beta', 'ml O2/(ml·mmHg)', 'dissolved O2 solubility', 'placenta', 'A(12)'),
    'O_HBF': ('Hb_f', 'g/ml', 'fetal haemoglobin concentration', 'heart', 'A(12)'),
    'O_HBM': ('Hb_m', 'g/ml', 'maternal haemoglobin concentration', 'mother', 'A(12)'),
    'O_C1F': ('c1,f', 'mmHg³', 'fetal O2 dissociation curve constant (left-shifted)', 'heart', 'A(13)'),
    'O_C1M': ('c1,m', 'mmHg³', 'maternal O2 dissociation curve constant', 'mother', 'A(13)'),
    'O_C2F': ('c2,f', 'mmHg²', 'fetal O2 dissociation curve constant', 'heart', 'A(13)'),
    'O_C2M': ('c2,m', 'mmHg²', 'maternal O2 dissociation curve constant', 'mother', 'A(13)'),
    'G_CREF': ('cO2,a,ref', 'ml O2/ml', 'brain autoregulation: reference arterial O2 content', 'heart', 'A(19)'),
    'G_GAMMA': ('gamma', '–', 'brain autoregulation strength', 'heart', 'A(19)'),
    'G_TAURC': ('tau_Rcera', 's', 'brain autoregulation time constant', 'heart', 'A(18)'),
    'G_ETA1B': ('eta1,b', '–', 'baroreceptor curve offset, ln 0.5 exact from B(3)', 'baro', 'A(21), B(3)'),
    'G_ETA2B': ('eta2,b', '1/mmHg', 'baroreceptor curve slope', 'baro', 'A(21), B(4)'),
    'G_RBMIN': ('r_b,min', '–', 'baroreceptor output minimum', 'baro', 'A(21)'),
    'G_RBMAX': ('r_b,max', '–', 'baroreceptor output maximum', 'baro', 'A(21)'),
    'G_TAURB': ('tau_rb', 's', 'baroreceptor low-pass time constant', 'baro', 'A(20)'),
    'G_YREFB': ('y_ref,b', 'mmHg', 'baroreceptor set point: fetal arterial transmural pressure', 'baro', 'A(21)'),
    'G_ETA1C': ('eta1,c', '–', 'chemoreceptor curve offset, ln(0.001/1.001) exact from B(3)', 'chemo', 'A(21), B(3)'),
    'G_ETA2C': ('eta2,c', '1/mmHg', 'chemoreceptor curve slope', 'chemo', 'A(21), B(4)'),
    'G_RCMIN': ('r_c,min', '–', 'chemoreceptor output minimum', 'chemo', 'A(21)'),
    'G_RCMAX': ('r_c,max', '–', 'chemoreceptor output maximum', 'chemo', 'A(21)'),
    'G_TAURCH': ('tau_rc', 's', 'chemoreceptor low-pass time constant', 'chemo', 'A(20)'),
    'G_YREFC': ('y_ref,c', 'mmHg', 'chemoreceptor set point: fetal arterial pO2', 'chemo', 'A(21)'),
    'G_WBS': ('w_b,s', '–', 'sympathetic weight of the baroreceptor (inhibitory)', 'baro', 'A(22)'),
    'G_WCS': ('w_c,s', '–', 'sympathetic weight of the chemoreceptor', 'chemo', 'A(22)'),
    'G_WBV': ('w_b,v', '–', 'vagal weight of the baroreceptor', 'baro', 'A(22)'),
    'G_WCV': ('w_c,v', '–', 'vagal weight of the chemoreceptor', 'chemo', 'A(22)'),
    'G_DTV': ('D_Tv', 's', 'vagal heart-period effector: delay', 'heart', 'A(23), A(25)'),
    'G_TAUTV': ('tau_Tv', 's', 'vagal heart-period effector: time constant', 'heart', 'A(23), A(25)'),
    'G_KTV': ('k_Tv', '–', 'vagal heart-period effector: gain (x g)', 'heart', 'A(23), A(25)'),
    'G_DTS': ('D_Ts', 's', 'sympathetic heart-period effector: delay', 'heart', 'A(23), A(25)'),
    'G_TAUTS': ('tau_Ts', 's', 'sympathetic heart-period effector: time constant', 'heart', 'A(23), A(25)'),
    'G_KTS': ('k_Ts', '–', 'sympathetic heart-period effector: gain (x g)', 'heart', 'A(23), A(25)'),
    'G_DE': ('D_Emax', 's', 'contractility effector: delay', 'heart', 'A(23)-(24)'),
    'G_TAUE': ('tau_Emax', 's', 'contractility effector: time constant', 'heart', 'A(23)-(24)'),
    'G_KE': ('k_Emax', '–', 'contractility effector: gain (x g)', 'heart', 'A(23)-(24)'),
    'G_DR': ('D_Rtisa', 's', 'peripheral-resistance effector: delay', 'heart', 'A(23)-(24)'),
    'G_TAUR': ('tau_Rtisa', 's', 'peripheral-resistance effector: time constant', 'heart', 'A(23)-(24)'),
    'G_KR': ('k_Rtisa', '–', 'peripheral-resistance effector: gain (x g)', 'heart', 'A(23)-(24)'),
    'G_DV': ('D_Vven,0', 's', 'venous-volume effector: delay', 'heart', 'A(23)-(24)'),
    'G_TAUV': ('tau_Vven,0', 's', 'venous-volume effector: time constant', 'heart', 'A(23)-(24)'),
    'G_KV': ('k_Vven,0', '–', 'venous-volume effector: gain (x g)', 'heart', 'A(23)-(24)'),
    'S_PREST': ('p_rest', 'mmHg', 'resting uterine pressure in the model', 'uterus', 'A(8), B(2)'),
    'S_MRT': ('m_Rtisa', '–', 'thesis noise on R_tisa: off (no variability in fitting)', 'heart', 'Th Eq. 4.1-4.2'),
    'S_MRB': ('m_rb', '–', 'thesis noise on baroreceptor output: off', 'baro', 'Th Eq. 4.1-4.2'),
    'S_MEV': ('m_ev', '–', 'thesis noise on vagal efferent: off', 'heart', 'Th Eq. 4.1-4.2'),
}
OURS = {'O_THTISM', 'S_MRT', 'S_MRB', 'S_MEV'}


def _ref(src):
    s = src.split(':')[0].split(' (')[0].strip()
    return _TAB.get(s, s.replace('A s', 'A §').replace('B s', 'B §'))


def parameter_table():
    """Every parameter of the model and of our data pipeline, as a list of dicts:
    name, symbol, value or range, unit, meaning, part (where it sits in the body), equations, class (1/2/3),
    source (Jongen table/section) and note. Built from jongen.PARAMS and the constants of this module."""
    rows = []
    for name, val, src in J.PARAMS:
        if name in ('S_WUM', 'S_GAIN'):
            continue
        sym, unit, meaning, part, eq = META[name]
        cls = 2 if name in OURS else 1
        note = ''
        if src.startswith('DERIVED'):
            note = 'Not printed in any table; derived from the paper\'s stated rule (B p.2475).'
        if name == 'O_DPLAC':
            note = ('Set by the paper\'s rule (steady-state fetal arterial pO2 = 18.5 mmHg, B §2.2.2). The table prints '
                    '0.082, which gives 18.25 mmHg in our implementation.')
        if name.startswith('G_ETA1'):
            note = 'Exact value from B(3) so the receptor output is 0 at its set point; the table prints it rounded.'
        if name == 'O_THTISM':
            note = 'Not stated in the papers.'
        rows.append(dict(name=name, symbol=sym, value=float(val), unit=unit, meaning=meaning, part=part, equations=eq,
                         cls=cls, source=_ref(src) if cls == 1 else 'ours', note=note))
    ours = [
        ('tone', 'p_tone', TONE, 'mmHg', 'pressure between pulses in our TOCO-derived input; set to the paper\'s p_rest because TOCO has no absolute scale', 'toco', 'A(8)'),
        ('p_con_ref', 'p_con', P_CON, 'mmHg', 'reference pulse height above the tone; the fitted a multiplies it (a = 1 is the paper\'s reference contraction, peak 85 mmHg)', 'toco', 'A(8)'),
        ('sigma', 'sigma', SIGMA, 's', 'width of each Gaussian pulse (about 60 s from start to end, like the paper\'s 60 s sin² contraction)', 'toco', 'replaces A(8)'),
        ('peak_medfilt', '', PEAK['medfilt'] * DT, 's', 'TOCO median filter before peak detection', 'toco', ''),
        ('peak_prominence', '', PEAK['prominence'], 'TOCO units', 'minimum peak prominence to count as a contraction', 'toco', ''),
        ('peak_distance', '', PEAK['distance'] * DT, 's', 'minimum time between two contraction peaks', 'toco', ''),
        ('peak_width', '', PEAK['width'] * DT, 's', 'minimum peak width', 'toco', ''),
        ('uc_carry', '', 10, 's', 'TOCO gaps up to this length are carried forward', 'toco', ''),
        ('uc_flat', '', FLAT_S, 's', 'TOCO = 0 is a valid reading, except an exactly constant value for at least this long (sensor off)', 'toco', ''),
        ('fhr_range', '', '50–210', 'bpm', 'FHR samples outside are missing', 'heart', ''),
        ('fhr_jump', '', JUMP, 'bpm', 'FHR cleaning: jumps larger than this from the last accepted value are artefacts unless sustained > 30 s', 'heart', ''),
        ('fhr_interp', '', INTERP_S, 's', 'FHR gaps up to this length are interpolated', 'heart', ''),
        ('window', '', '20 + 2 + 10', 'min', 'fit on the first 20 min, skip 2 min, compare on the next 10 min', 'toco', ''),
        ('uc_history', '', MEMORY * DT, 's', 'a row is used only if TOCO has no gap in the previous 120 s', 'toco', ''),
        ('initial_state', '', 'steady state', '', 'every window starts from the model\'s resting steady state (A §2.4 criterion)', 'heart', ''),
        ('dt_input', '', DU, 's', 'input grid for p_ut; the model itself steps at 0.02 s (A §2.4)', 'toco', ''),
        ('observation', '', '2-s mean', '', 'model FHR averaged over each 2-s CTG bin', 'heart', 'A(25)'),
    ]
    for name, sym, val, unit, meaning, part, eq in ours:
        rows.append(dict(name=name, symbol=sym, value=val, unit=unit, meaning=meaning, part=part, equations=eq, cls=2,
                         source='ours', note=''))
    method = 'grid, then bounded Nelder-Mead (≤ 150 evaluations) on the training MSE'
    rows += [
        dict(name='a', symbol='a', value=None, range=list(A_BOUNDS), grid=list(GRID['a']), unit='× 70 mmHg',
             meaning='contraction strength: pulse height = a × 70 mmHg above the 15 mmHg tone', part='uterus',
             equations='A(8)', cls=3, source='ours', method=method + ' (log scale)', note=''),
        dict(name='S_WUM', symbol='w_um', value=None, range=list(W_BOUNDS), grid=list(GRID['w_um']), unit='–',
             meaning='how much of the uterine pressure squeezes the cord: p_um = w_um (p_ut − p_rest). 0 = Part A, 1 = Part B',
             part='cord', equations='B(1)-(2)', cls=3, source='B(2), B §2.3 (paper uses 0 or 1)', method=method, note=''),
        dict(name='S_GAIN', symbol='g', value=None, range=list(G_BOUNDS), grid=list(GRID['g']), unit='–',
             meaning='common multiplier of all five effector gains k; g = 1 is the published regulation', part='heart',
             equations='A(23)-(25)', cls=3, source='ours (not in the paper)', method=method + ' (log scale)', note=''),
        dict(name='b', symbol='b', value=None, range=None, unit='bpm',
             meaning='baseline FHR: model FHR is shifted so its mean residual on the training rows is 0', part='heart',
             equations='', cls=3, source='ours', method='closed form (least squares)', note=''),
    ]
    return rows


def export(path):
    """Write the parameter table to `path` (JSON). Returns the counts per class."""
    rows = parameter_table()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    counts = {str(k): sum(r['cls'] == k for r in rows) for k in CLASS}
    out = dict(generated_by='jongen2016/model.py parameter_table()', classes=CLASS, counts=counts, parameters=rows)
    path.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    return counts


if __name__ == '__main__':
    print(export(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'docs/parameters.json'))
