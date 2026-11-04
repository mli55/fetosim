"""TOCO -> uterine pressure p_ut(t), the model's only input.

Default ("clean pulses"): TOCO has no absolute pressure scale, so only the timing of contractions is used. Peaks of the
median-filtered TOCO each get one Gaussian pulse on a 15 mmHg resting tone; the pulse height is a x 70 mmHg, where
a = 1 is the reference contraction of Jongen 2016 A §2.4 / A(8) (peak 85 mmHg, about 60 s long).
Two earlier mappings are kept for comparison: linear_input (p_ut = 15 + a max(0, TOCO - rest)) and the per-contraction
sin^2 pulses of sin2_contractions / sin2_pressure.

    from fetosim.input import toco_to_pressure as TP
    p_ut = TP.pulses(TP.detect_peaks(toco_bins), TP.time_grid(len(toco_bins)), height=70.0)
"""
import numpy as np
from scipy.signal import find_peaks, medfilt

from fetosim.input.ctg_signals import DT, NB, TRAIN, carry_uc, toco_tone
from fetosim.model.simulate import DU, P_CON, TONE  # noqa: F401  (P_CON: reference height used with a)

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
