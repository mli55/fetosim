"""CTG signals: 4 Hz FHR / TOCO -> 2-s bins, the "TOCO = 0 is a valid reading" rule, row masks for fitting, FHR
artefact cleaning and the moving TOCO tone. Nothing here is from Jongen 2016; these are our signal-handling choices.

    from fetosim.input import ctg_signals as CS
    fhr_bins, toco_bins, _ = CS.bin_signals(raw_4hz, uc_zero_valid=True)
"""
import numpy as np
import pandas as pd

DT = 2                          # 2-s bins (8 samples at 4 Hz)
TRAIN, GAP, TEST = 600, 60, 300  # bins: 20 min fit, 2 min unused gap, 10 min comparison
NB = TRAIN + GAP + TEST
MEMORY = 60                     # bins: a row needs no TOCO gap in the previous 120 s
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


# ---- FHR cleaning and TOCO tone --------------------------------------------------------------------------------
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
