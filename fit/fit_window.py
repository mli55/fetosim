"""Fit the model to one CTG window: a (contraction strength), w_um (cord-compression weight, Jongen 2016 B(2)) and g
(common multiplier of the effector gains of A(23)) by a grid search followed by bounded Nelder-Mead on the training
MSE; the baseline FHR b in closed form. Also assembles the full parameter table: the model's parameter set plus the
input and fitting choices, each with its class (1 fixed at the paper's value, 2 fixed by our choice, 3 fitted).

    from fetosim.fit import fit_window as FW
    res = FW.fit(fhr_bins, toco_bins)        # dict: a, w_um, g, b, fhr_model, p_ut, peaks_s, ...
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

if __name__ == '__main__':                  # run as a script: put the folder holding fetosim/ on sys.path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from fetosim.input.ctg_signals import DT, FLAT_S, INTERP_S, JUMP, MEMORY, TRAIN, valid_rows
from fetosim.input.toco_to_pressure import PEAK, SIGMA, detect_peaks, pulses, time_grid
from fetosim.model.simulate import DU, ENTRIES, P_CON, TONE, simulate as _simulate


def simulate(put, w_um, g, nb):
    """Model FHR on 2-s CTG bins (see model/simulate.simulate)."""
    return _simulate(put, w_um, g, nb, DT)


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


CLASS = {1: 'fixed at the paper\'s value', 2: 'fixed by our choice', 3: 'fitted per 20-min window'}
_CLS = {'paper': 1, 'our choice': 2, 'fitted': 3}


def parameter_table():
    """Every parameter of the model and of the input / fitting pipeline, as a list of dicts: name, symbol, value or
    range, unit, meaning, part (where it sits in the body), equations, cls (1/2/3), source and note. Model entries come
    from the parameter set (model/parameters/jongen2016.json); input and fitting choices from the constants here."""
    rows = []
    for e in ENTRIES:
        if e['class'] == 'fitted':
            continue
        rows.append(dict(name=e['name'], symbol=e['symbol'], value=float(e['value']), unit=e['unit'], meaning=e['meaning'],
                         part=e['part'], equations=e['equations'], cls=_CLS[e['class']], source=e['source'], note=e['note']))
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
    fitted = {e['name']: e for e in ENTRIES if e['class'] == 'fitted'}

    def model_row(name, bounds, grid, how):
        e = fitted[name]
        return dict(name=name, symbol=e['symbol'], value=None, range=list(bounds), grid=list(grid), unit=e['unit'],
                    meaning=e['meaning'], part=e['part'], equations=e['equations'], cls=3, source=e['source'],
                    method=how, note='')
    rows += [
        dict(name='a', symbol='a', value=None, range=list(A_BOUNDS), grid=list(GRID['a']), unit='× 70 mmHg',
             meaning='contraction strength: pulse height = a × 70 mmHg above the 15 mmHg tone', part='uterus',
             equations='A(8)', cls=3, source='ours', method=method + ' (log scale)', note=''),
        model_row('S_WUM', W_BOUNDS, GRID['w_um'], method),
        model_row('S_GAIN', G_BOUNDS, GRID['g'], method + ' (log scale)'),
        dict(name='b', symbol='b', value=None, range=None, unit='bpm',
             meaning='baseline FHR: model FHR is shifted so its mean residual on the training rows is 0', part='heart',
             equations='', cls=3, source='ours', method='closed form (least squares)', note=''),
    ]
    return rows


def export(path):
    """Write the full parameter table to `path` (JSON). Returns the counts per class."""
    rows = parameter_table()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    counts = {str(k): sum(r['cls'] == k for r in rows) for k in CLASS}
    out = dict(generated_by='fit/fit_window.py parameter_table()', classes=CLASS, counts=counts, parameters=rows)
    path.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    return counts


if __name__ == '__main__':
    print(export(sys.argv[1] if len(sys.argv) > 1 else 'parameter_table.json'))
