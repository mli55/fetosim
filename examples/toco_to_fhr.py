#!/usr/bin/env python3
"""Shortest path from TOCO to FHR: three functions, one per stage, built on the project's own modules,
plus a plot of every intermediate signal (one panel per block of the TOCO -> FHR diagram).

    toco_bins(toco_4hz)            4 Hz TOCO -> 2-s bins
    uterine_pressure(bins, a)      bins -> peaks -> p_ut(t) = 15 + 70a sum_k exp(-(t - t_k)^2 / 2 sigma^2)
    fhr(p_ut, nb, w_um, g)         p_ut -> Jongen 2016 ODEs (core.c) -> FHR = 60 / T, in 2-s bins

usage: python examples/toco_to_fhr.py [--csv FILE] [--a A] [--w-um W] [--g G] [--out PNG] [--no-show]

  --csv FILE   4 Hz CTG as CSV with a 'toco' column and optionally an 'fhr' column (overlaid on the model FHR).
               Without it, a synthetic 20-min TOCO with 4 contractions is used.
  --a A        contraction strength, pulse height = A x 70 mmHg (default 1)
  --w-um W     cord-compression weight (default 1; 0 = placenta route only)
  --g G        common reflex-gain multiplier (default 1)
  --out PNG    figure file (default outputs/toco_to_fhr.png, relative to the current folder)
  --no-show    save the figure without opening a window
"""
import sys

if '--help' in sys.argv or '-h' in sys.argv:
    print(__doc__)
    sys.exit(0)

import argparse
from pathlib import Path

ap = argparse.ArgumentParser(add_help=False)
ap.add_argument('--csv', type=Path)
ap.add_argument('--a', type=float, default=1.0)
ap.add_argument('--w-um', type=float, default=1.0)
ap.add_argument('--g', type=float, default=1.0)
ap.add_argument('--out', type=Path, default=Path('outputs/toco_to_fhr.png'))
ap.add_argument('--no-show', action='store_true')
ARGS = ap.parse_args() if __name__ == '__main__' else None
if ARGS and ARGS.csv and not ARGS.csv.is_file():
    sys.exit(f'error: {ARGS.csv} not found')

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # folder holding fetosim/; no install needed

from fetosim.input.toco_to_pressure import detect_peaks, pulses, time_grid  # noqa: E402
from fetosim.model import simulate as S  # noqa: E402


def toco_bins(toco_4hz: np.ndarray) -> np.ndarray:
    """Average 4 Hz TOCO into 2-s bins (8 samples each)."""
    n = len(toco_4hz) // 8
    return np.asarray(toco_4hz[:n * 8], float).reshape(n, 8).mean(1)


def uterine_pressure(bins: np.ndarray, a: float = 1.0) -> np.ndarray:
    """One Gaussian pulse of height a * 70 mmHg on a 15 mmHg tone at each TOCO peak, on the model's input grid."""
    return pulses(detect_peaks(bins), time_grid(len(bins)), a * S.P_CON)


def fhr(p_ut: np.ndarray, nb: int, w_um: float = 1.0, g: float = 1.0) -> np.ndarray:
    """Run the model on p_ut and return FHR (bpm) averaged over nb 2-s bins. w_um = cord compression, g = reflex gain."""
    f = S.simulate(p_ut, w_um, g, nb)
    if f is None:
        raise RuntimeError('model run became non-finite')
    return f


def synthetic_toco() -> np.ndarray:
    """20 min of 4 Hz TOCO: tone 10, four Gaussian contractions of height 60 at 3, 8, 13 and 18 min."""
    t = np.arange(20 * 60 * 4) / 4
    return 10 + sum(60 * np.exp(-0.5 * ((t - c) / 20) ** 2) for c in (180, 480, 780, 1080))


def plot_blocks(toco_4hz, fhr_4hz=None, a=1.0, w_um=1.0, g=1.0, out=None, show=True):
    """One panel per diagram block: TOCO, p_ut, placenta inflow / cord flow, fetal pO2, fetal pressure, receptors,
    vagal / sympathetic signals, heart-period effectors and FHR (model, plus measured if given)."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams['font.family'] = 'Helvetica'
    c = sns.color_palette('deep')

    bins = toco_bins(toco_4hz)
    peaks = detect_peaks(bins)
    put = uterine_pressure(bins, a)
    o, _, ok = S.run(put, w_um=w_um, g=g, du=S.DU)
    if not ok:
        raise RuntimeError('model run became non-finite')
    f_bins = fhr(put, len(bins), w_um, g)
    t4, tb, tm = np.arange(len(toco_4hz)) / 4 / 60, (np.arange(len(bins)) + 1) * 2 / 60, o['t'] / 60

    panels = [
        ('TOCO (input)', 'TOCO units', [(t4, toco_4hz, '4 Hz', '-'), (tb, bins, '2-s bins', '--')]),
        ('Uterine pressure', 'mmHg', [(time_grid(len(bins)) / 60, put, 'p_ut', '-')]),
        ('Placenta inflow / cord flow', 'ml/min', [(tm, 60 * o['q_uta'], 'q_uta (uterine artery)', '-'),
                                                   (tm, 60 * o['q_uma'], 'q_uma (umbilical artery)', '--')]),
        ('Fetal oxygen', 'mmHg', [(tm, o['po2_a'], 'pO2_a', '-')]),
        ('Fetal pressure', 'mmHg', [(tm, o['ptm_a'], 'p_a', '-')]),
        ('Receptors', '(-1 to 1)', [(tm, o['r_c'], 'r_c chemo', '-'), (tm, o['r_b'], 'r_b baro', '--')]),
        ('Vagal and sympathetic', '', [(tm, o['e_v'], 'e_v', '-'), (tm, o['e_s'], 'e_s', '--')]),
        ('Effectors on heart period', '', [(tm, o['z_tv'], 'z_Tv', '-'), (tm, o['z_ts'], 'z_Ts', '--')]),
        ('FHR (output)', 'bpm', [(tm, o['fhr'], 'model, 0.1 s', '-'), (tb, f_bins, 'model, 2-s bins', '--')]),
    ]
    if fhr_4hz is not None:
        fm = np.where((fhr_4hz >= 50) & (fhr_4hz <= 210), fhr_4hz, np.nan)
        panels[-1][2].append((t4, fm, 'measured', ':'))

    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 2 * len(panels)), sharex=True)
    for ax, (title, unit, lines) in zip(axes, panels):
        for i, (x, y, lab, ls) in enumerate(lines):
            ax.plot(x, y, ls, color=c[i], lw=1.4, label=lab)
        for p in peaks:
            ax.axvline(p / 60, color='0.6', lw=0.6, ls=':')
        ax.set_title(title, loc='left', fontsize=10, fontweight='bold')
        ax.set_ylabel(unit)
        ax.legend(loc='upper right', fontsize=8, frameon=False, ncol=len(lines))
        ax.spines[['top', 'right']].set_visible(False)
    axes[-1].set_xlabel('time (min); dotted lines = detected TOCO peaks')
    fig.suptitle(f'TOCO to FHR, block by block (a = {a}, w_um = {w_um}, g = {g})', x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150)
        print(f'saved {out}')
    if show:
        plt.show()


if __name__ == '__main__':
    if ARGS.csv:
        import pandas as pd
        df = pd.read_csv(ARGS.csv)
        if 'toco' not in df:
            sys.exit(f"error: {ARGS.csv} has no 'toco' column (columns: {', '.join(df.columns)})")
        toco, fhr_meas = df['toco'].to_numpy(float), (df['fhr'].to_numpy(float) if 'fhr' in df else None)
    else:
        toco, fhr_meas = synthetic_toco(), None
    plot_blocks(toco, fhr_meas, ARGS.a, ARGS.w_um, ARGS.g, ARGS.out, not ARGS.no_show)
