"""Randall et al. 2019 baroreflex + respiratory sinus arrhythmia (RSA) heart-rate model, ported to Python.

Source: Randall EB, Billeschou A, Brinth LS, Mehlsen J, Olufsen MS. A model-based analysis of autonomic nervous
function in response to the Valsalva maneuver. J Appl Physiol 127:1386-1402, 2019. doi:10.1152/japplphysiol.00015.2019
(PMC6879835). Authors' code: Zenodo doi:10.5281/zenodo.3234934 (MIT license), files driver_baroreflex.f (FCN),
load_global.m, initialconditions.m, spline.f, splint.f. The equations below follow that code line by line.

Open loop: inputs are systolic blood pressure SBP(t) (mmHg) and thoracic pressure Pth(t) (mmHg); output is heart rate
H(t) (bpm). Six states: carotid and aortic baroreceptor strain (eb_c, eb_a), baroreflex parasympathetic (Tpb), respiratory
parasympathetic (Tpr) and sympathetic (Ts) outflow, and heart rate H.

Two points copied from the authors' Fortran rather than the paper text:
- The sympathetic delay enters as dTs/dt = (-Ts(t - Ds) + Ks*Gs(t)) / tau_s (driver_baroreflex.f, F(5)).
  The history before t0 is the initial value (PHI).
- Inputs are interpolated with Numerical Recipes cubic splines with zero end slopes (spline(...,0,0,...)), held
  constant outside the data range (splint).

The authors used the stiff DDE solver RADAR5 (rtol = atol = 1e-8). This port uses fixed-step RK4 (default h = 0.5 ms)
with linear interpolation of the delayed Ts history.
"""
import math

import numpy as np

NAMES = ['A', 'B', 'Kb', 'Kpb', 'Kpr', 'Ks', 'taub', 'taupb', 'taupr', 'taus', 'tauH',
         'qw', 'qpb', 'qpr', 'qs', 'sw', 'spb', 'spr', 'ss', 'HI', 'Hpb', 'Hpr', 'Hs', 'Ds']
IDX = {n: i for i, n in enumerate(NAMES)}
STATE = ['eb_c', 'eb_a', 'Tpb', 'Tpr', 'Ts', 'H']


def nominal_parameters(Pbar, Pthbar, Hbar, HminR, HmaxR, age):
    """Nominal parameter vector (load_global.m) from resting means of SBP, Pth and HR, resting HR range and age."""
    A, B, Kb, Kpb, Kpr, Ks = 5.0, 0.5, 0.1, 5.0, 1.0, 5.0
    taub, taupb, taupr, taus, tauH = 0.9, 1.8, 6.0, 10.0, 0.5
    qw, qpb, qpr, qs = 0.04, 10.0, 1.0, 10.0
    Ds = 3.0
    sw, spr = Pbar, Pthbar
    HI = max(118 - 0.57 * age, Hbar)                 # intrinsic HR (age), not below resting mean
    HM = 208 - 0.7 * age                             # maximal HR (age)
    Hs = (1 / Ks) * (HM / HI - 1)
    ewc = _wall_strain(Pbar, A, qw, sw)
    ewa = _wall_strain(Pbar - Pthbar, A, qw, sw)
    n_ss = B * (ewc - Kb * ewc) + (1 - B) * (ewa - Kb * ewa)
    Tpb_ss, Ts_ss = 0.8, 0.2                         # 80 % parasympathetic, 20 % sympathetic at rest
    spb = n_ss + math.log(Kpb / Tpb_ss - 1) / qpb
    ss = n_ss - math.log(Ks / Ts_ss - 1) / qs
    Tpr_ss = Kpr / (1 + math.exp(qpr * (Pthbar - spr)))
    Hpr = (HmaxR - HminR) / HI / Tpr_ss
    Hpb = (1 - Hbar / HI + Hpr * Tpr_ss + Hs * Ts_ss) / Tpb_ss
    return np.array([A, B, Kb, Kpb, Kpr, Ks, taub, taupb, taupr, taus, tauH, qw, qpb, qpr, qs,
                     sw, spb, spr, ss, HI, Hpb, Hpr, Hs, Ds], float)


def _wall_strain(P, A, qw, sw):
    e = math.exp(-qw * (P - sw))
    return 1 - math.sqrt((1 + e) / (A + e))


def initial_conditions(p, P0, Pth0, H0):
    """Steady state for the first input samples (initialconditions.m); H starts at the measured H0."""
    A, B, Kb, Kpb, Kpr, Ks = p[0:6]
    qw, qpb, qpr, qs, sw, spb, spr, ss = p[11:19]
    ewc = _wall_strain(P0, A, qw, sw); ewa = _wall_strain(P0 - Pth0, A, qw, sw)
    n0 = B * (ewc - Kb * ewc) + (1 - B) * (ewa - Kb * ewa)
    return np.array([Kb * ewc, Kb * ewa,
                     Kpb / (1 + math.exp(-qpb * (n0 - spb))),
                     Kpr / (1 + math.exp(qpr * (Pth0 - spr))),
                     Ks / (1 + math.exp(qs * (n0 - ss))),
                     H0], float)


def nr_spline(x, y, yp1=0.0, ypn=0.0):
    """Second derivatives of a cubic spline with end slopes yp1, ypn (Numerical Recipes spline, as in spline.f)."""
    n = len(x); y2 = np.zeros(n); u = np.zeros(n)
    y2[0] = -0.5; u[0] = (3 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)
    for i in range(1, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]); pp = sig * y2[i - 1] + 2
        y2[i] = (sig - 1) / pp
        u[i] = (6 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / pp
    qn = 0.5; un = (3 / (x[-1] - x[-2])) * (ypn - (y[-1] - y[-2]) / (x[-1] - x[-2]))
    y2[-1] = (un - qn * u[-2]) / (qn * y2[-2] + 1)
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


def nr_splint(xa, ya, y2a, x):
    """Evaluate the spline at x (array), constant outside [xa[0], xa[-1]] (splint.f)."""
    x = np.asarray(x, float)
    khi = np.clip(np.searchsorted(xa, x, side='right'), 1, len(xa) - 1); klo = khi - 1
    h = xa[khi] - xa[klo]; a = (xa[khi] - x) / h; b = (x - xa[klo]) / h
    y = a * ya[klo] + b * ya[khi] + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * h ** 2 / 6
    return np.where(x <= xa[0], ya[0], np.where(x >= xa[-1], ya[-1], y))


def simulate(t, sbp, pth, p, y0, h=5e-4):
    """Integrate from t[0] to t[-1] with RK4 step h. Returns the states (len(t) x 6) at the input times t.
    sbp and pth are sampled at t; they are spline-interpolated as in the authors' code."""
    t = np.asarray(t, float); t0 = t[0]
    nsteps = int(round((t[-1] - t0) / h))
    tg = t0 + 0.5 * h * np.arange(2 * nsteps + 1)            # full and half steps
    P = nr_splint(t, sbp, nr_spline(t, np.asarray(sbp, float)), tg)
    PT = nr_splint(t, pth, nr_spline(t, np.asarray(pth, float)), tg)
    (A, B, Kb, Kpb, Kpr, Ks, taub, taupb, taupr, taus, tauH,
     qw, qpb, qpr, qs, sw, spb, spr, ss, HI, Hpb, Hpr, Hs, Ds) = [float(v) for v in p]
    ts_hist = np.empty(nsteps + 1); ts0 = float(y0[4])
    exp, sqrt = math.exp, math.sqrt

    def ts_lag(tt):                                           # Ts(tt - Ds), history = initial value
        s = (tt - Ds - t0) / h
        if s <= 0:
            return ts0
        i = int(s); w = s - i
        return ts_hist[i] if w == 0 else ts_hist[i] * (1 - w) + ts_hist[i + 1] * w

    def f(y, Pc, Pth, y5l):
        e = exp(-qw * (Pc - sw)); ewc = 1 - sqrt((1 + e) / (A + e))
        e = exp(-qw * (Pc - Pth - sw)); ewa = 1 - sqrt((1 + e) / (A + e))
        nb = B * (ewc - y[0]) + (1 - B) * (ewa - y[1])
        Gpb = 1 / (1 + exp(-qpb * (nb - spb))); Gs = 1 / (1 + exp(qs * (nb - ss)))
        Gpr = 1 / (1 + exp(qpr * (Pth - spr)))
        Ht = HI * (1 - Hpb * y[2] + Hpr * y[3] + Hs * y[4])
        return ((-y[0] + Kb * ewc) / taub, (-y[1] + Kb * ewa) / taub, (-y[2] + Kpb * Gpb) / taupb,
                (-y[3] + Kpr * Gpr) / taupr, (-y5l + Ks * Gs) / taus, (-y[5] + Ht) / tauH)

    y = [float(v) for v in y0]; ts_hist[0] = y[4]
    out_idx = np.rint((t - t0) / h).astype(int); out = np.empty((len(t), 6)); j = 0
    for k in range(nsteps + 1):
        while j < len(t) and out_idx[j] == k:
            out[j] = y; j += 1
        if k == nsteps:
            break
        tk = t0 + k * h
        l0, l1, l2 = ts_lag(tk), ts_lag(tk + 0.5 * h), ts_lag(tk + h)
        k1 = f(y, P[2 * k], PT[2 * k], l0)
        k2 = f([a + 0.5 * h * b for a, b in zip(y, k1)], P[2 * k + 1], PT[2 * k + 1], l1)
        k3 = f([a + 0.5 * h * b for a, b in zip(y, k2)], P[2 * k + 1], PT[2 * k + 1], l1)
        k4 = f([a + h * b for a, b in zip(y, k3)], P[2 * k + 2], PT[2 * k + 2], l2)
        y = [a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(y, k1, k2, k3, k4)]
        ts_hist[k + 1] = y[4]
    return out


def heart_rate(t, sbp, pth, p, y0, h=5e-4):
    """Model heart rate H(t) in bpm at the input times."""
    return simulate(t, sbp, pth, p, y0, h)[:, 5]
