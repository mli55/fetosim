"""PhysioNet rvmh1 (Response to Valsalva Maneuver in Humans, 27 healthy adults): reader and beat detector.
normNNN_BP: finger arterial pressure, 200 Hz; normNNN_P: expiratory pressure, 10 Hz (used as thoracic pressure).
The data are not in the repository: set RVMH1_DIR to the folder holding the .hea / .dat files
(https://physionet.org/content/response-to-valsalva-maneuver-in-humans/1.0.0/)."""
import os
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks


def root():
    d = os.environ.get('RVMH1_DIR')
    assert d, 'set RVMH1_DIR to the rvmh1 data folder'
    return Path(d)


def read(rec):
    """Signal (physical units) and sampling rate of one WFDB format-16 record."""
    h = (root() / (rec + '.hea')).read_text().splitlines()
    fs = float(h[0].split()[2]); f = h[1].split(); assert f[1] == '16', f
    gain, base = f[2].split('/')[0].split('('); gain = float(gain); base = float(base.rstrip(')'))
    x = np.fromfile(root() / f[0], dtype='<i2').astype(float)
    return (x - base) / gain, fs


def beats(b, fs):
    """Systolic peaks: time (s), systolic pressure (mmHg), heart rate (bpm) from the interval to the previous peak."""
    i, _ = find_peaks(b, distance=int(.3 * fs), prominence=15)
    t = i / fs; hr = np.r_[np.nan, 60 / np.diff(t)]
    return t, b[i], hr
