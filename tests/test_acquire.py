# SPDX-License-Identifier: GPL-3.0-or-later
"""Acquisition on a synthetic sky: every satellite that is there is found at its code phase and
Doppler, nothing that is not there is reported."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import acquire, synth  # noqa: E402

FS = 2.048e6


def test_finds_every_synthetic_satellite_and_no_other():
    birds = [dict(prn=3, doppler_hz=1234.0, code_phase_samples=700, cn0_dbhz=45),
             dict(prn=7, doppler_hz=-2510.0, code_phase_samples=1500, cn0_dbhz=42),
             dict(prn=12, doppler_hz=90.0, code_phase_samples=10, cn0_dbhz=40),
             dict(prn=25, doppler_hz=4020.0, code_phase_samples=2000, cn0_dbhz=38)]
    x, _ = synth.sky(FS, 0.2, birds)
    x = x.astype(np.complex128)
    found = acquire.sky(x, FS, threshold=2.5, n_noncoh=40)
    got = {r["prn"]: r for r in found}
    assert set(got) == {b["prn"] for b in birds}, sorted(got)
    for b in birds:
        r = got[b["prn"]]
        assert abs(r["code_phase"] - b["code_phase_samples"]) < 1.0, (b, r)
        assert abs(r["doppler_hz"] - b["doppler_hz"]) < 15.0, (b, r)          # refined to ~10 Hz
    # the metric separates: every real bird well above every absent PRN
    allm = acquire.search(x, FS, n_noncoh=40)
    absent = max(v["metric"] for p, v in allm.items() if p not in got)
    present = min(v["metric"] for p, v in allm.items() if p in got)
    assert present > 1.5 * absent, (present, absent)
