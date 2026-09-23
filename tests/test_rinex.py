# SPDX-License-Identifier: GPL-3.0-or-later
"""RINEX 3 output: GPS time to calendar against the known epochs, the observation file round
trip, and - when georinex is installed - both files read back by an independent parser."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import pvt, rinex  # noqa: E402
from test_pvt import RX_LLH, constellation, llh_to_ecef, observe, visible  # noqa: E402


def test_gps_time_to_calendar_at_the_known_epochs():
    assert rinex.gps_to_calendar(0, 0.0)[:3] == (1980, 1, 6)
    assert rinex.gps_to_calendar(1024, 0.0)[:3] == (1999, 8, 22)                 # the first rollover
    assert rinex.gps_to_calendar(2048, 0.0)[:3] == (2019, 4, 7)                  # the second
    y, mo, d, h, mi, s = rinex.gps_to_calendar(2048, 3 * 86400 + 12 * 3600 + 34 * 60 + 56.25)
    assert (y, mo, d, h, mi) == (2019, 4, 10, 12, 34) and abs(s - 56.25) < 1e-9
    assert rinex.full_week(381) == 2429


def _observables_at(t_rx, rx, ephs, cyc0):
    vis = visible(rx, ephs, t_rx)
    entries = [dict(prn=e["prn"], eph=e, t_sv=observe(rx, e, t_rx), carrier_cycles=cyc0[e["prn"]] - 1000.0 * (t_rx % 7),
                    carrier_hz=-1000.0, cn0_db=44.0) for e in vis]
    fx = pvt.solve(entries)
    obs = [dict(sys="GPS", prn=e["prn"], pr_m=pvt.C * (fx["t_rx"] - e["t_sv"]), phase_cyc=-e["carrier_cycles"],
                doppler_hz=e["carrier_hz"], cn0_db=e["cn0_db"]) for e in entries]
    return fx, obs, vis


def test_observation_file_round_trip(tmp_path):
    """Two epochs of the synthetic constellation written and read back: the pseudoranges, phases,
    Dopplers and C/N0 survive to the millimetre, and the pseudoranges are the ranges a solver
    wants (their spread across satellites matches the geometry, not a clock)."""
    t0 = 302400.0
    rx = llh_to_ecef(*RX_LLH)
    ephs = constellation(t0, rx)
    cyc0 = {e["prn"]: 12345.0 * e["prn"] for e in ephs}
    path = str(tmp_path / "test.obs")
    wr = rinex.ObsWriter(path)
    written = []
    for t_rx in (t0 + 100.0, t0 + 101.0):
        fx, obs, _ = _observables_at(t_rx, rx, ephs, cyc0)
        wr.add_epoch(2429, fx["t_rx"], obs)
        written.append((fx["t_rx"], obs))
    wr.close()
    back = rinex.read_obs(path)
    assert len(back) == 2
    for (t_rx, obs), (sod, rows) in zip(written, back):
        assert abs(sod - (t_rx % 86400)) < 1e-6
        assert len(rows) == len(obs)
        by_prn = {r[1]: r for r in rows}
        for o in obs:
            r = by_prn[o["prn"]]
            assert r[0] == "G"
            assert abs(r[2] - o["pr_m"]) < 1e-3 and abs(r[3] - o["phase_cyc"]) < 1e-3
            assert abs(r[4] - o["doppler_hz"]) < 1e-3 and abs(r[5] - o["cn0_db"]) < 1e-3
        prs = np.array([o["pr_m"] for o in obs])
        assert 19e6 < prs.min() and prs.max() < 27e6                             # ranges to the GPS shell


def test_navigation_file_has_every_field(tmp_path):
    t0 = 302400.0
    ephs = constellation(t0, llh_to_ecef(*RX_LLH))
    for e in ephs:
        e.update(WN=381, IODE2=17, IODC=17, TGD=-1.2e-8, health=0, URA=2, sys="GPS")
    path = str(tmp_path / "test.nav")
    rinex.write_nav(path, ephs)
    lines = open(path).read().splitlines()
    body = lines[lines.index(next(l for l in lines if l.endswith("END OF HEADER"))) + 1:]
    assert len(body) == 8 * len(ephs)
    assert body[0].startswith(f"G{ephs[0]['prn']:02d} 2026 ")
    for k in range(len(ephs)):
        assert all(len(body[8 * k + i]) >= 4 + 19 * 4 for i in range(1, 7)), body[8 * k + 1]


def test_georinex_reads_both_files(tmp_path):
    gr = pytest.importorskip("georinex")
    t0 = 302400.0
    rx = llh_to_ecef(*RX_LLH)
    ephs = constellation(t0, rx)
    for e in ephs:
        e.update(WN=381, IODE2=17, IODC=17, TGD=0.0, health=0, URA=2, sys="GPS")
    cyc0 = {e["prn"]: 0.0 for e in ephs}
    obs_path, nav_path = str(tmp_path / "t.obs"), str(tmp_path / "t.nav")
    wr = rinex.ObsWriter(obs_path)
    fx, obs, _ = _observables_at(t0 + 100.0, rx, ephs, cyc0)
    wr.add_epoch(2429, fx["t_rx"], obs)
    wr.close()
    rinex.write_nav(nav_path, ephs)
    o = gr.load(obs_path)
    assert "C1C" in o and "L1C" in o and o.sizes["time"] == 1
    assert sorted(int(s[1:]) for s in o.sv.values) == sorted(x["prn"] for x in obs)
    n = gr.load(nav_path)
    assert "sqrtA" in n and n.sizes["sv"] == len(ephs)
    assert np.allclose(sorted(float(v) for v in n["sqrtA"].values.ravel() if np.isfinite(v)), sorted(e["sqrtA"] for e in ephs))
