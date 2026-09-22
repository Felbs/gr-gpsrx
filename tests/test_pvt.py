# SPDX-License-Identifier: GPL-3.0-or-later
"""The PVT solver on a synthetic constellation: put a receiver somewhere, invent six orbits
with a spread of geometry, compute what each channel WOULD observe (transmit time at a common
receive instant, plus SV clock offsets), and require the solver to get the receiver back.
Gate 3 proper (against numpy-gps on real air) needs a capture; this is the solver's own truth."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import pvt  # noqa: E402

FS = 2.048e6
# a plausible receiver: mid-latitude, near sea level (an arbitrary point in the North Atlantic)
RX_LLH = (45.0, -30.0, 120.0)


def llh_to_ecef(lat, lon, h):
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    la, lo = np.radians(lat), np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(la) ** 2)
    return np.array([(N + h) * np.cos(la) * np.cos(lo), (N + h) * np.cos(la) * np.sin(lo), (N * (1 - e2) + h) * np.sin(la)])


def constellation(t0, rx):
    """Six GPS-like orbits chosen so the satellites are ABOVE the receiver at t0: for each, pick
    an argument of latitude that puts the satellite roughly over a target direction and a node
    that lines the plane up, then keep it if it is really visible. Inventing elements blind put
    two of six above the horizon (measured); this places them on purpose."""
    base = dict(af1=0.0, af2=0.0, Crs=0.0, dn=0.0, Cuc=0.0, Cus=0.0, Cic=0.0, Cis=0.0, Crc=0.0, IDOT=0.0,
                sqrtA=5153.7, e=0.005, toe=t0, toc=t0, WN=345, i0=0.96, OmegaDot=-8e-9)
    lat, lon, _ = pvt.ecef_to_llh(rx)
    eph, rng = [], np.random.default_rng(7)
    k = 0
    while len(eph) < 6 and k < 400:
        k += 1
        # a satellite over a point within ~35 degrees of the receiver is well above its horizon
        dlat, dlon = rng.uniform(-30, 30), rng.uniform(-40, 40)
        sub_lat, sub_lon = np.radians(lat + dlat), np.radians(lon + dlon)
        # for inclination i0 the sub-satellite latitude fixes the argument of latitude u:
        # sin(sub_lat) = sin(i0) sin(u); the node then follows from the longitude
        if abs(np.sin(sub_lat)) > np.sin(base["i0"]):
            continue
        u = np.arcsin(np.sin(sub_lat) / np.sin(base["i0"]))
        if rng.random() < 0.5:
            u = np.pi - u
        # longitude of the ascending node in ECEF terms at t0: lon_sub = Om_ecef + atan2(cos i sin u, cos u)
        Om_ecef = sub_lon - np.arctan2(np.cos(base["i0"]) * np.sin(u), np.cos(u))
        Om0 = Om_ecef + pvt.OMEGA_E * t0                       # sat_ecef subtracts OMEGA_E * toe
        om = rng.uniform(0, 2 * np.pi)
        e = dict(base, prn=len(eph) + 1, Omega0=Om0, M0=u - om, omega=om, af0=rng.uniform(-3e-4, 3e-4))
        _, el = pvt.az_el(rx, pvt.sat_ecef(e, t0))
        if el > np.radians(15):
            eph.append(e)
    return eph


def visible(rx, ephs, t):
    return [e for e in ephs if pvt.az_el(rx, pvt.sat_ecef(e, t))[1] > np.radians(10)]


def observe(rx, eph, t_rx):
    """What a channel would report at receive time t_rx (GPS): the SV transmit time in the SV's
    own clock, i.e. t_tx(gps) + clock offset - found by iterating the light-time equation with
    the Earth rotating under the signal (Sagnac)."""
    t_tx = t_rx - 0.07
    for _ in range(10):
        sp = pvt.sat_ecef(eph, t_tx)
        tau = t_rx - t_tx
        th = pvt.OMEGA_E * tau
        rot = np.array([[np.cos(th), np.sin(th), 0], [-np.sin(th), np.cos(th), 0], [0, 0, 1]])
        _, el = pvt.az_el(rx, rot @ sp)
        t_tx = t_rx - (np.linalg.norm(rot @ sp - rx) + pvt.tropo_delay(el)) / pvt.C     # the air adds delay
    return t_tx + pvt.clock_corr(eph, t_tx)        # to first order the SV-clock reading at transmit


def test_solver_recovers_the_receiver_from_six_satellites():
    t0 = 302400.0
    rx = llh_to_ecef(*RX_LLH)
    ephs = constellation(t0, rx)
    assert len(ephs) == 6, len(ephs)
    t_rx = t0 + 100.0
    ephs = visible(rx, ephs, t_rx)
    assert len(ephs) >= 5, len(ephs)
    entries = [dict(prn=e["prn"], eph=e, t_sv=observe(rx, e, t_rx)) for e in ephs]
    fx = pvt.solve(entries)
    err = np.linalg.norm(np.array(fx["ecef"]) - rx)
    assert err < 0.5, (err, fx["rms_m"])
    assert fx["rms_m"] < 0.5
    assert fx["altitude_plausible"]


def test_transmit_time_from_a_channel_count_and_common_instant():
    """Two channels anchored at different periods whose epochs arrived at different samples must
    be referred to one instant: an epoch that arrived 2048 samples earlier is slid 1 ms forward."""
    anchor_a = (1000, 302400.0)
    anchor_b = (3500, 302406.0)
    obs_a = {"epochs": 1000 + 8000, "epoch_sample": 20_480_000.25, "carrier_hz": 0.0}
    obs_b = {"epochs": 3500 + 5500, "epoch_sample": 20_480_000.25 - 2048, "carrier_hz": 0.0}
    ta = pvt.transmit_time_sv(obs_a, anchor_a, FS)
    tb = pvt.transmit_time_sv(obs_b, anchor_b, FS)
    assert abs(ta - 302408.0) < 1e-9 and abs(tb - 302411.5) < 1e-9
    s_ref, out = pvt.refer_to_common_sample([dict(prn=1, eph=None, t_sv=ta, epoch_sample=obs_a["epoch_sample"]),
                                             dict(prn=2, eph=None, t_sv=tb, epoch_sample=obs_b["epoch_sample"])], FS)
    assert s_ref == obs_a["epoch_sample"]
    assert abs(out[1]["t_sv"] - (302411.5 + 1e-3)) < 1e-12


def test_observable_reports_the_epochs_fractional_arrival_sample():
    """The channel's period ends on a whole sample a little AFTER the code epoch; the observable
    must give the epoch's own arrival, fractional. A third of a chip is 100 m: this matters."""
    from gpsrx import synth
    from gpsrx.track import Channel
    x, tr = synth.satellite(9, FS, 0.3, doppler_hz=500.0, code_phase_samples=1000)
    ch = Channel(9, FS, 500.0, 1000)
    xd = x.astype(np.complex128)
    pos = 0
    for _ in range(250):
        n = ch.samples_needed()
        ch.step(xd[pos:pos + n])
        pos += n
    o = ch.observable()
    assert o["samples_in"] - 2 < o["epoch_sample"] <= o["samples_in"]
    # the synthetic code started at sample 1000 and runs at tr['code_rate']: epoch k arrives at
    # 1000 + k * 1023 / rate * fs. The tracked channel's epoch count is offset by however many
    # periods the handover cost; check the arrival against the nearest true epoch.
    true = 1000 + np.arange(0, 400) * 1023.0 / tr["code_rate"] * FS
    assert np.min(np.abs(true - o["epoch_sample"])) < 0.3, np.min(np.abs(true - o["epoch_sample"]))   # < 0.15 chip


def test_a_wrong_clock_correction_order_would_be_kilometres_off():
    """Law 1, as a number: applying af0 = 2.1e-4 s BEFORE the integer-ms assembly would move that
    satellite's pseudorange by c * af0 = 63 km. The solver applies it after; the truth generator
    applies it consistently, so a correct fix means the order is right."""
    assert abs(pvt.C * 2.1e-4 - 62956) < 1
