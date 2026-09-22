#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""Position, velocity (not yet), time: pseudoranges from the channels' observables, satellite
positions from the ephemerides, one least-squares solve.

The orbit and clock math is numpy-gps's (IS-GPS-200 Table 20-IV and §20.3.3.3.3.1), unchanged.
What differs from numpy-gps is how the transmit time is assembled - and it is simpler here.
numpy-gps decodes offline and must EXTRAPOLATE each satellite's time from a fitted anchor,
which leaves a whole-millisecond ambiguity per satellite it then searches over. A streaming
channel COUNTS: every code period since the anchored subframe is exactly 1 ms of that
satellite's clock, so

    t_sv(transmit, at the period boundary)  =  TOW_anchor + (epochs - epochs_anchor) * 1 ms

with no ambiguity at all. The three assembly laws numpy-gps paid for still hold:
  1. assemble in SV time, THEN apply the SV clock correction (af0 alone can be 0.5 ms = 150 km);
  2. a code phase is the offset to the NEXT code epoch (a channel reports the phase at its
     period boundary, where it is zero by construction, so this one costs nothing here);
  3. all satellites must be referred to ONE receive instant - the channels stamp every period
     boundary with the absolute input sample index, so each channel's transmit time is
     interpolated to a common sample.

Nothing here prints or stores a position. The `fix` dict carries ECEF metres and the caller
decides what to do with them; the panel shows counts, rms, DOP and a plausible-altitude flag."""
import numpy as np

MU = 3.986005e14
OMEGA_E = 7.2921151467e-5
C = 299792458.0
F_REL = -4.442807633e-10


# ---- orbit and clock (numpy-gps, verbatim in substance) ----------------------------------
def _wrap_week(t):
    if t > 302400:
        return t - 604800
    if t < -302400:
        return t + 604800
    return t


def ecc_anomaly(eph, t):
    A = eph["sqrtA"] ** 2
    tk = _wrap_week(t - eph["toe"])
    M = eph["M0"] + (np.sqrt(MU / A ** 3) + eph.get("dn", 0.0)) * tk
    E = M
    for _ in range(15):
        E = M + eph["e"] * np.sin(E)
    return E


def sat_ecef(eph, t):
    """Satellite ECEF (m) at GPS time-of-week t (IS-GPS-200 Table 20-IV)."""
    A = eph["sqrtA"] ** 2
    tk = _wrap_week(t - eph["toe"])
    E = ecc_anomaly(eph, t)
    e = eph["e"]
    nu = np.arctan2(np.sqrt(1 - e * e) * np.sin(E), np.cos(E) - e)
    phi = nu + eph["omega"]
    s2, c2 = np.sin(2 * phi), np.cos(2 * phi)
    u = phi + eph["Cus"] * s2 + eph["Cuc"] * c2
    r = A * (1 - e * np.cos(E)) + eph["Crs"] * s2 + eph["Crc"] * c2
    i = eph["i0"] + eph["IDOT"] * tk + eph["Cis"] * s2 + eph["Cic"] * c2
    xo, yo = r * np.cos(u), r * np.sin(u)
    Om = eph["Omega0"] + (eph["OmegaDot"] - OMEGA_E) * tk - OMEGA_E * eph["toe"]
    return np.array([xo * np.cos(Om) - yo * np.cos(i) * np.sin(Om),
                     xo * np.sin(Om) + yo * np.cos(i) * np.cos(Om),
                     yo * np.sin(i)])


def clock_corr(eph, t_sv):
    """SV clock offset (s) incl. the relativistic term. t_gps = t_sv - clock_corr."""
    dt = _wrap_week(t_sv - eph["toc"])
    dtr = F_REL * eph["e"] * eph["sqrtA"] * np.sin(ecc_anomaly(eph, t_sv))
    return eph["af0"] + eph["af1"] * dt + eph["af2"] * dt * dt + dtr


def ecef_to_llh(p):
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    x, y, z = p
    lon = np.arctan2(y, x)
    r = np.hypot(x, y)
    lat = np.arctan2(z, r * (1 - e2))
    for _ in range(8):
        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = r / np.cos(lat) - N
        lat = np.arctan2(z, r * (1 - e2 * N / (N + h)))
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.degrees(lat), np.degrees(lon), r / np.cos(lat) - N


def az_el(rx, sp):
    lat, lon, _ = ecef_to_llh(rx)
    la, lo = np.radians(lat), np.radians(lon)
    d = sp - rx
    e = np.array([-np.sin(lo), np.cos(lo), 0.0]) @ d
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)]) @ d
    u = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]) @ d
    return np.arctan2(e, n) % (2 * np.pi), np.arcsin(u / np.linalg.norm(d))


def tropo_delay(el):
    return 2.47 / (np.sin(el) + 0.0121)


def klobuchar(a, b, lat_deg, lon_deg, az, el, t_gps):
    E = el / np.pi
    psi = 0.0137 / (E + 0.11) - 0.022
    phi_i = np.clip(lat_deg / 180.0 + psi * np.cos(az), -0.416, 0.416)
    lam_i = lon_deg / 180.0 + psi * np.sin(az) / np.cos(phi_i * np.pi)
    phi_m = phi_i + 0.064 * np.cos((lam_i - 1.617) * np.pi)
    t = (4.32e4 * lam_i + t_gps) % 86400.0
    amp = max(sum(a[k] * phi_m ** k for k in range(4)), 0.0)
    per = max(sum(b[k] * phi_m ** k for k in range(4)), 72000.0)
    x = 2.0 * np.pi * (t - 50400.0) / per
    F = 1.0 + 16.0 * (0.53 - E) ** 3
    return F * (5e-9 + amp * (1.0 - x * x / 2.0 + x ** 4 / 24.0)) if abs(x) < 1.57 else F * 5e-9


# ---- transmit time from a channel's observable -------------------------------------------
def transmit_time_sv(obs, anchor, fs):
    """SV-clock transmit time (s of week) of the signal arriving at obs['sample_abs'].

    anchor = (epochs_anchor, tow_anchor): the code period on which an anchored subframe's
    first bit began, and that subframe's TOW. The channel's period boundary at obs is
    (obs['epochs'] - epochs_anchor) periods later - exactly that many milliseconds of SV time.
    The channel reports the code phase AT the boundary, where it is 0 by construction; a
    non-zero value (a period cut mid-way) is the offset to the NEXT epoch: subtract it (law 2)."""
    e_anchor, tow_anchor = anchor
    frac = obs.get("code_phase", 0.0) / 1023.0             # fraction of a period, [0, 1)
    return tow_anchor + (obs["epochs"] - e_anchor) * 1e-3 - frac * 1e-3


def refer_to_common_sample(entries, fs):
    """entries: list of dicts {prn, eph, t_sv, sample_abs, rate_hz}. Each channel's transmit time is
    at its own boundary sample; slide each to the latest boundary among them by its own code
    rate (SV time advances at (1 + doppler/L1) x receiver time). Returns t_common and the list
    with t_sv adjusted (law 3)."""
    s_ref = max(e["sample_abs"] for e in entries)
    out = []
    for e in entries:
        dt_rx = (s_ref - e["sample_abs"]) / fs
        drift = 1.0 + e.get("carrier_hz", 0.0) / 1575.42e6
        out.append(dict(e, t_sv=e["t_sv"] + dt_rx * drift))
    return s_ref, out


# ---- the solve ----------------------------------------------------------------------------
def solve_ls(sats, weights=None):
    """sats = [(ecef_xyz, pseudorange_m)] -> (x, y, z, c*dt), Gauss-Newton."""
    SP = np.asarray([sp for sp, _ in sats], float)
    PR = np.asarray([pr for _, pr in sats], float)
    n = len(PR)
    x = np.zeros(4)
    A = np.empty((n, 4))
    A[:, 3] = 1.0
    sw = None if weights is None else np.sqrt(np.asarray(weights, float))[:, None]
    for _ in range(12):
        d = x[:3] - SP
        rng = np.sqrt((d * d).sum(axis=1))
        A[:, :3] = d / rng[:, None]
        res = PR - (rng + x[3])
        dx = np.linalg.lstsq(A if sw is None else A * sw, res if sw is None else res * sw[:, 0], rcond=None)[0]
        x = x + dx
        if np.linalg.norm(dx[:3]) < 1e-3:
            break
    return x, A


def solve(entries, iono=None, weights=None):
    """entries: [{prn, eph, t_sv}] all referred to one receive instant (SV clock transmit
    times). Returns a fix dict: ecef, llh, rms_m, dop, n, residuals, clock_bias_s, t_rx."""
    prs = []
    for e in entries:
        t_gps = e["t_sv"] - clock_corr(e["eph"], e["t_sv"])              # law 1: correct AFTER assembly
        prs.append((e["prn"], e["eph"], t_gps))
    t_rx = max(t for _, _, t in prs) + 0.075
    x = np.zeros(4)
    for it in range(8):
        sats = []
        for prn, eph, t_tx in prs:
            sp = sat_ecef(eph, t_tx)
            tau = max(t_rx - t_tx, 0.0)
            th = OMEGA_E * tau                                            # Sagnac
            rot = np.array([[np.cos(th), np.sin(th), 0], [-np.sin(th), np.cos(th), 0], [0, 0, 1]])
            pr = C * (t_rx - t_tx)
            if it >= 2 and np.linalg.norm(x[:3]) > 1e6:                     # atmosphere, once we know where we are
                az, el = az_el(x[:3], rot @ sp)
                pr -= tropo_delay(el)
                if iono:
                    la, lo, _ = ecef_to_llh(x[:3])
                    pr -= C * klobuchar(iono["a"], iono["b"], la, lo, az, el, t_rx)
            sats.append((rot @ sp, pr))
        x, A = solve_ls(sats, weights)
        t_rx -= x[3] / C
    res = np.array([pr - (np.linalg.norm(x[:3] - sp) + x[3]) for sp, pr in sats])
    lat, lon, h = ecef_to_llh(x[:3])
    try:
        Q = np.linalg.inv(A.T @ A)
        pdop = float(np.sqrt(np.trace(Q[:3, :3])))
    except np.linalg.LinAlgError:
        pdop = float("nan")
    return {"ecef": x[:3].tolist(), "llh": (float(lat), float(lon), float(h)), "clock_bias_s": float(x[3] / C),
            "t_rx": float(t_rx), "rms_m": float(np.sqrt(np.mean(res ** 2))), "pdop": pdop, "n": len(sats),
            "residuals_m": res.tolist(), "prns": [p for p, _, _ in prs],
            "altitude_plausible": bool(-500 < h < 9000)}


def fix_from_channels(channels, fs, iono=None):
    """channels: list of {prn, eph, anchor:(epochs, tow), obs:{epochs, code_phase, sample_abs,
    carrier_hz}}. The whole path: transmit times, common instant, solve."""
    entries = []
    for c in channels:
        if c.get("eph") is None or c.get("anchor") is None or c.get("obs") is None:
            continue
        t_sv = transmit_time_sv(c["obs"], c["anchor"], fs)
        entries.append(dict(prn=c["prn"], eph=c["eph"], t_sv=t_sv, sample_abs=c["obs"]["sample_abs"],
                            carrier_hz=c["obs"].get("carrier_hz", 0.0)))
    if len(entries) < 4:
        return None
    s_ref, entries = refer_to_common_sample(entries, fs)
    fx = solve(entries, iono=iono)
    fx["sample_abs"] = int(s_ref)
    return fx
