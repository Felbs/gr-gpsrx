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
import os

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
    M = eph["M0"] + (np.sqrt(eph.get("mu", MU) / A ** 3) + eph.get("dn", 0.0)) * tk
    E = M
    for _ in range(15):
        E = M + eph["e"] * np.sin(E)
    return E


def sat_ecef(eph, t):
    """Satellite ECEF (m) at time-of-week t (IS-GPS-200 Table 20-IV; Galileo's ICD 5.1.1 is the same
    algorithm with its own mu, carried in the ephemeris as 'mu')."""
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
    """SV-clock transmit time (s of week) of the code epoch `obs['epochs']`, which arrived at the
    receiver at fractional input sample obs['epoch_sample'].

    anchor = (epochs_anchor, tow_anchor): the code period on which an anchored subframe's first
    bit began, and that subframe's TOW. Epoch `epochs` is (epochs - epochs_anchor) periods later:
    exactly that many milliseconds of SV time. No fraction: the channel reports the epoch's own
    arrival sample, so the range fraction lives in the receive time, not here."""
    e_anchor, tow_anchor = anchor
    return tow_anchor + (obs["epochs"] - e_anchor) * obs.get("period_s", 1e-3)   # 1 ms GPS, 4 ms Galileo


def refer_to_common_sample(entries, fs):
    """entries: [{prn, eph, t_sv, epoch_sample, carrier_hz}]: each channel's transmit time is that of
    the code epoch arriving at its own (fractional) sample. Slide every one to the latest of those
    samples: the satellite's clock advances by the receiver interval scaled by (1 + Doppler/L1).
    Law 3: one receive instant for all. Returns (s_ref, entries)."""
    s_ref = max(e["epoch_sample"] for e in entries)
    out = []
    for e in entries:
        dt_rx = (s_ref - e["epoch_sample"]) / fs
        drift = 1.0 + e.get("carrier_hz", 0.0) / 1575.42e6
        out.append(dict(e, t_sv=e["t_sv"] + dt_rx * drift,
                        carrier_cycles=e.get("carrier_cycles", 0.0) + e.get("carrier_hz", 0.0) * dt_rx))
    return s_ref, out


# ---- the solve ----------------------------------------------------------------------------
def solve_ls(sats, weights=None, isb=None):
    """sats = [(ecef_xyz, pseudorange_m)] -> (x, y, z, c*dt[, c*isb]), Gauss-Newton. isb: a 0/1
    vector marking the satellites of a second system (Galileo), which get their own clock
    unknown - the inter-system bias, GST minus GPS time as this receiver sees it."""
    SP = np.asarray([sp for sp, _ in sats], float)
    PR = np.asarray([pr for _, pr in sats], float)
    n = len(PR)
    k = 5 if isb is not None and np.any(isb) and not np.all(isb) else 4
    x = np.zeros(k)
    A = np.empty((n, k))
    A[:, 3] = 1.0
    if k == 5:
        A[:, 4] = np.asarray(isb, float)
    sw = None if weights is None else np.sqrt(np.asarray(weights, float))[:, None]
    for _ in range(12):
        d = x[:3] - SP
        rng = np.sqrt((d * d).sum(axis=1))
        A[:, :3] = d / rng[:, None]
        res = PR - (rng + A[:, 3:] @ x[3:])
        dx = np.linalg.lstsq(A if sw is None else A * sw, res if sw is None else res * sw[:, 0], rcond=None)[0]
        x = x + dx
        if np.linalg.norm(dx[:3]) < 1e-3:
            break
    return x, A


# chi-square 99% quantiles by degrees of freedom, for RAIM (no scipy in the engine)
CHI2_99 = {1: 6.63, 2: 9.21, 3: 11.34, 4: 13.28, 5: 15.09, 6: 16.81, 7: 18.48, 8: 20.09, 9: 21.67, 10: 23.21}
SIGMA0_M = 3.0                  # pseudorange noise at zenith and strong signal, metres


def measurement_variance(el, cn0_db):
    """The weight of one pseudorange, RTKLIB's shape: sigma^2 = sigma0^2 (1 + 1/sin^2 el), so a
    satellite at 10 degrees counts a thirtieth of one at zenith (more air, more multipath),
    de-weighted further by 10^((42 - C/N0)/10) below 42 dB-Hz (a 32 dB-Hz signal counts a tenth)."""
    se = max(np.sin(el), 0.05)
    v = SIGMA0_M ** 2 * (1.0 + 1.0 / (se * se))
    if cn0_db and cn0_db > 0:
        v *= 10 ** (max(0.0, 42.0 - cn0_db) / 10.0)
    return v


def _solve_once(prs, iono, t_rx0, weighted=True):
    t_rx = t_rx0
    x = np.zeros(4)
    weights = None
    isb = np.array([1.0 if len(p) > 4 and p[4] == "GAL" else 0.0 for p in prs])
    for it in range(8):
        sats, vars_ = [], []
        for p in prs:
            prn, eph, t_tx, cn0 = p[:4]
            sp = sat_ecef(eph, t_tx)
            tau = max(t_rx - t_tx, 0.0)
            th = OMEGA_E * tau                                            # Sagnac
            rot = np.array([[np.cos(th), np.sin(th), 0], [-np.sin(th), np.cos(th), 0], [0, 0, 1]])
            pr = C * (t_rx - t_tx)
            # atmosphere only once the estimate is near the Earth's surface (6.3-6.5 Mm from the
            # centre) and the satellite is above the horizon; a Klobuchar term evaluated at a
            # nonsense elevation returned NaN and the least squares then failed to converge
            rn = np.linalg.norm(x[:3])
            el = None
            if it >= 2 and 6.2e6 < rn < 6.6e6:
                az, el = az_el(x[:3], rot @ sp)
                if el > 0.02:
                    pr -= tropo_delay(el)
                    if iono:
                        la, lo, _ = ecef_to_llh(x[:3])
                        d_i = C * klobuchar(iono["a"], iono["b"], la, lo, az, el, t_rx)
                        if np.isfinite(d_i):
                            pr -= d_i
            sats.append((rot @ sp, pr))
            vars_.append(measurement_variance(el if el is not None else np.pi / 4, cn0))
        # weights only once the elevations are known (the first passes are unweighted)
        weights = [1.0 / v for v in vars_] if (weighted and it >= 3) else None
        x, A = solve_ls(sats, weights, isb=isb)
        t_rx -= x[3] / C
    res = np.array([pr - (np.linalg.norm(x[:3] - sp) + A[i, 3:] @ x[3:]) for i, (sp, pr) in enumerate(sats)])
    w = np.array(weights) if weights is not None else np.ones(len(sats))
    return x, A, sats, res, w, t_rx


def solve(entries, iono=None, weights=None, raim=True, mask_deg=0.0):
    """entries: [{prn, eph, t_sv, cn0_db?}] all referred to one receive instant (SV clock transmit
    times). Weighted least squares (elevation + C/N0), then RAIM: a chi-square test on the
    weighted residuals; if it fails with six or more satellites, drop one at a time and keep the
    solution with the smallest normalised residual (gnss-sdr / RTKLIB's fault detection and
    exclusion). mask_deg > 0: satellites below that elevation (from a first solve) are left out
    when enough remain - low satellites carry the most multipath and troposphere error, and the
    weighting only demotes them. Returns a fix dict: ecef, llh, rms_m, pdop, n, residuals,
    clock_bias_s, t_rx, raim (test statistic, threshold, excluded PRN), masked."""
    prs = []
    for e in entries:
        t_gps = e["t_sv"] - clock_corr(e["eph"], e["t_sv"])              # law 1: correct AFTER assembly
        prs.append((e["prn"], e["eph"], t_gps, float(e.get("cn0_db", 0.0) or 0.0), e.get("sys", "GPS")))
    t_rx0 = max(p[2] for p in prs) + 0.075
    weighted = weights is None or weights is not False

    def run(subset):
        x, A, sats, res, w, t_rx = _solve_once(subset, iono, t_rx0, weighted)
        dof = len(subset) - (5 if 0 < sum(1 for p in subset if p[4] == "GAL") < len(subset) else 4)
        stat = float(np.sum(w * res * res)) if dof > 0 else 0.0
        return dict(x=x, A=A, sats=sats, res=res, w=w, t_rx=t_rx, stat=stat, dof=dof, prs=subset)

    best = run(prs)
    masked = []
    if mask_deg > 0 and len(prs) > 4:
        els = [np.degrees(az_el(best["x"][:3], sp)[1]) for sp, _ in best["sats"]]
        keep = [p for p, el in zip(prs, els) if el >= mask_deg]
        # keep one more than the solve needs, so RAIM keeps its redundancy: masking a 5-satellite
        # sky down to 4 moved the fix 26 m and quadrupled the scatter (Smith Point capture) -
        # the low satellite was worth more than its error
        n_min = 6 if 0 < sum(1 for p in keep if p[4] == "GAL") < len(keep) else 5
        if len(keep) < len(prs) and len(keep) >= n_min:
            masked = [("E" if p[4] == "GAL" else "") + str(p[0]) for p, el in zip(prs, els) if el < mask_deg]
            prs = keep
            best = run(prs)
    excluded = None
    thr = CHI2_99.get(best["dof"], 6.63 + 2.3 * max(best["dof"] - 1, 0))
    raim_pass = best["dof"] <= 0 or best["stat"] <= thr
    if raim and not raim_pass and len(prs) >= 6:
        cands = []
        for k in range(len(prs)):
            sub = prs[:k] + prs[k + 1:]
            r = run(sub)
            cands.append((r["stat"] / max(r["dof"], 1), prs[k][0], r))
        cands.sort(key=lambda c: c[0])
        norm, prn_out, r = cands[0]
        thr_sub = CHI2_99.get(r["dof"], 6.63 + 2.3 * max(r["dof"] - 1, 0))
        if r["stat"] <= thr_sub:                      # excluding this one makes the rest consistent
            best, excluded, raim_pass = r, prn_out, True
    x, A, sats, res, w, t_rx, prs = best["x"], best["A"], best["sats"], best["res"], best["w"], best["t_rx"], best["prs"]
    lat, lon, h = ecef_to_llh(x[:3])
    azel = {}
    for p, (sp, _) in zip(prs, sats):
        az, el = az_el(x[:3], sp)
        azel[("E" if p[4] == "GAL" else "") + str(p[0])] = (float(np.degrees(az)) % 360.0, float(np.degrees(el)))
    try:
        Aw = A * np.sqrt(w)[:, None]
        Q = np.linalg.inv(Aw.T @ Aw) * float(np.mean(1.0 / w))      # in metres^2, comparable to unweighted
        pdop = float(np.sqrt(np.trace(np.linalg.inv(A.T @ A)[:3, :3])))
        cov = Q[:3, :3]
    except np.linalg.LinAlgError:
        pdop, cov = float("nan"), np.full((3, 3), np.nan)
    return {"ecef": x[:3].tolist(), "llh": (float(lat), float(lon), float(h)), "clock_bias_s": float(x[3] / C),
            "t_rx": float(t_rx), "rms_m": float(np.sqrt(np.mean(res ** 2))), "pdop": pdop, "n": len(sats),
            "residuals_m": res.tolist(), "prns": [("E" if p[4] == "GAL" else "") + str(p[0]) for p in prs], "azel": azel,
            "isb_m": float(x[4]) if len(x) > 4 else None,
            "isb_s": float(-x[4] / C) if len(x) > 4 else None,       # GST - GPS time as this receiver sees it
            "weights": (w / w.max()).tolist(), "cov_ecef": cov.tolist(),
            "raim": {"stat": best["stat"], "threshold": thr, "pass": bool(raim_pass), "excluded": excluded},
            "masked": masked,
            "altitude_plausible": bool(-500 < h < 9000),
            # valid: on the planet AND the residuals are consistent (or there were too few
            # satellites to tell). A detected fault that exclusion could not isolate is reported,
            # but not averaged, not filtered, not believed (a drive leg: six satellites, 150 m rms,
            # no single exclusion consistent - two bad anchors at once, measured)
            "valid": bool(-500 < h < 9000 and raim_pass)}


LAMBDA_L1 = C / 1575.42e6           # 0.1903 m per carrier cycle


def hatch_smooth(entries, state, fs, s_ref, M=100, slip_m=30.0):
    """Carrier smoothing of the code pseudorange (Hatch 1982; gnss-sdr's enable_carrier_smoothing):
    the code range is noisy but unbiased, the carrier's CHANGE between epochs is exact to
    millimetres. Blend: pr_s = pr/n + (n-1)/n (pr_s_prev + delta_carrier), n growing to M.

    One correction the textbook omits and this radio needs: the LO synthesizer sits a few Hz off
    its nominal, so the carrier's rate differs from the code's by a constant that is COMMON to
    every satellite (-3.1 m/s here, measured: the same on all seven). Harmless while every
    filter has the same age; when one satellite restarts it becomes a differential error of
    tens of metres (28 m scatter, measured). So the median code-minus-carrier rate across the
    satellites is removed from the carrier increment each epoch: what is left is per-satellite
    noise, which is what the filter is for. A jump of more than slip_m between the code range
    and the corrected prediction is a cycle slip or a re-assignment: that satellite restarts.
    `state` is per PRN and persists across epochs; entries need t_sv (slid to s_ref),
    carrier_cycles (slid) and prn."""
    t_ref = s_ref / fs
    raw = []
    for e in entries:
        pr = C * (t_ref - e["t_sv"])                          # metres, with the receiver clock in it
        phase_m = -LAMBDA_L1 * e["carrier_cycles"]           # positive Doppler = closing = range falling
        st = state.get((e.get("sys", "GPS"), e["prn"]))
        d = (pr - st["pr_raw"]) - (phase_m - st["phase_m"]) if st is not None else None
        raw.append((e, pr, phase_m, st, d))
    ds = [d for _, _, _, _, d in raw if d is not None and abs(d) < 3 * slip_m]
    common = float(np.median(ds)) if len(ds) >= 3 else 0.0
    out = []
    for e, pr, phase_m, st, d in raw:
        if st is not None and d is not None and abs(d - common) < slip_m:
            n = min(st["n"] + 1, M)
            pr_s = pr / n + (n - 1) / n * (st["pr_s"] + phase_m - st["phase_m"] + common)
        else:
            n, pr_s = 1, pr
        state[(e.get("sys", "GPS"), e["prn"])] = {"n": n, "pr_s": pr_s, "phase_m": phase_m, "pr_raw": pr}
        out.append(dict(e, t_sv=t_ref - pr_s / C, smoothed=n))
    return out


def fix_from_channels(channels, fs, iono=None, hatch=None, hatch_m=100, mask_deg=0.0):
    """channels: list of {prn, eph, anchor:(epochs, tow), obs:{epochs, code_phase, sample_abs,
    carrier_hz}}. The whole path: transmit times, common instant, solve."""
    entries = []
    for c in channels:
        if c.get("eph") is None or c.get("anchor") is None or c.get("obs") is None:
            continue
        t_sv = transmit_time_sv(c["obs"], c["anchor"], fs)
        o = c["obs"]
        entries.append(dict(prn=c["prn"], eph=c["eph"], t_sv=t_sv, carrier_hz=o.get("carrier_hz", 0.0),
                            epoch_sample=o.get("epoch_sample", o.get("sample_abs")), cn0_db=o.get("cn0_db", 0.0),
                            carrier_cycles=o.get("carrier_cycles", 0.0), sys=o.get("sys", c.get("sys", "GPS"))))
    if len(entries) < 4:
        return None
    s_ref, entries = refer_to_common_sample(entries, fs)
    raw = entries
    if hatch is not None and hatch_m > 1 and all("carrier_cycles" in e for e in entries):
        entries = hatch_smooth(entries, hatch, fs, s_ref, M=hatch_m)
    fx = solve(entries, iono=iono, mask_deg=mask_deg)
    fx["smoothed"] = {("E" if e.get("sys") == "GAL" else "") + str(e["prn"]): e.get("smoothed", 0) for e in entries}
    fx["epoch_sample"] = float(s_ref)
    # the RAW observables at this instant, as a RINEX file wants them: the code pseudorange with
    # the SV clock, ionosphere and troposphere still in it (t_rx from the solve is the epoch,
    # so the receiver clock term is ~0), the carrier phase in cycles (increasing with range),
    # the Doppler (positive = closing) and C/N0. rinex.py writes them.
    fx["observables"] = [dict(sys=e.get("sys", "GPS"), prn=int(e["prn"]),
                              pr_m=float(C * (fx["t_rx"] - e["t_sv"])),
                              phase_cyc=float(-e.get("carrier_cycles", 0.0)),
                              doppler_hz=float(e.get("carrier_hz", 0.0)),
                              cn0_db=float(e.get("cn0_db", 0.0) or 0.0)) for e in raw]
    return fx


LOST_THRESHOLD_S = 100e-6       # a dropped buffer is >= one USB transfer (~0.5 ms); the clock is steady to ns


def samples_lost(prev, now, threshold=LOST_THRESHOLD_S, max_dt_s=30.0):
    """The stream-continuity watchdog ("samples == wall x fs, or void"). `prev` and `now` are
    (epoch_sample / fs, clock_offset_s, drift) of two valid fixes, clock_offset_s being the sample
    clock's offset from GPS time (sample / fs - t_rx) and drift its fitted rate (s/s, or None). The
    offset moves smoothly (-800 ppb on an RSPdx: 0.8 us per second); a JUMP of more than the threshold
    between fixes means samples went missing in between - every channel counted code periods over
    a gap that was not there, and its count is wrong against its anchor by exactly the gap. Returns
    the gap in seconds (negative: samples lost), or None."""
    if prev is None or now is None:
        return None
    dt = now[0] - prev[0]
    if dt <= 0 or dt > max_dt_s:
        return None
    predicted = prev[1] + (prev[2] or 0.0) * dt
    jump = now[1] - predicted
    return float(jump) if abs(jump) > threshold else None


class PositionFilter:
    """Kalman filter on the position solution (gnss-sdr's enable_pvt_kf): state [x y z vx vy vz]
    in ECEF, constant-velocity model, the least-squares fix as the measurement with its own
    covariance. vel_sd is the process noise on velocity per sqrt(second): ~0.05 for a fixed
    antenna, ~1 for a car. The filter is reset when a fix is more than `reset_m` from its
    prediction (a jump the model cannot explain: a re-acquisition, a bad epoch)."""

    def __init__(self, vel_sd=0.5, pos_sd=0.05, meas_sd_min=1.0, reset_m=200.0):
        self.q_v, self.q_p = float(vel_sd), float(pos_sd)
        self.meas_min = float(meas_sd_min)
        self.reset_m = float(reset_m)
        self.x = None
        self.P = None
        self.t = None

    def update(self, t, ecef, cov=None):
        z = np.asarray(ecef, float)
        R = np.asarray(cov, float) if cov is not None and np.all(np.isfinite(cov)) else np.eye(3) * 100.0
        R = R + np.eye(3) * self.meas_min ** 2
        if self.x is None or self.t is None or t < self.t:
            self.x = np.concatenate([z, np.zeros(3)])
            self.P = np.diag([100.0] * 3 + [10.0] * 3)
            self.t = t
            return self.x.copy(), False
        dt = max(t - self.t, 1e-3)
        F = np.eye(6)
        F[:3, 3:] = np.eye(3) * dt
        Q = np.zeros((6, 6))
        Q[:3, :3] = np.eye(3) * (self.q_p ** 2 * dt + self.q_v ** 2 * dt ** 3 / 3)
        Q[:3, 3:] = Q[3:, :3] = np.eye(3) * (self.q_v ** 2 * dt ** 2 / 2)
        Q[3:, 3:] = np.eye(3) * (self.q_v ** 2 * dt)
        xp = F @ self.x
        Pp = F @ self.P @ F.T + Q
        H = np.zeros((3, 6))
        H[:, :3] = np.eye(3)
        y = z - H @ xp
        if np.linalg.norm(y) > self.reset_m:
            self.x = np.concatenate([z, np.zeros(3)])
            self.P = np.diag([100.0] * 3 + [10.0] * 3)
            self.t = t
            return self.x.copy(), True
        S = H @ Pp @ H.T + R
        K = Pp @ H.T @ np.linalg.inv(S)
        self.x = xp + K @ y
        self.P = (np.eye(6) - K @ H) @ Pp
        self.t = t
        return self.x.copy(), False
