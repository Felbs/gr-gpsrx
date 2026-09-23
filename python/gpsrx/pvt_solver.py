#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx PVT: observables + navigation messages in, a position out.

Message-only. Every Channel's 'obs' (once a second: its epoch count and the epoch's arrival
sample) and every Nav Decoder's 'nav' (ephemeris and the timing anchor) come in; whenever
four or more channels have all three - an observable, an ephemeris and an anchor - the block
refers every observable to one receive sample, forms the transmit times from the epoch
counts, and solves. At most once a second. The Klobuchar terms come from any channel that
decoded page 18, or from `iono_file` (an archived broadcast) if none has yet.

'fix' carries the solution and its quality: n, prns, rms_m, pdop, altitude_plausible,
ecef, llh, clock_bias_s, epoch_sample, and a running mean/scatter over the last `average`
fixes. THE POSITION IS IN THE MESSAGE. Nothing in this module prints it; the Sky Panel
shows quality only; write `fix_file` (default none) to keep it, in a directory that is not
in a repository.

The engine is python/gpsrx/pvt.py, tested against a synthetic constellation."""
import json
import os
import threading
import time

import numpy as np
import pmt
from gnuradio import gr

from . import pvt, rinex

EPH_VALID_S = 2.0 * 3600.0      # a broadcast ephemeris is fitted for +-2 h about its toe


class pvt_solver(gr.basic_block):
    def __init__(self, samp_rate=2.048e6, average=15, iono_file="", fix_file="", min_interval_s=1.0, eph_file="",
                 smoothing=100, kf_vel_sd=4.0, rinex_file="", mask_deg=0.0):
        gr.basic_block.__init__(self, name="gpsrx_pvt", in_sig=None, out_sig=None)
        self.fs = float(samp_rate)
        self.average = int(average)
        self.fix_file = fix_file or ""
        # eph_file: a warm start. Every complete ephemeris decoded is kept there by PRN, and a
        # channel that has a timing anchor but has not finished its own decode borrows the stored
        # orbit if the anchor's TOW is within EPH_VALID_S of its toe (the broadcast ephemeris is
        # good for hours; the same satellites are overhead for the whole drive). The timing is
        # always this capture's own. Measured need: on a 90 s capture from a moving car, seven
        # satellites acquired strongly and only three finished decoding - fades break subframes.
        self.eph_file = eph_file or ""
        self.eph_store = {}
        if self.eph_file and os.path.exists(self.eph_file):
            try:
                raw = json.load(open(self.eph_file))
                self.eph_store = {(k if k[:1] in "GE" else "G" + k): v for k, v in raw.items()}   # old files: bare GPS PRNs
            except (ValueError, OSError):
                self.eph_store = {}
        self.n_borrowed = 0
        # carrier smoothing (Hatch) of the code pseudoranges, per satellite, M epochs; 0 = off
        self.hatch_m = int(smoothing)
        self.hatch = {}
        # Kalman filter on the fixes (constant velocity); 0 = off
        self.kf = pvt.PositionFilter(vel_sd=kf_vel_sd) if kf_vel_sd > 0 else None
        self.min_interval = float(min_interval_s)
        self.message_port_register_in(pmt.intern("obs"))
        self.message_port_register_in(pmt.intern("nav"))
        self.message_port_register_in(pmt.intern("status"))
        self.set_msg_handler(pmt.intern("obs"), self.on_obs)
        self.set_msg_handler(pmt.intern("nav"), self.on_nav)
        self.set_msg_handler(pmt.intern("status"), self.on_status)
        # an observable older than this (stream time) than the newest is a dead channel's:
        # a lost satellite's last report poisoned every later fix on a field capture
        # (rms 3 m -> 50 m after one of seven channels dropped, measured 9/22)
        self.max_age_s = 2.5
        self.message_port_register_out(pmt.intern("fix"))
        self._lock = threading.Lock()
        self.obs = {}                       # slot -> latest observable (with prn)
        self.obs_prev = {}                  # slot -> the one before
        self.nav = {}                       # (slot, prn) -> {anchor, eph}
        self.iono, self.iono_src = None, "none"
        if iono_file and os.path.exists(iono_file):
            d = json.load(open(iono_file))
            self.iono, self.iono_src = {"a": d["iono_a"], "b": d["iono_b"]}, "archived"
        self.fixes = []
        self.history = []
        self.n_fixes = 0
        self._last = -1.0
        self._pending = None
        # the timing product: (sample, GPS time) pairs from recent fixes -> clock drift, next PPS
        self._time_hist = []
        self.n_resync = 0                   # stream discontinuities (samples lost) detected and recovered from
        self._clock_prev = None             # (epoch_sample, clock_offset_s, drift) of the last valid fix
        # RINEX 3 output (rinex.py): the observables of every valid fix to <rinex_file>, the decoded
        # ephemerides to the same name with .nav - for RTKLIB and anyone who wants to check our numbers
        self.rinex_file = rinex_file or ""
        self.rinex = rinex.ObsWriter(self.rinex_file) if self.rinex_file else None
        self._rinex_eph = {}                # skey -> the latest complete ephemeris seen (GPS and Galileo)
        self.mask_deg = float(mask_deg)     # elevation mask (0 = weighting only)

    def _samples_lost(self, fx):
        """The stream-continuity watchdog (pvt.samples_lost): did the sample clock's offset from GPS
        time jump since the last valid fix? Keeps (sample time, offset, fitted drift) per fix."""
        s_ref, t_rx = float(fx["epoch_sample"]), float(fx["t_rx"])
        drift = None
        if len(self._time_hist) >= 5:
            S = np.array([p[0] for p in self._time_hist]) / self.fs
            Tg = np.array([p[1] for p in self._time_hist])
            drift = float(np.polyfit(S - S[0], (S - S[0]) - (Tg - Tg[0]), 1)[0])
        now = (s_ref / self.fs, s_ref / self.fs - t_rx, drift)
        prev, self._clock_prev = self._clock_prev, now
        return pvt.samples_lost(prev, now)

    def on_nav(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict):
            return
        with self._lock:
            key = (int(d["slot"]), int(d["prn"]))
            ent = self.nav.setdefault(key, {})
            ent["sys"] = str(d.get("system", "GPS"))
            if d["anchor"][1] >= 0:
                ent["anchor"] = (int(d["anchor"][0]), float(d["anchor"][1]))
            skey = ("E" if ent["sys"] == "GAL" else "G") + str(int(d["prn"]))
            if d.get("complete") and "eph" in d:
                ent["eph"] = dict(d["eph"], prn=int(d["prn"]))
                ent["borrowed"] = False
                if self.rinex is not None:
                    e = dict(ent["eph"], sys=ent["sys"])
                    if ent["sys"] == "GAL" and d.get("wn") is not None:
                        e["WN"] = int(d["wn"])
                    self._rinex_eph[skey] = e
                    try:
                        rinex.write_nav(os.path.splitext(self.rinex_file)[0] + ".nav", list(self._rinex_eph.values()))
                    except OSError:
                        pass
                if self.eph_file:
                    self.eph_store[skey] = ent["eph"]
                    try:
                        os.makedirs(os.path.dirname(os.path.abspath(self.eph_file)), exist_ok=True)
                        with open(self.eph_file, "w") as fh:
                            json.dump({str(k): v for k, v in self.eph_store.items()}, fh)
                    except OSError:
                        pass
            elif "eph" not in ent and skey in self.eph_store and "anchor" in ent:
                cand = self.eph_store[skey]
                tow = float(d["anchor"][1])
                if abs(((tow - cand["toe"]) + 302400) % 604800 - 302400) < EPH_VALID_S:
                    ent["eph"] = cand
                    ent["borrowed"] = True
                    self.n_borrowed += 1
            if "iono" in d and ent["sys"] == "GPS":
                self.iono, self.iono_src = {"a": list(d["iono"]["a"]), "b": list(d["iono"]["b"])}, "decoded"

    def _timing(self, fx):
        """GPS time at the fix's receive sample, the sample clock's drift, and the sample index of the
        next whole GPS second - a 1PPS on the sample clock. The solve's t_rx already has the receiver
        clock bias taken out, so sample s_ref IS GPS time t_rx to the solve's precision (a few ns for
        a metre-level fix). Drift from a straight-line fit of (sample/fs - t_rx) over the last 30 s."""
        s_ref, t_rx = float(fx["epoch_sample"]), float(fx["t_rx"])
        wn = None
        for key, ent in self.nav.items():
            if "eph" in ent and "WN" in ent["eph"]:
                wn = int(ent["eph"]["WN"])
                break
        self._time_hist.append((s_ref, t_rx))
        self._time_hist = [p for p in self._time_hist if s_ref - p[0] < 30.0 * self.fs]
        drift_ppb = None
        if len(self._time_hist) >= 5:
            S = np.array([p[0] for p in self._time_hist]) / self.fs
            Tg = np.array([p[1] for p in self._time_hist])
            # (sample time - GPS time) grows at the clock's rate error; unwrap any week roll first
            off = S - S[0] - (Tg - Tg[0])
            drift_ppb = float(np.polyfit(S - S[0], off, 1)[0] * 1e9)
        rate = 1.0 + (drift_ppb or 0.0) * 1e-9              # sample-clock seconds per GPS second
        next_sec = np.ceil(t_rx + 1e-9)
        return dict(gps_tow=t_rx, gps_week_mod1024=wn, clock_offset_s=s_ref / self.fs - t_rx,
                    clock_drift_ppb=drift_ppb, next_pps_sample=float(s_ref + (next_sec - t_rx) * self.fs * rate),
                    next_pps_tow=float(next_sec))

    def on_status(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            slot = int(d["slot"])
            if d.get("what") in ("lost", "idle", "tracking"):
                self.obs.pop(slot, None)
                self.obs_prev.pop(slot, None)                     # its last observable is not a measurement any more
            if d.get("what") in ("tracking", "lost", "idle", "slip"):
                # a new assignment restarts the channel's epoch count at zero: the old TIMING ANCHOR
                # for this slot is dead (new counts against an old anchor put one satellite 24 s -
                # 7e9 m - off on a drive leg, measured). The ephemeris is per satellite and stays.
                # 'slip': the channel's bit grid moved - whole code periods vanished from the stream -
                # so its count against the old anchor is off by those periods: the same treatment.
                for key in [k for k in self.nav if k[0] == slot]:
                    self.nav[key].pop("anchor", None)
                    self.hatch.pop(key[1], None)
                if d.get("what") == "slip":
                    self.n_resync += 1
                    self._time_hist = []
                    self._clock_prev = None

    def on_obs(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            slot = int(d["slot"])
            self.obs_prev[slot] = self.obs.get(slot)         # two deep: the solve picks the one at or before T
            self.obs[slot] = d
            # Once per second of STREAM time (the receiver's clock is the sample counter), and only
            # when every live channel has delivered its observable for that second: the fix then
            # depends on the samples alone, not on which channel's message happened to arrive
            # first (that gave ~1 m of run-to-run difference on the same file, measured). A
            # channel's observable is "for" second T once its epoch is past (T-1) s.
            now = float(d.get("epoch_sample", 0.0)) / self.fs
            T = int(now // self.min_interval) * self.min_interval
            if T > self._last:
                self._pending = T
            if self._pending is None:
                return
            newest = max(float(o.get("epoch_sample", 0.0)) for o in self.obs.values())
            T_s = self._pending * self.fs
            live = {k: o for k, o in self.obs.items() if newest - float(o.get("epoch_sample", 0.0)) <= self.max_age_s * self.fs}
            if any(float(o.get("epoch_sample", 0.0)) < T_s for o in live.values()):
                return                                        # someone has not reached T yet: wait
            chans = []
            for slot, o in live.items():
                # the observable AT OR BEFORE T: this one if it is, else the previous one
                if float(o.get("epoch_sample", 0.0)) > T_s:
                    o = self.obs_prev.get(slot)
                    if o is None or int(o.get("prn", 0)) != int(live[slot].get("prn", 0)):
                        continue
                if newest - float(o.get("epoch_sample", 0.0)) > self.max_age_s * self.fs:
                    continue
                ent = self.nav.get((slot, int(o["prn"])))
                if ent and "eph" in ent and "anchor" in ent:
                    chans.append(dict(prn=int(o["prn"]), eph=ent["eph"], anchor=ent["anchor"], obs=o,
                                      borrowed=ent.get("borrowed", False), sys=ent.get("sys", "GPS")))
            if len(chans) < 4:
                return
            self._last = self._pending
            self._pending = None
            try:
                fx = pvt.fix_from_channels(chans, self.fs, iono=self.iono,
                                           hatch=self.hatch if self.hatch_m > 1 else None, hatch_m=self.hatch_m,
                                           mask_deg=self.mask_deg)
            except Exception as e:           # a bad ephemeris must not take the flowgraph down
                fx = None
                err = repr(e)
            if fx is None:
                self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(dict(ok=False, n=len(chans))))
                return
            self.n_fixes += 1
            lost_s = self._samples_lost(fx) if fx["valid"] else None
            if lost_s is not None:
                # THE STREAM LOST SAMPLES (a SoapySDR overflow, a dropped buffer): every channel went
                # on counting code periods against samples that never arrived, so every count is
                # now wrong against its anchor by the same amount, the solve's receiver clock jumped
                # by that amount - and nothing else in the receiver can tell. The counting observable
                # is only as good as the stream's continuity ("samples == wall x fs, or void"). Drop
                # every anchor and filter: the nav decoders re-anchor on their next subframe / page.
                for ent in self.nav.values():
                    ent.pop("anchor", None)
                self.hatch.clear()
                self._time_hist = []
                if self.kf is not None:
                    self.kf = pvt.PositionFilter(vel_sd=self.kf.q_v)
                self.fixes = []
                self.n_resync += 1
                self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(dict(
                    ok=True, valid=False, resync=True, samples_lost_s=float(lost_s), n=fx["n"], prns=fx["prns"],
                    rms_m=fx["rms_m"], pdop=fx["pdop"], altitude_plausible=fx["altitude_plausible"],
                    residuals_m=fx["residuals_m"], ecef=fx["ecef"], llh=list(fx["llh"]), clock_bias_s=fx["clock_bias_s"],
                    epoch_sample=fx["epoch_sample"], iono=self.iono_src, count=self.n_fixes, raim=fx["raim"],
                    weights=fx["weights"], azel={str(p): list(v) for p, v in fx["azel"].items()}, isb_s=fx.get("isb_s"))))
                return
            time_out = self._timing(fx) if fx["valid"] else {}
            if self.rinex is not None and fx["valid"] and time_out.get("gps_week_mod1024") is not None:
                try:
                    self.rinex.add_epoch(rinex.full_week(time_out["gps_week_mod1024"]), fx["t_rx"], fx.get("observables", []))
                except OSError:
                    pass
            kf_out = None
            if self.kf is not None and fx["valid"]:
                xk, reset = self.kf.update(fx["epoch_sample"] / self.fs, fx["ecef"], np.array(fx.get("cov_ecef")))
                kf_out = dict(ecef_kf=xk[:3].tolist(), llh_kf=list(pvt.ecef_to_llh(xk[:3])),
                              vel_ecef=xk[3:].tolist(), speed_mps=float(np.linalg.norm(xk[3:])), kf_reset=bool(reset))
            if fx["valid"]:
                self.fixes.append(fx["ecef"])
                self.fixes = self.fixes[-self.average:]
            P = np.array(self.fixes) if self.fixes else None
            out = dict(ok=True, n=fx["n"], prns=fx["prns"], rms_m=fx["rms_m"], pdop=fx["pdop"],
                       altitude_plausible=fx["altitude_plausible"], valid=fx["valid"], residuals_m=fx["residuals_m"],
                       ecef=fx["ecef"], llh=list(fx["llh"]), clock_bias_s=fx["clock_bias_s"],
                       epoch_sample=fx["epoch_sample"], iono=self.iono_src, count=self.n_fixes,
                       raim=fx["raim"], weights=fx["weights"], smoothed=fx.get("smoothed", {}), **(kf_out or {}), **time_out,
                       borrowed=[c["prn"] for c in chans if c.get("borrowed")],
                       azel={str(p): list(v) for p, v in fx["azel"].items()},     # the sky, for the panel
                       isb_s=fx.get("isb_s"))
            if P is not None and len(P) >= 3:
                mean = P.mean(axis=0)
                out["ecef_mean"] = mean.tolist()
                out["llh_mean"] = list(pvt.ecef_to_llh(mean))
                out["averaged"] = int(len(P))
                out["scatter_m"] = float(np.sqrt(np.mean(np.sum((P - mean) ** 2, axis=1))))
            self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(out))
            if self.fix_file:
                # the latest fix, plus every fix so far (quality and ECEF): a run's history is
                # what tells a good stretch from a bad one afterwards
                self.history.append(dict(count=self.n_fixes, n=fx["n"], prns=fx["prns"], rms_m=fx["rms_m"], excluded=fx["raim"]["excluded"],
                                         pdop=fx["pdop"], ecef=fx["ecef"], t_stream_s=fx["epoch_sample"] / self.fs,
                                         ecef_kf=(kf_out or {}).get("ecef_kf"), speed_mps=(kf_out or {}).get("speed_mps"),
                                         altitude_plausible=fx["altitude_plausible"], valid=fx["valid"],
                                         clock_offset_s=time_out.get("clock_offset_s"), resyncs=self.n_resync))
                os.makedirs(os.path.dirname(os.path.abspath(self.fix_file)), exist_ok=True)
                with open(self.fix_file, "w") as fh:
                    json.dump(dict(out, history=self.history), fh, indent=1)
