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

from . import pvt

EPH_VALID_S = 2.0 * 3600.0      # a broadcast ephemeris is fitted for +-2 h about its toe


class pvt_solver(gr.basic_block):
    def __init__(self, samp_rate=2.048e6, average=15, iono_file="", fix_file="", min_interval_s=1.0, eph_file="",
                 smoothing=100, kf_vel_sd=4.0):
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
                self.eph_store = {int(k): v for k, v in json.load(open(self.eph_file)).items()}
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
        self.nav = {}                       # (slot, prn) -> {anchor, eph}
        self.iono, self.iono_src = None, "none"
        if iono_file and os.path.exists(iono_file):
            d = json.load(open(iono_file))
            self.iono, self.iono_src = {"a": d["iono_a"], "b": d["iono_b"]}, "archived"
        self.fixes = []
        self.history = []
        self.n_fixes = 0
        self._last = 0.0

    def on_nav(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict):
            return
        with self._lock:
            key = (int(d["slot"]), int(d["prn"]))
            ent = self.nav.setdefault(key, {})
            ent["anchor"] = (int(d["anchor"][0]), float(d["anchor"][1]))
            if d.get("complete") and "eph" in d:
                ent["eph"] = dict(d["eph"], prn=int(d["prn"]))
                ent["borrowed"] = False
                if self.eph_file:
                    self.eph_store[int(d["prn"])] = ent["eph"]
                    try:
                        os.makedirs(os.path.dirname(os.path.abspath(self.eph_file)), exist_ok=True)
                        with open(self.eph_file, "w") as fh:
                            json.dump({str(k): v for k, v in self.eph_store.items()}, fh)
                    except OSError:
                        pass
            elif "eph" not in ent and int(d["prn"]) in self.eph_store:
                cand = self.eph_store[int(d["prn"])]
                tow = float(d["anchor"][1])
                if abs(((tow - cand["toe"]) + 302400) % 604800 - 302400) < EPH_VALID_S:
                    ent["eph"] = cand
                    ent["borrowed"] = True
                    self.n_borrowed += 1
            if "iono" in d:
                self.iono, self.iono_src = {"a": list(d["iono"]["a"]), "b": list(d["iono"]["b"])}, "decoded"

    def on_status(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            slot = int(d["slot"])
            if d.get("what") in ("lost", "idle", "tracking"):
                self.obs.pop(slot, None)                     # its last observable is not a measurement any more
            if d.get("what") in ("tracking", "lost", "idle"):
                # a new assignment restarts the channel's epoch count at zero: the old TIMING ANCHOR
                # for this slot is dead (new counts against an old anchor put one satellite 24 s -
                # 7e9 m - off on a drive leg, measured). The ephemeris is per satellite and stays.
                for key in [k for k in self.nav if k[0] == slot]:
                    self.nav[key].pop("anchor", None)
                    self.hatch.pop(key[1], None)

    def on_obs(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            self.obs[int(d["slot"])] = d
            # once per `min_interval` of STREAM time: the receiver's clock is the sample counter
            # (a wall clock would give one fix per wall second, five per second of capture in replay)
            now = float(d.get("epoch_sample", 0.0)) / self.fs
            if now - self._last < self.min_interval:
                return
            newest = max(float(o.get("epoch_sample", 0.0)) for o in self.obs.values())
            chans = []
            for slot, o in self.obs.items():
                if newest - float(o.get("epoch_sample", 0.0)) > self.max_age_s * self.fs:
                    continue
                ent = self.nav.get((slot, int(o["prn"])))
                if ent and "eph" in ent and "anchor" in ent:
                    chans.append(dict(prn=int(o["prn"]), eph=ent["eph"], anchor=ent["anchor"], obs=o,
                                      borrowed=ent.get("borrowed", False)))
            if len(chans) < 4:
                return
            self._last = now
            try:
                fx = pvt.fix_from_channels(chans, self.fs, iono=self.iono,
                                           hatch=self.hatch if self.hatch_m > 1 else None, hatch_m=self.hatch_m)
            except Exception as e:           # a bad ephemeris must not take the flowgraph down
                fx = None
                err = repr(e)
            if fx is None:
                self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(dict(ok=False, n=len(chans))))
                return
            self.n_fixes += 1
            kf_out = None
            if self.kf is not None and fx["altitude_plausible"] and fx["raim"]["pass"]:
                xk, reset = self.kf.update(fx["epoch_sample"] / self.fs, fx["ecef"], np.array(fx.get("cov_ecef")))
                kf_out = dict(ecef_kf=xk[:3].tolist(), llh_kf=list(pvt.ecef_to_llh(xk[:3])),
                              vel_ecef=xk[3:].tolist(), speed_mps=float(np.linalg.norm(xk[3:])), kf_reset=bool(reset))
            if fx["altitude_plausible"]:
                self.fixes.append(fx["ecef"])
                self.fixes = self.fixes[-self.average:]
            P = np.array(self.fixes) if self.fixes else None
            out = dict(ok=True, n=fx["n"], prns=fx["prns"], rms_m=fx["rms_m"], pdop=fx["pdop"],
                       altitude_plausible=fx["altitude_plausible"], residuals_m=fx["residuals_m"],
                       ecef=fx["ecef"], llh=list(fx["llh"]), clock_bias_s=fx["clock_bias_s"],
                       epoch_sample=fx["epoch_sample"], iono=self.iono_src, count=self.n_fixes,
                       raim=fx["raim"], weights=fx["weights"], smoothed=fx.get("smoothed", {}), **(kf_out or {}),
                       borrowed=[c["prn"] for c in chans if c.get("borrowed")],
                       azel={str(p): list(v) for p, v in fx["azel"].items()})     # the sky, for the panel
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
                                         altitude_plausible=fx["altitude_plausible"]))
                os.makedirs(os.path.dirname(os.path.abspath(self.fix_file)), exist_ok=True)
                with open(self.fix_file, "w") as fh:
                    json.dump(dict(out, history=self.history), fh, indent=1)
