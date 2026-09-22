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


class pvt_solver(gr.basic_block):
    def __init__(self, samp_rate=2.048e6, average=15, iono_file="", fix_file="", min_interval_s=1.0):
        gr.basic_block.__init__(self, name="gpsrx_pvt", in_sig=None, out_sig=None)
        self.fs = float(samp_rate)
        self.average = int(average)
        self.fix_file = fix_file or ""
        self.min_interval = float(min_interval_s)
        self.message_port_register_in(pmt.intern("obs"))
        self.message_port_register_in(pmt.intern("nav"))
        self.set_msg_handler(pmt.intern("obs"), self.on_obs)
        self.set_msg_handler(pmt.intern("nav"), self.on_nav)
        self.message_port_register_out(pmt.intern("fix"))
        self._lock = threading.Lock()
        self.obs = {}                       # slot -> latest observable (with prn)
        self.nav = {}                       # (slot, prn) -> {anchor, eph}
        self.iono, self.iono_src = None, "none"
        if iono_file and os.path.exists(iono_file):
            d = json.load(open(iono_file))
            self.iono, self.iono_src = {"a": d["iono_a"], "b": d["iono_b"]}, "archived"
        self.fixes = []
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
            if "iono" in d:
                self.iono, self.iono_src = {"a": list(d["iono"]["a"]), "b": list(d["iono"]["b"])}, "decoded"

    def on_obs(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            self.obs[int(d["slot"])] = d
            now = time.time()
            if now - self._last < self.min_interval:
                return
            chans = []
            for slot, o in self.obs.items():
                ent = self.nav.get((slot, int(o["prn"])))
                if ent and "eph" in ent and "anchor" in ent:
                    chans.append(dict(prn=int(o["prn"]), eph=ent["eph"], anchor=ent["anchor"], obs=o))
            if len(chans) < 4:
                return
            self._last = now
            try:
                fx = pvt.fix_from_channels(chans, self.fs, iono=self.iono)
            except Exception as e:           # a bad ephemeris must not take the flowgraph down
                fx = None
                err = repr(e)
            if fx is None:
                self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(dict(ok=False, n=len(chans))))
                return
            self.n_fixes += 1
            if fx["altitude_plausible"]:
                self.fixes.append(fx["ecef"])
                self.fixes = self.fixes[-self.average:]
            P = np.array(self.fixes) if self.fixes else None
            out = dict(ok=True, n=fx["n"], prns=fx["prns"], rms_m=fx["rms_m"], pdop=fx["pdop"],
                       altitude_plausible=fx["altitude_plausible"], residuals_m=fx["residuals_m"],
                       ecef=fx["ecef"], llh=list(fx["llh"]), clock_bias_s=fx["clock_bias_s"],
                       epoch_sample=fx["epoch_sample"], iono=self.iono_src, count=self.n_fixes,
                       azel={str(p): list(v) for p, v in fx["azel"].items()})     # the sky, for the panel
            if P is not None and len(P) >= 3:
                mean = P.mean(axis=0)
                out["ecef_mean"] = mean.tolist()
                out["llh_mean"] = list(pvt.ecef_to_llh(mean))
                out["averaged"] = int(len(P))
                out["scatter_m"] = float(np.sqrt(np.mean(np.sum((P - mean) ** 2, axis=1))))
            self.message_port_pub(pmt.intern("fix"), pmt.to_pmt(out))
            if self.fix_file:
                os.makedirs(os.path.dirname(os.path.abspath(self.fix_file)), exist_ok=True)
                with open(self.fix_file, "w") as fh:
                    json.dump(out, fh, indent=1)
