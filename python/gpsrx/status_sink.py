#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx Status: the receiver's state as text, one table every `every_s` seconds.

Takes 'sky', 'status', 'obs', 'nav' and 'fix' messages from the other blocks and prints what
an operator wants: which satellites, lock and Doppler per channel, subframes decoded, and the
fix's QUALITY - satellites used, residual rms, PDOP, scatter. It never prints the position;
the Sky Panel (Qt) draws the same state. `log_file` appends each table (no position either)."""
import threading
import time

import pmt
from gnuradio import gr


class status_sink(gr.basic_block):
    def __init__(self, every_s=2.0, log_file="", quiet=False):
        gr.basic_block.__init__(self, name="gpsrx_status", in_sig=None, out_sig=None)
        self.every = float(every_s)
        self.log_file = log_file or ""
        self.quiet = bool(quiet)
        for port in ("sky", "status", "obs", "nav", "fix"):
            self.message_port_register_in(pmt.intern(port))
        self.set_msg_handler(pmt.intern("sky"), self.on_sky)
        self.set_msg_handler(pmt.intern("status"), self.on_status)
        self.set_msg_handler(pmt.intern("obs"), self.on_obs)
        self.set_msg_handler(pmt.intern("nav"), self.on_nav)
        self.set_msg_handler(pmt.intern("fix"), self.on_fix)
        self._lock = threading.Lock()
        self.sky, self.chan, self.nav, self.fix = None, {}, {}, None
        self.events = []
        self.t0 = time.time()
        self._stop = threading.Event()
        self._thr = threading.Thread(target=self._loop, daemon=True)
        self._thr.start()

    def _d(self, msg):
        d = pmt.to_python(msg)
        return d if isinstance(d, dict) else None

    def on_sky(self, msg):
        d = self._d(msg)
        if d:
            with self._lock:
                self.sky = d
                self.events.append(f"sky: {len(d['birds'])} satellites in {d['seconds']:.1f} s: " + ", ".join(
                    f"PRN{b['prn']} {b['metric']:.1f}" for b in d["birds"]))

    def on_status(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self._lock:
                s = int(d["slot"])
                c = self.chan.setdefault(s, {})
                c["what"] = d.get("what")
                c["prn"] = int(d.get("prn", 0))
                if d.get("what") == "lost":
                    self.events.append(f"slot {s}: PRN{c['prn']} lost after {d.get('periods', 0) / 1000:.0f} s")
                    c["obs"] = None
                    self.nav.pop(s, None)

    def on_obs(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self._lock:
                self.chan.setdefault(int(d["slot"]), {})["obs"] = d

    def on_nav(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self._lock:
                self.nav[int(d["slot"])] = d
                if d.get("n_subframes") == 1:
                    self.events.append(f"slot {d['slot']}: PRN{d['prn']} first subframe (id {d['subframe']}) - timing anchored")
                if d.get("complete") and self.nav.get(("done", d["slot"])) != d["prn"]:
                    self.nav[("done", d["slot"])] = d["prn"]
                    self.events.append(f"slot {d['slot']}: PRN{d['prn']} ephemeris complete")

    def on_fix(self, msg):
        d = self._d(msg)
        if d:
            with self._lock:
                if d.get("ok") and (self.fix is None or not self.fix.get("ok")):
                    self.events.append(f"FIRST FIX: {d['n']} satellites, rms {d['rms_m']:.1f} m, PDOP {d['pdop']:.1f}")
                self.fix = d

    def table(self):
        with self._lock:
            lines = [f"--- gpsrx {time.time() - self.t0:6.0f} s ---"]
            for s in sorted(k for k in self.chan if isinstance(k, int)):
                c = self.chan[s]
                o, nv = c.get("obs"), self.nav.get(s)
                if c.get("what") != "tracking" or not c.get("prn"):
                    lines.append(f" slot {s}: idle")
                    continue
                txt = f" slot {s}: PRN{c['prn']:2d}"
                if o:
                    txt += f"  lock {o['lock']:.2f}  C/N0 {o.get('cn0_db', 0):4.1f}  Doppler {o['carrier_hz']:+7.1f} Hz  {o['epochs'] / 1000:5.0f} s"
                if nv:
                    txt += f"  subframes {nv['n_subframes']}" + ("  EPH" if nv.get("complete") else "")
                    if nv.get("rejected"):
                        txt += f"  (rejected {nv['rejected']})"
                lines.append(txt)
            f = self.fix
            if f and f.get("ok"):
                txt = (f" FIX #{f['count']}: {f['n']} satellites {f['prns']}, rms {f['rms_m']:.1f} m, PDOP {f['pdop']:.1f}, "
                       f"altitude {'plausible' if f['altitude_plausible'] else 'NOT plausible'}, "
                       f"RAIM {'ok' if f.get('valid', True) else 'FAULT - not believed'}"
                       + (f" (PRN {f['raim']['excluded']} excluded)" if f.get('raim', {}).get('excluded') else "") + f", iono {f['iono']}")
                if "scatter_m" in f:
                    txt += f", {f['averaged']}-fix scatter {f['scatter_m']:.1f} m"
                lines.append(txt)
                if f.get("gps_tow") is not None:
                    lines.append(f" TIME: GPS TOW {f['gps_tow']:.6f} s (week {f.get('gps_week_mod1024')} mod 1024) at sample {f['epoch_sample']:.1f}; "
                                 f"sample clock {(f['clock_drift_ppb'] or 0):+.1f} ppb; next PPS at sample {f['next_pps_sample']:.1f}"
                                 + (" (drift not yet fitted)" if f.get("clock_drift_ppb") is None else ""))
            elif f:
                lines.append(f" solve failed with {f.get('n')} channels")
            ev, self.events = self.events, []
        return "\n".join(lines + [" * " + e for e in ev])

    def _loop(self):
        while not self._stop.wait(self.every):
            t = self.table()
            if not self.quiet:
                print(t, flush=True)
            if self.log_file:
                with open(self.log_file, "a") as fh:
                    fh.write(t + "\n")

    def stop(self):
        self._stop.set()
        return True
