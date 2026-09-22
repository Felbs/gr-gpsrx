#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx Acquisition: finds the satellites and hands each one to a Channel.

Complex baseband in, nothing out on the stream. Every `interval_s` seconds it copies a
snapshot of `snapshot_ms` and searches it on its own thread (the search takes seconds; the
stream must not wait), publishes what it found on 'sky', and assigns every new satellite to a
free Channel slot on 'assign'. Channels report back on 'status' (tracking / lost / idle) so
the slot table stays true.

The assignment carries the absolute input sample at which the code started in the snapshot;
by the time the Channel gets it that sample is seconds old, and the Channel wraps it forward
with the Doppler-shifted code period. (Educational point: acquisition only ever tells you
where the code WAS.)

The engine is python/gpsrx/acquire.py, tested without GNU Radio."""
import threading
import time

import numpy as np
import pmt
from gnuradio import gr

from . import acquire


class acquisition(gr.sync_block):
    def __init__(self, samp_rate=2.048e6, n_slots=8, snapshot_ms=110, interval_s=20.0, threshold=2.5,
                 n_noncoh=100, prns=None, hold=False, settle_s=0.5):
        gr.sync_block.__init__(self, name="gpsrx_acquisition", in_sig=[np.complex64], out_sig=None)
        self.fs = float(samp_rate)
        # hold: stop the stream while a search runs. For REPLAY: a file source runs as fast as its
        # readers allow, and idle channels read fast, so 30 s of capture streamed past during one
        # 7 s search (measured) and were never tracked. Live, the radio paces the stream and the
        # search simply costs its own duration of samples - as it does in every receiver.
        self.hold = bool(hold)
        self.n_slots = int(n_slots)
        self.snap_n = int(round(self.fs * snapshot_ms * 1e-3))
        self.interval = float(interval_s)
        self.threshold = float(threshold)
        self.n_noncoh = int(n_noncoh)
        self.prns = list(prns) if prns else acquire.ALL_PRNS
        self.message_port_register_in(pmt.intern("status"))
        self.set_msg_handler(pmt.intern("status"), self.on_status)
        self.message_port_register_out(pmt.intern("sky"))
        self.message_port_register_out(pmt.intern("assign"))
        self._lock = threading.Lock()
        self.slots = {s: 0 for s in range(self.n_slots)}          # slot -> prn (0 = free)
        self._buf, self._buf_start, self._filled = None, 0, 0
        self._next_at = float(settle_s)                            # stream time (s) of the next snapshot: a tuner is still settling at 0
        self._worker = None
        self.n_searches = 0

    # ---- slot table from the channels' reports --------------------------------------------
    def on_status(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "slot" not in d:
            return
        with self._lock:
            s = int(d["slot"])
            if d.get("what") in ("lost", "idle"):
                self.slots[s] = 0
            elif d.get("what") == "tracking":
                self.slots[s] = int(d.get("prn", 0))

    # ---- the stream: copy a snapshot now and then -------------------------------------------
    def work(self, input_items, output_items):
        x = input_items[0]
        n = len(x)
        start = self.nitems_read(0)
        if self.hold and self._worker is not None and self._worker.is_alive():
            time.sleep(0.02)                                   # replay: let the search finish first
            return 0
        if self._buf is None and start / self.fs >= self._next_at and (self._worker is None or not self._worker.is_alive()):
            self._buf = np.empty(self.snap_n, np.complex64)
            self._buf_start, self._filled = start, 0
        if self._buf is not None:
            take = min(n, self.snap_n - self._filled)
            self._buf[self._filled:self._filled + take] = x[:take]
            self._filled += take
            if self._filled >= self.snap_n:
                snap, s0 = self._buf, self._buf_start
                self._buf = None
                self._next_at = (start + n) / self.fs + self.interval
                self._worker = threading.Thread(target=self._search, args=(snap, s0), daemon=True)
                self._worker.start()
        return n

    def _search(self, snap, s0):
        t = time.time()
        found = acquire.sky(snap.astype(np.complex128), self.fs, threshold=self.threshold,
                            n_noncoh=self.n_noncoh, prns=self.prns)
        self.n_searches += 1
        for r in found:
            r["sample"] = float(s0 + r["code_phase"])          # absolute: the flowgraph's one clock
        self.message_port_pub(pmt.intern("sky"), pmt.to_pmt(dict(
            sample0=int(s0), seconds=float(time.time() - t), birds=found, search=self.n_searches)))
        with self._lock:
            tracked = {p for p in self.slots.values() if p}
            for r in found:
                if r["prn"] in tracked:
                    continue
                free = [s for s, p in self.slots.items() if p == 0]
                if not free:
                    break
                s = free[0]
                self.slots[s] = r["prn"]                        # claimed until the channel says otherwise
                self.message_port_pub(pmt.intern("assign"), pmt.to_pmt(dict(
                    slot=s, prn=int(r["prn"]), doppler_hz=float(r["doppler_hz"]), sample=float(r["sample"]),
                    metric=float(r["metric"]))))
