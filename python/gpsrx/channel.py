#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx Channel: one satellite's tracking loop as a GNU Radio block.

Complex baseband in; per code period (1 ms) the prompt correlation out as a stream item
(complex, 1 kHz: the navigation bits are its sign), and once a second an 'obs' message with
what the solver needs (code-period count, code phase, carrier Doppler, the absolute input
sample index of the period boundary, C/N0, lock).

Idle until 'assign' says which PRN, Doppler and code phase to track - the Acquisition block
sends that. 'assign' with prn 0 idles the channel again. A Channel Bank pre-allocates twelve
of these because a running flowgraph cannot grow; an idle channel costs one memcpy-free
`return` per work() call.

The engine is python/gpsrx/track.py, tested without GNU Radio. This block only does the
stream bookkeeping: consume exactly one code period at a time (the period length changes as
the DLL tracks the code rate), and stamp the absolute sample index of every boundary, which is
the receiver's only clock."""
import threading

import numpy as np
import pmt
from gnuradio import gr

from .cacode import CODE_LEN, CODE_RATE, L1_HZ
from .track import Channel as Engine

LOST_LOCK = 0.2            # PLL lock indicator below this for LOST_PERIODS periods = lost
LOST_PERIODS = 500


class channel(gr.basic_block):
    def __init__(self, samp_rate=2.048e6, slot=0, pll_bw=18.0, dll_bw=2.0, obs_every_ms=1000, batch_ms=3):
        gr.basic_block.__init__(self, name="gpsrx_channel", in_sig=[np.complex64], out_sig=[np.complex64])
        self.fs = float(samp_rate)
        self.slot = int(slot)
        self.pll_bw, self.dll_bw = float(pll_bw), float(dll_bw)
        self.obs_every = int(obs_every_ms)
        self.eng = None
        self.prn = 0
        self._lock = threading.Lock()
        self._lost_run = 0
        self.n_periods = 0
        self.message_port_register_in(pmt.intern("assign"))
        self.set_msg_handler(pmt.intern("assign"), self.on_assign)
        self.message_port_register_out(pmt.intern("obs"))
        self.message_port_register_out(pmt.intern("status"))
        # The output stream is 1 kHz: with a default-sized output buffer the scheduler called
        # general_work once PER PERIOD (measured: 998 calls for 1 s, 0.13x real time for 8 channels)
        # and the Python round trip cost more than the correlations. Ask for periods in batches -
        # small ones: the upstream buffer (8191 items by default, not settable from Python) holds
        # only ~4 periods, and a batch that cannot be filled stalls the flowgraph.
        self.batch = int(batch_ms)
        self.set_output_multiple(self.batch)
        self.set_min_output_buffer(16 * self.batch)

    # ---- assignment from the Acquisition block ----------------------------------------
    def on_assign(self, msg):
        d = pmt.to_python(msg)
        if not isinstance(d, dict):
            return
        if "slot" in d and int(d["slot"]) != self.slot:
            return                                         # every channel hears every assignment
        with self._lock:
            prn = int(d.get("prn", 0))
            if prn <= 0:
                self.eng, self.prn = None, 0
                self._status("idle")
                return
            # 'sample' = the absolute input sample index at which the code starts (acquisition's
            # code phase, made absolute); the engine wants it relative to the first sample it sees,
            # so remember the absolute index and convert in work()
            self._pending = (prn, float(d.get("doppler_hz", 0.0)), float(d.get("sample", 0)))
            self.eng = None
            self.prn = prn

    def _status(self, what, **kw):
        self.message_port_pub(pmt.intern("status"), pmt.to_pmt(dict(slot=self.slot, prn=self.prn, what=what, **kw)))

    # ---- the stream --------------------------------------------------------------------
    def general_work(self, input_items, output_items):
        x = input_items[0]
        out = output_items[0]
        n_in = len(x)
        with self._lock:
            if self.prn == 0:
                self.consume(0, n_in)                     # idle: eat the input, produce nothing
                return 0
            if self.eng is None:
                prn, dop, sample = self._pending
                start = self.nitems_read(0)
                # the engine's sample 0 is THIS call's first sample; acquisition's code start is
                # `sample` absolute, usually seconds in the past by the time the search finishes.
                # The code repeats every 1023 chips AT ITS DOPPLER-SHIFTED RATE: extrapolating with
                # the nominal period would be 3 chips/s wrong at 5 kHz of Doppler (22 chips over
                # a 7 s search). Wrap with the true period.
                period = self.fs * CODE_LEN / (CODE_RATE * (1.0 + dop / L1_HZ))
                rel = (sample - start) % period
                self.eng = Engine(prn, self.fs, dop, rel, pll_bw=self.pll_bw, dll_bw=self.dll_bw)
                self._t0_abs = start                       # engine sample k == absolute start + k
                self._lost_run = 0
                # the prompt stream carries the assignment as a tag on its first item: the Nav
                # Decoder downstream restarts its period count there, in step with the engine's
                self.add_item_tag(0, self.nitems_written(0), pmt.intern("gpsrx_assign"),
                                  pmt.to_pmt(dict(prn=prn, slot=self.slot)))
                self._status("tracking", doppler_hz=dop)
            eng = self.eng
            consumed = produced = 0
            # With set_output_multiple(batch) the scheduler only accepts whole batches: if the input
            # cannot cover one, consume nothing and wait (forecast asked for enough, but a period
            # can run long while the DLL is pulling in).
            want = min(len(out), self.batch)
            if n_in < want * (self._nominal() + 8):
                self.consume(0, 0)
                return 0
            while produced < want:
                need = eng.samples_needed()
                if consumed + need > n_in:
                    break
                ip, qp = eng.step(x[consumed:consumed + need].astype(np.complex128))
                out[produced] = ip + 1j * qp
                produced += 1
                consumed += need
                self.n_periods += 1
                if eng.s.lock < LOST_LOCK:
                    self._lost_run += 1
                else:
                    self._lost_run = 0
                if self.n_periods % self.obs_every == 0:
                    o = eng.observable()
                    # the engine counts samples from its own start; the solver needs the epoch's
                    # arrival on the flowgraph's one clock, nitems_read: make both absolute
                    o["sample_abs"] = int(self._t0_abs + eng.s.samples_in)    # boundary of the period just done
                    o["epoch_sample"] = float(self._t0_abs + o["epoch_sample"])
                    o["slot"] = self.slot
                    self.message_port_pub(pmt.intern("obs"), pmt.to_pmt(o))
                if self._lost_run >= LOST_PERIODS:
                    self._status("lost", periods=self.n_periods)
                    self.eng, self.prn = None, 0
                    break
        self.consume(0, consumed)
        if produced % self.batch:                          # a short batch (end of input, or lost): pad it
            out[produced:produced + (-produced % self.batch)] = 0
            produced += -produced % self.batch
        return produced

    def _nominal(self):
        return int(round(self.fs * 1e-3))

    def forecast(self, noutput_items, ninputs):
        # a period is ~1 ms of samples, and can run a few samples long while the DLL pulls in
        return [int(noutput_items * (self._nominal() + 8))]
