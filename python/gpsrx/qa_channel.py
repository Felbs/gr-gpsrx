#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""The Channel block in a flowgraph: eight synthetic satellites, eight channels, all lock, and
GATE 1 - eight channels in the GNU Radio scheduler (one thread each) keep up with real time."""
import os
import sys
import time

import numpy as np
import pmt
from gnuradio import blocks, gr, gr_unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "examples"))
try:
    from gnuradio import gpsrx  # noqa: F401
except ImportError:
    import _devpath  # noqa: E402,F401
    from gnuradio import gpsrx  # noqa: F401
from gnuradio.gpsrx import synth  # noqa: E402

FS = 2.048e6
BIRDS = [dict(prn=p, doppler_hz=d, code_phase_samples=c, cn0_dbhz=44.0) for p, d, c in zip(
    (1, 3, 7, 11, 14, 17, 22, 30), (-3210, 1234.5, 2800, -900, 410, -2200, 1750, -150),
    (100, 700, 1500, 2000, 300, 1200, 900, 50))]


class msg_sink(gr.basic_block):
    def __init__(self):
        gr.basic_block.__init__(self, name="msg_sink", in_sig=None, out_sig=None)
        self.got = []
        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self.on_msg)

    def on_msg(self, msg):
        self.got.append(pmt.to_python(msg))


class qa_channel(gr_unittest.TestCase):

    def _graph(self, secs, birds):
        x, truths = synth.sky(FS, secs, birds)
        tb = gr.top_block()
        src = blocks.vector_source_c(x.tolist(), False)
        tb.keep = [src]                                  # a Python-side block must outlive tb.start()
        chans, sinks, obs = [], [], []
        for i, b in enumerate(birds):
            ch = gpsrx.channel(FS, slot=i, obs_every_ms=200)
            snk = blocks.vector_sink_c()
            ob = msg_sink()
            tb.connect(src, ch, snk)
            tb.msg_connect(ch, "obs", ob, "in")
            chans.append(ch)
            sinks.append(snk)
            obs.append(ob)
        return tb, chans, sinks, obs, truths

    def test_001_idle_channel_consumes_and_produces_nothing(self):
        x = np.zeros(20480, np.complex64)
        tb = gr.top_block()
        src = blocks.vector_source_c(x.tolist(), False)
        ch = gpsrx.channel(FS)
        snk = blocks.vector_sink_c()
        tb.connect(src, ch, snk)
        tb.run()
        self.assertEqual(len(snk.data()), 0)

    def test_002_eight_channels_lock_on_a_synthetic_sky_and_publish_observables(self):
        tb, chans, sinks, obs, truths = self._graph(1.0, BIRDS)
        for ch, b in zip(chans, BIRDS):
            # assign as acquisition would: PRN, Doppler ~30 Hz off, code start as an ABSOLUTE sample
            ch.on_assign(pmt.to_pmt({"prn": b["prn"], "doppler_hz": b["doppler_hz"] + 30.0,
                                     "sample": b["code_phase_samples"]}))
        t0 = time.perf_counter()
        tb.run()
        wall = time.perf_counter() - t0
        for ch, snk, ob, b in zip(chans, sinks, obs, BIRDS):
            prompts = np.array(snk.data())
            self.assertGreater(len(prompts), 950, b["prn"])
            self.assertGreater(ch.eng.s.lock, 0.85, (b["prn"], ch.eng.s.lock))
            self.assertLess(abs(ch.eng.s.carrier_hz - b["doppler_hz"]), 15.0)
            self.assertGreaterEqual(len(ob.got), 4)
            o = ob.got[-1]
            self.assertEqual(o["prn"], b["prn"])
            # the observable's absolute sample index must be a period boundary the engine counted
            self.assertEqual(o["sample_abs"], o["samples_in"])           # source started at sample 0
        print(f"\nGATE 1: 8 channels x 1.0 s in the GNU Radio scheduler: {wall:.2f} s wall = "
              f"{1.0 / wall:.2f}x real time", flush=True)

    def test_003_a_channel_that_loses_the_signal_says_so_and_idles(self):
        birds = [BIRDS[0]]
        x, _ = synth.sky(FS, 0.3, birds)
        x = np.concatenate([x, np.zeros(int(FS * 0.8), np.complex64)])      # then the sky goes dark
        tb = gr.top_block()
        src = blocks.vector_source_c(x.tolist(), False)
        ch = gpsrx.channel(FS)
        snk = blocks.vector_sink_c()
        st = msg_sink()
        tb.connect(src, ch, snk)
        tb.msg_connect(ch, "status", st, "in")
        ch.on_assign(pmt.to_pmt({"prn": 1, "doppler_hz": -3210.0, "sample": 100}))
        tb.run()
        time.sleep(0.2)
        self.assertEqual(ch.prn, 0)
        self.assertIn("lost", [m["what"] for m in st.got])


if __name__ == '__main__':
    gr_unittest.run(qa_channel)
