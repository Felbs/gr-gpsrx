#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""The C++ Channel against the Python Channel: the same synthetic sky through both, the same
observables out (epoch count identical, code phase and Doppler within the loops' noise), the
same tags and messages - and GATE 1: eight C++ channels in the scheduler beat real time."""
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

    def on_msg(self, m):
        self.got.append(pmt.to_python(m))


def run_channels(make, secs, birds, obs_every=200):
    x, truths = synth.sky(FS, secs, birds)
    tb = gr.top_block()
    src = blocks.vector_source_c(x.tolist(), False)
    chans, sinks, obs, stat = [], [], [], []
    for i, b in enumerate(birds):
        ch = make(i, obs_every)
        snk = blocks.vector_sink_c()
        ob, st = msg_sink(), msg_sink()
        tb.connect(src, ch, snk)
        tb.msg_connect(ch, "obs", ob, "in")
        tb.msg_connect(ch, "status", st, "in")
        chans.append(ch)
        sinks.append(snk)
        obs.append(ob)
        stat.append(st)
    tb.keep = [src] + chans + sinks + obs + stat
    # assign BEFORE start: messages queue until the scheduler runs, and an idle C++ channel eats
    # two seconds of samples faster than a sleep can hand it the assignment
    for i, b in enumerate(birds):
        chans[i].to_basic_block()._post(pmt.intern("assign"), pmt.to_pmt(
            dict(slot=i, prn=b["prn"], doppler_hz=b["doppler_hz"] + 20.0, sample=float(b["code_phase_samples"]))))
    t0 = time.time()
    tb.start()
    tb.wait()
    return time.time() - t0, [s.data() for s in sinks], [o.got for o in obs], [s.got for s in stat], truths


class qa_channel_cc(gr_unittest.TestCase):

    def test_001_cpp_matches_python_channel(self):
        secs = 2.0
        birds = BIRDS[:4]
        _, out_py, obs_py, st_py, _ = run_channels(
            lambda i, e: gpsrx.channel(FS, slot=i, obs_every_ms=e), secs, birds)
        _, out_cc, obs_cc, st_cc, _ = run_channels(
            lambda i, e: gpsrx.channel_cc(FS, i, 18.0, 2.0, e), secs, birds)
        for i in range(len(birds)):
            self.assertEqual([d["what"] for d in st_py[i]], [d["what"] for d in st_cc[i]])
            self.assertEqual(len(obs_py[i]), len(obs_cc[i]))
            for a, b in zip(obs_py[i], obs_cc[i]):
                self.assertEqual(a["epochs"], b["epochs"])
                self.assertEqual(a["prn"], b["prn"])
                self.assertLess(abs(a["carrier_hz"] - b["carrier_hz"]), 2.0, (a, b))       # Hz
                self.assertLess(abs(a["code_phase"] - b["code_phase"]), 0.02, (a, b))      # chips
                self.assertLess(abs(a["epoch_sample"] - b["epoch_sample"]), 0.05, (a, b))  # samples
                self.assertGreater(b["lock"], 0.85, b)
            # the prompt streams agree in sign (the data bits) once both loops have settled
            n = min(len(out_py[i]), len(out_cc[i]))
            p, c = np.array(out_py[i][300:n]), np.array(out_cc[i][300:n])
            agree = np.mean(np.sign(p.real) == np.sign(c.real))
            self.assertGreater(agree, 0.99, agree)

    def test_002_gate_1_eight_cpp_channels_beat_real_time(self):
        secs = 4.0
        wall, out, obs, _, _ = run_channels(lambda i, e: gpsrx.channel_cc(FS, i, 18.0, 2.0, e), secs, BIRDS)
        for i in range(len(BIRDS)):
            self.assertGreater(obs[i][-1]["lock"], 0.85, obs[i][-1])
        rate = secs / wall
        print(f"GATE 1 (C++): 8 channels x {secs:.0f} s in the GNU Radio scheduler: {wall:.2f} s wall = {rate:.1f}x real time")
        self.assertGreater(rate, 1.0)


if __name__ == "__main__":
    gr_unittest.run(qa_channel_cc)
