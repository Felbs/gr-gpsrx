#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""The whole receiver in a flowgraph on a synthetic sky: Acquisition finds every satellite and
assigns it, every Channel locks, every Nav Decoder anchors on the real navigation message the
synthetic satellites carry, and PVT is fed. (The synthetic sky has no consistent geometry, so
the solve itself is tested against a synthetic constellation in tests/test_pvt.py.)"""
import os
import sys

import numpy as np
import pmt
from gnuradio import blocks, gr, gr_unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "examples"))
try:
    from gnuradio import gpsrx  # noqa: F401
except ImportError:
    import _devpath  # noqa: E402,F401
    from gnuradio import gpsrx  # noqa: F401
from gnuradio.gpsrx import navgen, synth  # noqa: E402

FS = 2.048e6


class msg_sink(gr.basic_block):
    def __init__(self):
        gr.basic_block.__init__(self, name="msg_sink", in_sig=None, out_sig=None)
        self.got = []
        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self.on_msg)      # by name: the gateway cannot find a lambda

    def on_msg(self, m):
        self.got.append(pmt.to_python(m))


class qa_receiver(gr_unittest.TestCase):

    def test_001_sky_to_anchored_channels(self):
        secs = 22.0
        bits = np.where(navgen.frame_bits(navgen.EXAMPLE_EPH, tow_count0=50400, n_subframes=4) > 0, 1.0, -1.0)
        birds = [dict(prn=p, doppler_hz=d, code_phase_samples=c, cn0_dbhz=45.0, bits=bits) for p, d, c in zip(
            (3, 7, 12, 25), (1234.0, -2510.0, 90.0, 4020.0), (700, 1500, 10, 2000))]
        x, _ = synth.sky(FS, secs, birds)
        tb = gr.top_block()
        src = blocks.vector_source_c(x.tolist(), False)
        rx = gpsrx.receiver(FS, n_channels=4, interval_s=100.0, n_noncoh=40, hold=True)
        sinks = {p: msg_sink() for p in ("sky", "status", "obs", "nav", "fix")}
        tb.connect(src, rx)
        for p, s in sinks.items():
            tb.msg_connect(rx, p, s, "in")
        tb.keep = [src, rx] + list(sinks.values())
        tb.run()
        sky = sinks["sky"].got
        self.assertEqual(len(sky), 1)
        self.assertEqual(sorted(b["prn"] for b in sky[0]["birds"]), [3, 7, 12, 25])
        tracking = {d["prn"] for d in sinks["status"].got if d.get("what") == "tracking"}
        self.assertEqual(tracking, {3, 7, 12, 25})
        self.assertFalse([d for d in sinks["status"].got if d.get("what") == "lost"])
        # every channel locked and reported an absolute observable on the flowgraph's clock
        last = {}
        for o in sinks["obs"].got:
            last[o["prn"]] = o
        self.assertEqual(sorted(last), [3, 7, 12, 25])
        for o in last.values():
            self.assertGreater(o["lock"], 0.85, o)
            self.assertGreater(o["epoch_sample"], 1.0 * FS)                     # seconds into the stream, absolute
        # every decoder anchored: TOW 302400 + 6k on a period 6000 k after the first
        by_prn = {}
        for d in sinks["nav"].got:
            by_prn.setdefault(d["prn"], []).append(d)
        self.assertEqual(sorted(by_prn), [3, 7, 12, 25])
        for prn, msgs in by_prn.items():
            anchors = [m["anchor"] for m in msgs]
            self.assertGreaterEqual(len(anchors), 2, prn)
            self.assertTrue(all(np.diff([a[0] for a in anchors]) == 6000), (prn, anchors))
            self.assertTrue(all(abs(np.diff([a[1] for a in anchors]) - 6.0) < 1e-9), (prn, anchors))
            self.assertTrue(all(abs((a[1] - 302400.0) % 6.0) < 1e-9 for a in anchors), (prn, anchors))


if __name__ == "__main__":
    gr_unittest.run(qa_receiver)
