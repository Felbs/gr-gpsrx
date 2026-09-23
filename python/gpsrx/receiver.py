#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx Receiver: the whole receiver as one hierarchical block, for people who want the
position and not the canvas.

  complex baseband in -> Acquisition -> N x (Channel -> Nav Decoder) -> PVT -> 'fix' out
                                                                    \\-> 'status' text (optional)

Every message port of the inner blocks is exposed too ('sky', 'status', 'obs', 'nav', 'fix'),
so the same block can feed a Sky Panel or a message debug. The canvas version, with the
blocks laid out one by one, is examples/gpsrx_replay.grc: same blocks, same wires."""
import pmt
from gnuradio import gr

from .acquisition import acquisition
from .channel import channel
from .nav_decoder import nav_decoder
from .pvt_solver import pvt_solver


class receiver(gr.hier_block2):
    def __init__(self, samp_rate=2.048e6, n_channels=8, interval_s=20.0, threshold=2.5, n_noncoh=100,
                 pll_bw=18.0, dll_bw=2.0, iono_file="", fix_file="", average=15, hold=False, engine="auto",
                 eph_file="", pll_bw_narrow=15.0, dll_bw_narrow=0.5, coherent_ms=20, smoothing=100, kf_vel_sd=4.0, lo_search_hz=0.0, pll_order=3, n_galileo=0,
                 rinex_file=""):
        # engine: "cpp" (channel_cc, 15x real time for eight), "python" (channel, 0.17x: replay
        # only), "auto" = C++ if it was built. Same ports, tags and messages either way.
        try:
            from .gpsrx_python import channel_cc
        except ImportError:
            channel_cc = None
        if engine == "auto":
            engine = "cpp" if channel_cc is not None else "python"
        if engine == "cpp" and channel_cc is None:
            raise RuntimeError("gpsrx: the C++ Channel is not built (see util/build_win.cmd or cmake); engine='python' for replay")
        self.engine = engine
        gr.hier_block2.__init__(self, "gpsrx_receiver",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                                gr.io_signature(0, 0, 0))
        for port in ("sky", "status", "obs", "nav", "fix"):
            self.message_port_register_hier_out(port)
        self.acq = acquisition(samp_rate, n_slots=n_channels, interval_s=interval_s, threshold=threshold,
                               n_noncoh=n_noncoh, hold=hold, lo_search_hz=lo_search_hz)
        # Galileo E1: its own Acquisition (4 ms BOC replicas) owning the slots after the GPS ones,
        # pilot-aided channels (C++ or Python, like the GPS ones), the same Nav Decoder block
        self.acq_gal = None
        if n_galileo:
            self.acq_gal = acquisition(samp_rate, n_slots=n_galileo, interval_s=interval_s, threshold=threshold,
                                       n_noncoh=n_noncoh, hold=hold, lo_search_hz=0.0, system="GAL", slot0=int(n_channels),
                                       settle_s=0.5 + 0.2)
        self.pvt = pvt_solver(samp_rate, average=average, iono_file=iono_file, fix_file=fix_file, eph_file=eph_file,
                              smoothing=smoothing, kf_vel_sd=kf_vel_sd, rinex_file=rinex_file)
        self.chans, self.decs = [], []
        self.connect(self, self.acq)
        if self.acq_gal is not None:
            self.connect(self, self.acq_gal)
            self.msg_connect(self.acq_gal, "sky", self, "sky")
        for s in range(int(n_channels) + int(n_galileo)):
            gal = s >= int(n_channels)
            acq = self.acq_gal if gal else self.acq
            if engine == "cpp":
                ch = channel_cc(float(samp_rate), s, float(pll_bw) if not gal else 12.0, float(dll_bw) if not gal else 1.0, 1000,
                                float(pll_bw_narrow), float(dll_bw_narrow), int(coherent_ms), int(pll_order),
                                "E1" if gal else "L1CA")
            else:
                ch = channel(samp_rate, slot=s, pll_bw=pll_bw if not gal else 12.0, dll_bw=dll_bw if not gal else 1.0,
                             pll_bw_narrow=pll_bw_narrow, dll_bw_narrow=dll_bw_narrow, coherent_ms=coherent_ms,
                             pll_order=pll_order, signal="E1" if gal else "L1CA",      # E1 = pilot-aided
                             obs_every_ms=250 if gal else 1000)      # an observable per second: 4 ms periods
            dec = nav_decoder(slot=s)
            self.connect(self, ch, dec)
            self.msg_connect(acq, "assign", ch, "assign")
            self.msg_connect(ch, "status", acq, "status")
            self.msg_connect(ch, "obs", self.pvt, "obs")
            self.msg_connect(ch, "status", self.pvt, "status")
            self.msg_connect(dec, "nav", self.pvt, "nav")
            self.msg_connect(ch, "status", self, "status")
            self.msg_connect(ch, "obs", self, "obs")
            self.msg_connect(dec, "nav", self, "nav")
            self.chans.append(ch)
            self.decs.append(dec)
        self.msg_connect(self.acq, "sky", self, "sky")
        self.msg_connect(self.pvt, "fix", self, "fix")
