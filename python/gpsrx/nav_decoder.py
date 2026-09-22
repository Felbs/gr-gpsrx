#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""gpsrx Nav Decoder: one Channel's prompt stream in, its navigation message out.

The prompts arrive at 1 kHz (one per code period); the sign of I is the 50 bit/s data after
bit synchronisation. Words, parity, subframes, ephemeris - and the thing the solver actually
needs from this block: the TIMING ANCHOR, which code period a subframe began on and what TOW
that subframe carries. From then on the Channel's period count IS the satellite's clock.

The stream carries a 'gpsrx_assign' tag on the first prompt of each assignment; the decoder
restarts there so its period index matches the Channel's epoch count. Publishes on 'nav'
after every accepted subframe: {slot, prn, subframe, tow, anchor, complete, eph?, iono?}.

The engine is python/gpsrx/nav.py, tested without GNU Radio."""
import numpy as np
import pmt
from gnuradio import gr

from . import nav, nav_gal


class nav_decoder(gr.sync_block):
    def __init__(self, slot=0):
        gr.sync_block.__init__(self, name="gpsrx_nav_decoder", in_sig=[np.complex64], out_sig=None)
        self.slot = int(slot)
        self.message_port_register_out(pmt.intern("nav"))
        self.dec = None
        self.prn = 0
        self.system = "GPS"
        self.k0 = 0                                   # stream index of the current assignment's period 0
        self._tag = pmt.intern("gpsrx_assign")

    def work(self, input_items, output_items):
        x = input_items[0]
        n = len(x)
        r0 = self.nitems_read(0)
        tags = self.get_tags_in_window(0, 0, n, self._tag)
        cuts = [(t.offset - r0, pmt.to_python(t.value)) for t in tags]
        pos = 0
        for i, (cut, d) in enumerate(cuts + [(n, None)]):
            if self.dec is not None:
                for j in range(pos, cut):
                    kept = self.dec.feed(float(x[j].real), r0 + j - self.k0)
                    if self.system == "GAL":
                        for wt, f in kept:
                            self._publish_gal(wt)
                    else:
                        for sf, tow_next, words, _, _ in kept:
                            self._publish(sf, tow_next)
            if d is not None:                             # a new assignment begins at `cut`
                self.prn = int(d.get("prn", 0))
                self.system = str(d.get("system", "GPS"))
                self.dec = None
                if self.prn:
                    self.dec = nav_gal.INavDecoder(self.prn) if self.system == "GAL" else nav.NavDecoder(self.prn)
                self.k0 = r0 + cut
            pos = cut
        return n

    def _publish_gal(self, wt):
        d = self.dec
        msg = dict(slot=self.slot, prn=self.prn, system="GAL", word=int(wt), n_pages=d.n_crc, crc_fail=d.n_crc_fail,
                   complete=bool(d.complete), n_subframes=d.n_crc, subframe=int(wt), tow=float(d.anchors[-1][1]) if d.anchors else -1.0,
                   anchor=[int(d.anchors[-1][0]), float(d.anchors[-1][1])] if d.anchors else [0, -1.0], rejected=d.n_crc_fail)
        if d.complete:
            msg["eph"] = {k: (float(v) if isinstance(v, (int, float)) else v) for k, v in d.eph.items() if k not in ("sys",)}
            msg["eph"]["sys"] = "GAL"
        if d.ggto:
            msg["ggto"] = {k: float(v) for k, v in d.ggto.items()}
        if d.wn is not None:
            msg["wn"] = int(d.wn)
        self.message_port_pub(pmt.intern("nav"), pmt.to_pmt(msg))

    def _publish(self, sf, tow_next):
        d = self.dec
        msg = dict(slot=self.slot, prn=self.prn, subframe=int(sf), tow=float(tow_next - 6.0),
                   anchor=[int(d.anchors[-1][0]), float(d.anchors[-1][1])], n_subframes=len(d.subframes),
                   rejected=d.n_rejected, complete=bool(d.complete))
        if d.complete:
            msg["eph"] = {k: float(v) for k, v in d.eph.items() if k not in ("iono_a", "iono_b")}
        if "iono_a" in d.eph:
            msg["iono"] = dict(a=[float(v) for v in d.eph["iono_a"]], b=[float(v) for v in d.eph["iono_b"]])
        self.message_port_pub(pmt.intern("nav"), pmt.to_pmt(msg))
