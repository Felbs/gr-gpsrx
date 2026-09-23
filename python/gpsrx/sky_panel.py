#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""Sky Panel (Qt): the receiver at a glance, and never the position.

Left: the sky as a polar plot (north up, zenith at the centre) with each tracked satellite
where the fix puts it, coloured by lock; before the first fix, the satellites sit on a ring
by Doppler. Right: one bar per channel - PRN, lock, Doppler, seconds tracked, subframes,
ephemeris - and the fix's quality line: satellites used, residual rms, PDOP, scatter.

The position is in the 'fix' message and this widget reads it only to place the satellites
(their azimuth and elevation come precomputed in the message). It draws no coordinate, no
map, no number that is one."""
import math
import threading
import time

import numpy as np
import pmt
from gnuradio import gr
from PyQt5 import QtCore, QtGui, QtWidgets


def _lock_colour(lock):
    if lock is None:
        return QtGui.QColor(110, 110, 110)
    if lock >= 0.85:
        return QtGui.QColor(40, 180, 90)
    if lock >= 0.5:
        return QtGui.QColor(215, 165, 40)
    return QtGui.QColor(200, 60, 50)


class _Sky(QtWidgets.QWidget):
    def __init__(self, panel):
        super().__init__()
        self.p = panel
        self.setMinimumSize(260, 260)

    def paintEvent(self, _):
        qp = QtGui.QPainter(self)
        qp.setRenderHint(QtGui.QPainter.Antialiasing)
        w, h = self.width(), self.height()
        qp.fillRect(0, 0, w, h, QtGui.QColor(18, 20, 24))
        r = min(w, h) / 2 - 18
        cx, cy = w / 2, h / 2
        qp.setPen(QtGui.QColor(60, 64, 72))
        for el in (0, 30, 60):
            rr = r * (90 - el) / 90
            qp.drawEllipse(QtCore.QPointF(cx, cy), rr, rr)
        qp.drawLine(int(cx - r), int(cy), int(cx + r), int(cy))
        qp.drawLine(int(cx), int(cy - r), int(cx), int(cy + r))
        qp.setPen(QtGui.QColor(140, 140, 140))
        if self.p.private:
            qp.drawText(int(cx - r), int(cy + r + 14), "sky turned, satellites unnamed (private view)")
        else:
            qp.drawText(int(cx - 4), int(cy - r - 4), "N")
            qp.drawText(int(cx + r + 3), int(cy + 4), "E")
        with self.p.lock:
            chans = dict(self.p.chan)
            azel = dict(self.p.azel)
        for slot, c in chans.items():
            prn = c.get("prn")
            if not prn or c.get("what") != "tracking":
                continue
            o = c.get("obs") or {}
            key = ("E" if o.get("sys") == "GAL" else "") + str(prn)     # the fix keys Galileo as E27
            if key in azel:
                az, el = azel[key]
                az += self.p.rotate
                rr = r * (90 - max(el, 0)) / 90
                x, y = cx + rr * math.sin(math.radians(az)), cy - rr * math.cos(math.radians(az))
            else:                                       # no fix yet: a ring, placed by Doppler
                dop = o.get("carrier_hz", 0.0)
                ang = math.radians(90.0 - dop / 7000.0 * 90.0 + self.p.rotate)
                x, y = cx + r * 0.92 * math.cos(ang), cy - r * 0.92 * math.sin(ang)
            col = _lock_colour(o.get("lock") if o else None)
            qp.setBrush(col)
            qp.setPen(col)
            qp.drawEllipse(QtCore.QPointF(x, y), 7, 7)
            if not self.p.private:
                qp.setPen(QtGui.QColor(230, 230, 230))
                qp.drawText(int(x + 9), int(y + 4), key)


class sky_panel(gr.basic_block, QtWidgets.QWidget):
    """
    in : messages 'sky', 'status', 'obs', 'nav', 'fix'
    """

    def __init__(self, label="GPS receiver", parent=None, private=False):
        gr.basic_block.__init__(self, name="gpsrx_sky_panel", in_sig=None, out_sig=None)
        QtWidgets.QWidget.__init__(self, parent)
        # private: for screenshots. A sky plot with PRN numbers at a known time can be inverted
        # to a rough position (the constellation is public). Hide the numbers and turn the whole
        # sky by an angle drawn at start and never shown: dots with no identity and no bearing
        # cannot be inverted. The picture is otherwise true.
        self.private = bool(private)
        self.rotate = float(np.random.default_rng().uniform(0, 360)) if self.private else 0.0
        self.lock = threading.Lock()
        self.chan, self.nav, self.azel = {}, {}, {}
        self.fix, self.sky, self.t0 = None, None, time.time()
        lay = QtWidgets.QHBoxLayout(self)
        self.skyw = _Sky(self)
        right = QtWidgets.QVBoxLayout()
        self.title = QtWidgets.QLabel(f"<b>{label}</b>")
        self.quality = QtWidgets.QLabel("no fix yet")
        self.quality.setStyleSheet("font-size: 13pt;")
        self.table = QtWidgets.QLabel("")
        self.table.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.table.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.table.setAlignment(QtCore.Qt.AlignTop)
        for wdg in (self.title, self.quality, self.table):
            right.addWidget(wdg)
        right.addStretch(1)
        lay.addWidget(self.skyw, 1)
        lay.addLayout(right, 2)
        for port in ("sky", "status", "obs", "nav", "fix"):
            self.message_port_register_in(pmt.intern(port))
        self.set_msg_handler(pmt.intern("sky"), self.on_sky)
        self.set_msg_handler(pmt.intern("status"), self.on_status)
        self.set_msg_handler(pmt.intern("obs"), self.on_obs)
        self.set_msg_handler(pmt.intern("nav"), self.on_nav)
        self.set_msg_handler(pmt.intern("fix"), self.on_fix)
        self.timer = QtCore.QTimer(self)               # handlers run on GNU Radio threads; paint on Qt's
        self.timer.timeout.connect(self.refresh)
        self.timer.start(500)

    def _d(self, msg):
        d = pmt.to_python(msg)
        return d if isinstance(d, dict) else None

    def on_sky(self, msg):
        d = self._d(msg)
        if d:
            with self.lock:
                self.sky = d

    def on_status(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self.lock:
                c = self.chan.setdefault(int(d["slot"]), {})
                c["what"], c["prn"] = d.get("what"), int(d.get("prn", 0))
                if d.get("what") != "tracking":
                    c["obs"] = None
                    self.nav.pop(int(d["slot"]), None)

    def on_obs(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self.lock:
                self.chan.setdefault(int(d["slot"]), {})["obs"] = d

    def on_nav(self, msg):
        d = self._d(msg)
        if d and "slot" in d:
            with self.lock:
                self.nav[int(d["slot"])] = d

    def on_fix(self, msg):
        d = self._d(msg)
        if d:
            with self.lock:
                self.fix = d
                if d.get("ok") and "azel" in d:
                    self.azel = d["azel"]

    def refresh(self):
        with self.lock:
            chans, nav, fix, sky = dict(self.chan), dict(self.nav), self.fix, self.sky
        lines = []
        if sky:
            lines.append(f"last search: {len(sky['birds'])} satellites in {sky['seconds']:.1f} s")
        for s in sorted(chans):
            c = chans[s]
            if c.get("what") != "tracking" or not c.get("prn"):
                lines.append(f"slot {s}: idle")
                continue
            o, nv = c.get("obs"), nav.get(s)
            sysl = "E" if (o or {}).get("sys") == "GAL" else "G"
            txt = f"slot {s}: {sysl}{c['prn']:02d}" if not self.private else f"slot {s}: {sysl}--"
            if o:
                txt += f"  lock {o['lock']:.2f}  C/N0 {o.get('cn0_db', 0):4.1f}" + (f"  {o['carrier_hz']:+7.1f} Hz" if not self.private else "") \
                    + f"  {o['epochs'] / 1000:4.0f} s" + (f"  S4 {o['s4']:.2f}" if o.get("s4") is not None else "")
            if nv:
                txt += f"  {'pages' if nv.get('system') == 'GAL' else 'subframes'} {nv['n_subframes']:2d}" + ("  ephemeris" if nv.get("complete") else "")
            lines.append(txt)
        self.table.setText("\n".join(lines) if lines else "waiting for the first search")
        if fix and fix.get("ok"):
            q = (f"FIX: {fix['n']} satellites, residual rms {fix['rms_m']:.1f} m, PDOP {fix['pdop']:.1f}"
                 + (f", scatter {fix['scatter_m']:.1f} m over {fix['averaged']}" if "scatter_m" in fix else "")
                 + ("" if fix["altitude_plausible"] else "  (altitude NOT plausible)"))
            self.quality.setText(q)
        elif fix:
            self.quality.setText(f"solve failed ({fix.get('n')} channels)")
        self.skyw.update()
