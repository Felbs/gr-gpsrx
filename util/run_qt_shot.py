#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""Runs a grcc-generated Qt flowgraph unattended for N seconds, saves a screenshot of its
window, and stops it cleanly (the radio is closed by the flowgraph, never by killing it).

  python run_qt_shot.py build/grc/gpsrx_radio_qt.py --set antenna="<port>" \
      --seconds 40 --png window.png [--call "src.set_gain(0,'IFGR',40)"]

`--set` values become the flowgraph's GRC Parameters. `--call` runs a method on a block
after construction (device-specific gain elements, for instance).
If the radio is shared, RXTUNE_LOCK (see gr-rxtune) names a site lock module and is honoured."""
import argparse
import contextlib
import importlib.util
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))
import _devpath  # noqa: E402,F401

from PyQt5 import Qt, QtCore  # noqa: E402


def site_lock(owner):
    if not os.environ.get("RXTUNE_LOCK"):
        return contextlib.nullcontext(None)
    from rxtune import lock
    return lock.from_env(owner=owner, priority=60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    ap.add_argument("--call", action="append", default=[], metavar="BLOCK.METHOD(ARGS)")
    ap.add_argument("--seconds", type=float, default=40.0)
    ap.add_argument("--png")
    ap.add_argument("--size", default="1400x950")
    a = ap.parse_args()
    kw = {}
    for item in a.set:
        k, v = item.split("=", 1)
        try:
            kw[k] = int(v)
        except ValueError:
            try:
                kw[k] = float(v)
            except ValueError:
                kw[k] = v

    spec = importlib.util.spec_from_file_location("fg", a.module)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cls = getattr(mod, os.path.splitext(os.path.basename(a.module))[0])

    with site_lock("gpsrx-grc") as lk:
        app = Qt.QApplication(sys.argv[:1])
        tb = cls(**kw)
        for call in a.call:
            eval("tb." + call, {"tb": tb})                     # noqa: S307  (the operator's own command line)
        w, h = (int(x) for x in a.size.split("x"))
        tb.resize(w, h)
        tb.start()
        tb.show()
        t0 = time.time()

        def finish():
            timer.stop()
            if a.png:
                tb.grab().save(a.png)
            tb.stop()
            tb.wait()
            app.quit()

        def tick():
            if lk is not None:
                lk.heartbeat()
                if lk.should_yield():
                    finish()
                    return
            if time.time() - t0 > a.seconds:
                finish()

        timer = QtCore.QTimer()
        timer.timeout.connect(tick)
        timer.start(500)
        app.exec_()
    p = tb.panel if hasattr(tb, "panel") else None
    if p is not None:
        f = p.fix
        if f and f.get("ok"):
            print(f"ran {time.time() - t0:.0f} s: fix from {f['n']} satellites, rms {f['rms_m']:.1f} m, PDOP {f['pdop']:.1f}")
            return 0
        print(f"ran {time.time() - t0:.0f} s: no fix ({len(p.chan)} channels heard from)")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
