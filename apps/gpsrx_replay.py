#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The receiver as a GNU Radio flowgraph, on a capture.

  python gpsrx_replay.py CAPTURE.cs16 [--rate 2.048e6] [--channels 8] [--fix-file lab_local/fix_replay.json]

File source -> Receiver hier block (Acquisition, Channels, Nav Decoders, PVT) -> Status sink.
The Status sink prints the receiver's state every 2 s: satellites, lock, subframes, and the
fix's quality. The position is written only to --fix-file (default: a gitignored directory).
Nothing here prints a coordinate.

Python channels run at ~0.2x real time with eight satellites (design gate 1): a 90 s capture
takes several minutes. The same flowgraph on the canvas is examples/gpsrx_replay.grc."""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))
from gnuradio import blocks, gr  # noqa: E402
try:
    from gnuradio import gpsrx
except ImportError:
    import _devpath  # noqa: E402,F401
    from gnuradio import gpsrx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, default=2.048e6)
    ap.add_argument("--start", type=float, default=0.5)
    ap.add_argument("--secs", type=float, default=0.0, help="0 = whole file")
    ap.add_argument("--channels", type=int, default=8)
    ap.add_argument("--interval", type=float, default=20.0, help="seconds of stream between searches")
    ap.add_argument("--fix-file", default=os.path.join(HERE, "..", "lab_local", "fix_replay.json"))
    ap.add_argument("--iono-file", default=os.path.join(HERE, "..", "lab_local", "iono_terms.json"))
    ap.add_argument("--log", default="")
    ap.add_argument("--eph-file", default="", help="warm start: ephemerides kept across runs (a private file)")
    ap.add_argument("--pll-narrow", type=float, default=15.0, help="PLL bandwidth after bit sync (Hz); 5 static, 15-20 moving; 0 = single stage")
    ap.add_argument("--dll-narrow", type=float, default=0.5)
    ap.add_argument("--coherent-ms", type=int, default=20)
    ap.add_argument("--smoothing", type=int, default=100, help="Hatch carrier smoothing window (epochs); 0 = off")
    ap.add_argument("--kf-vel-sd", type=float, default=4.0, help="Kalman PVT velocity process noise, m/s/sqrt(s): 0.05 static, 1 car, 0 = off")
    a = ap.parse_args()

    tb = gr.top_block()
    src = blocks.file_source(gr.sizeof_short, a.capture, False, int(a.start * a.rate) * 2,
                             int(a.secs * a.rate) * 2 if a.secs else 0)
    to_c = blocks.interleaved_short_to_complex(False, False, 32768.0)
    rx = gpsrx.receiver(a.rate, n_channels=a.channels, interval_s=a.interval,
                        iono_file=a.iono_file, fix_file=a.fix_file, hold=True, eph_file=a.eph_file,
                        pll_bw_narrow=a.pll_narrow, dll_bw_narrow=a.dll_narrow, coherent_ms=a.coherent_ms,
                        smoothing=a.smoothing, kf_vel_sd=a.kf_vel_sd)      # replay: stop time while searching
    st = gpsrx.status_sink(every_s=2.0, log_file=a.log)
    tb.connect(src, to_c, rx)
    for port in ("sky", "status", "obs", "nav", "fix"):
        tb.msg_connect(rx, port, st, port)
    tb.keep = [src, to_c, rx, st]
    t0 = time.time()
    print(f"gpsrx replay: {a.capture}, {a.channels} channels; position -> {a.fix_file} (not printed)", flush=True)
    tb.start()
    tb.wait()
    print(st.table())
    print(f"done in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
