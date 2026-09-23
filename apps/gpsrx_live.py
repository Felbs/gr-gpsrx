#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The receiver live on a SoapySDR radio, headless: for soaks and for boxes without a screen.

  python gpsrx_live.py --driver sdrplay --antenna "Antenna B" --bias-t --seconds 3600 \\
      --fix-file /somewhere/private/live.json --log /somewhere/private/live.log

Radio -> Receiver (C++ channels) -> Status. Prints the status table every 10 s and a one-line
summary at the end (satellites, fixes, rms, scatter, re-acquisitions). The position goes only
to --fix-file. --bias-t writes the device-level SDRplay setting (biasT_ctrl) and turns it off
at exit; other radios take --settings key=value for their own. Honours RXTUNE_LOCK if set."""
import argparse
import contextlib
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))
from gnuradio import gr, soapy  # noqa: E402
try:
    from gnuradio import gpsrx
except ImportError:
    import _devpath  # noqa: E402,F401
    from gnuradio import gpsrx


def site_lock(owner):
    if not os.environ.get("RXTUNE_LOCK"):
        return contextlib.nullcontext(None)
    from rxtune import lock
    return lock.from_env(owner=owner, priority=60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", default="sdrplay")
    ap.add_argument("--dev-args", default="")
    ap.add_argument("--antenna", default="")
    ap.add_argument("--rate", type=float, default=2.048e6)
    ap.add_argument("--gain", type=float, default=40.0)
    ap.add_argument("--agc", type=int, default=1)
    ap.add_argument("--bias-t", action="store_true", help="SDRplay biasT_ctrl on (an active antenna)")
    ap.add_argument("--settings", action="append", default=[], metavar="KEY=VALUE", help="device-level settings")
    ap.add_argument("--seconds", type=float, default=600.0)
    ap.add_argument("--channels", type=int, default=8)
    ap.add_argument("--galileo", type=int, default=0, help="Galileo E1 (pilot-aided) channels; needs --rate 4.096e6")
    ap.add_argument("--interval", type=float, default=20.0)
    ap.add_argument("--fix-file", default=os.path.join(HERE, "..", "lab_local", "live.json"))
    ap.add_argument("--rinex-file", default="", help="write RINEX 3 observations here (+ .nav with the ephemerides); a private file")
    ap.add_argument("--mask-deg", type=float, default=0.0, help="elevation mask (deg): satellites below it are left out when enough remain; 0 = weighting only")
    ap.add_argument("--eph-file", default=os.path.join(HERE, "..", "lab_local", "eph_live.json"))
    ap.add_argument("--iono-file", default=os.path.join(HERE, "..", "lab_local", "iono_terms.json"))
    ap.add_argument("--log", default="")
    ap.add_argument("--pll-narrow", type=float, default=15.0)
    ap.add_argument("--smoothing", type=int, default=100)
    ap.add_argument("--lo-search", type=float, default=0.0, help="RTL-SDR-class radios: wide first search for the LO offset, +-Hz (50000)")
    ap.add_argument("--pll-order", type=int, default=3, help="narrow-stage PLL order: 3 follows a Doppler rate (moving), 2 classic")
    ap.add_argument("--kf-vel-sd", type=float, default=4.0)
    ap.add_argument("--every", type=float, default=10.0)
    a = ap.parse_args()

    with site_lock("gpsrx-live") as lk:
        tb = gr.top_block()
        src = soapy.source(f"driver={a.driver}" + (f",{a.dev_args}" if a.dev_args else ""), "fc32", 1, "", "", [""], [""])
        src.set_sample_rate(0, a.rate)
        src.set_bandwidth(0, 2.5e6 if a.rate <= 2.5e6 else a.rate)   # BOC main lobes sit at +-1.023 MHz
        if a.antenna:
            src.set_antenna(0, a.antenna)
        src.set_frequency(0, 1575.42e6)
        src.set_gain_mode(0, bool(a.agc))
        src.set_gain(0, a.gain)
        src.set_min_output_buffer(int(16 * a.rate * 4e-3))     # whole code periods per work() call
        settings = list(a.settings) + (["biasT_ctrl=true"] if a.bias_t else [])
        for kv in settings:
            k, v = kv.split("=", 1)
            src.write_setting(k, v)
            print(f"setting {k} = {src.read_setting(k) if hasattr(src, 'read_setting') else v}", flush=True)
        rx = gpsrx.receiver(a.rate, n_channels=a.channels, interval_s=a.interval, iono_file=a.iono_file,
                            fix_file=a.fix_file, eph_file=a.eph_file, hold=False, engine="cpp",
                            pll_bw_narrow=a.pll_narrow, smoothing=a.smoothing, kf_vel_sd=a.kf_vel_sd, lo_search_hz=a.lo_search, pll_order=a.pll_order, n_galileo=a.galileo,
                            rinex_file=a.rinex_file, mask_deg=a.mask_deg)
        st = gpsrx.status_sink(every_s=a.every, log_file=a.log)
        tb.connect(src, rx)
        for port in ("sky", "status", "obs", "nav", "fix"):
            tb.msg_connect(rx, port, st, port)
        tb.keep = [src, rx, st]
        t0 = time.time()
        print(f"gpsrx live: {a.driver} {a.antenna or ''} for {a.seconds:.0f} s; position -> {a.fix_file} (not printed)", flush=True)
        tb.start()
        # the capture-integrity law, live: samples == wall x fs, or the stream is not to be trusted.
        # The radio's item count against the wall clock; the lag is constant once the pipeline has
        # filled, so a GROWING deficit is lost samples (a SoapySDR overflow prints 'O' and nothing else)
        # The item count advances in whole scheduler chunks (+-100 ms of jitter against the wall
        # clock, measured), so a loss is a STEP in the deficit's floor, not a single reading: the
        # minimum over the latest reports against the minimum over the first ones.
        n0, w0, deficits = None, None, []
        next_report = t0 + a.every
        try:
            while time.time() - t0 < a.seconds:
                time.sleep(1.0)
                if lk is not None:
                    lk.heartbeat()
                    if lk.should_yield():
                        print("yielding the radio", flush=True)
                        break
                now = time.time()
                if now >= next_report:
                    next_report += a.every
                    n = src.nitems_written(0)
                    if n0 is None and n > 0:
                        n0, w0 = n, now
                    elif n0 is not None:
                        deficits.append(((now - w0) * a.rate - (n - n0)) / a.rate * 1e3)   # ms the wall clock expected and did not get
                        txt = f" stream: {n} samples in {now - t0:.0f} s; wall x fs - samples = {deficits[-1]:+.0f} ms (jitter +-100)"
                        if len(deficits) >= 8:
                            lost = min(deficits[-4:]) - min(deficits[:4])
                            txt += f"; lost so far ~{lost:+.0f} ms" + (" <- SAMPLES LOST" if lost > 100.0 else "")
                        print(txt, flush=True)
        except KeyboardInterrupt:
            pass
        tb.stop()
        tb.wait()
        if a.bias_t:
            try:
                src.write_setting("biasT_ctrl", "false")
            except Exception:
                pass
        f = st.fix
        n_lost = sum(1 for c in st.chan.values() if c.get("what") == "lost")
        print(f"ran {time.time() - t0:.0f} s: " + (f"fix #{f['count']} from {f['n']} satellites, rms {f['rms_m']:.1f} m, "
              f"PDOP {f['pdop']:.1f}, {f.get('averaged', 0)}-fix scatter {f.get('scatter_m', float('nan')):.1f} m"
              if f and f.get("ok") else "no fix"), flush=True)
    return 0 if (st.fix and st.fix.get("ok")) else 1


if __name__ == "__main__":
    sys.exit(main())
