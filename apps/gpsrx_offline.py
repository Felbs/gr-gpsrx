#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The whole receiver, offline, on a capture: acquisition -> channels -> nav decode -> fix.

  python gpsrx_offline.py CAPTURE.cs16 [--secs 90] [--rate 2.048e6] [--out lab_local/fix.json]

Prints the sky, each channel's lock and ephemeris progress, and the fix's QUALITY (satellites,
rms, PDOP, altitude plausible). The position itself goes only to --out, which defaults to a
gitignored directory. Nothing here prints a coordinate.

Acquisition uses numpy-gps's acquire() when NUMPY_GPS_DIR is set (the same FFT search the
design's Acquisition block will carry); tracking, decoding and solving are this project's."""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from gpsrx import nav, pvt  # noqa: E402
from gpsrx.track import Channel  # noqa: E402


def load_cs16(path, fs, t0, secs):
    n = int(secs * fs)
    with open(path, "rb") as fh:
        fh.seek(int(t0 * fs) * 4)
        raw = np.frombuffer(fh.read(n * 4), np.int16).astype(np.float32) / 32768.0
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex128)


def acquire(x, fs, thr):
    sys.path.insert(0, os.environ["NUMPY_GPS_DIR"])
    import measure as gt
    acq = gt.acquire(x[:int(0.4 * fs)], fs, list(range(1, 33)), np.arange(-7000, 7001, 250.0), 400)
    sky = {p: r for p, r in acq.items() if r["metric"] > thr}
    # refine each Doppler to ~10 Hz: the coarse search's 250 Hz bins are wider than the PLL's pull-in
    # (PRN 22 and 30 failed to lock at +-125 Hz, measured). 100 ms of prompts at the coarse Doppler and
    # code phase; the residual tone's frequency is the Doppler error.
    for p, r in sky.items():
        pr = gt.prompts_ms(x, fs, p, r["dopp"], int(r["code_phase"]), 100)
        spec = np.abs(np.fft.fftshift(np.fft.fft(pr ** 2, 4096)))            # squaring removes the data bits
        fax = np.fft.fftshift(np.fft.fftfreq(4096, 1e-3))
        r["dopp"] = float(r["dopp"] + fax[int(np.argmax(spec))] / 2.0)
    return sky


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, default=2.048e6)
    ap.add_argument("--start", type=float, default=0.5)
    ap.add_argument("--secs", type=float, default=90.0)
    ap.add_argument("--thr", type=float, default=3.5)
    ap.add_argument("--max-channels", type=int, default=10)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "lab_local", "fix.json"))
    a = ap.parse_args()
    if not os.environ.get("NUMPY_GPS_DIR"):
        sys.exit("set NUMPY_GPS_DIR (acquisition is borrowed from numpy-gps in this phase)")
    fs = a.rate
    t0 = time.time()
    x = load_cs16(a.capture, fs, a.start, a.secs)
    print(f"loaded {len(x) / fs:.1f} s")
    sky = acquire(x, fs, a.thr)
    print(f"SKY: {len(sky)} satellites above metric {a.thr}: " + ", ".join(
        f"PRN{p} ({r['metric']:.1f}, {r['dopp']:+.0f} Hz)" for p, r in sorted(sky.items(), key=lambda kv: -kv[1]['metric'])))
    birds = sorted(sky.items(), key=lambda kv: -kv[1]["metric"])[:a.max_channels]

    # ---- track + decode every bird, one channel each (sequentially here; blocks run in parallel)
    chans = []
    for prn, r in birds:
        ch = Channel(prn, fs, float(r["dopp"]), float(r["code_phase"]))
        dec = nav.NavDecoder(prn)
        pos, k = 0, 0
        last_obs = None
        t1 = time.time()
        while True:
            n = ch.samples_needed()
            if pos + n > len(x):
                break
            ip, qp = ch.step(x[pos:pos + n])
            pos += n
            dec.feed(ip, k)
            k += 1
            if k % 1000 == 0:                                  # an observable a second, as the block publishes
                ch.obs_hist = getattr(ch, "obs_hist", []) + [ch.observable()]
        obs = ch.observable()
        anchor = dec.anchors[0] if dec.anchors else None
        missing = [key for key in nav.EPH_KEYS if key not in dec.eph]
        print(f"  PRN{prn:2d}: lock {ch.s.lock:.2f}  Doppler {ch.s.carrier_hz:+7.1f} Hz  bit sync {'yes' if dec.bits else 'NO'}  "
              f"subframes {len(dec.subframes)}  ephemeris {'COMPLETE' if dec.complete else 'missing ' + ','.join(missing)}  "
              f"({time.time() - t1:.0f} s)", flush=True)
        chans.append({"prn": prn, "eph": dec.eph if dec.complete else None, "anchor": anchor, "obs": obs,
                      "obs_hist": getattr(ch, "obs_hist", []), "lock": ch.s.lock, "subframes": len(dec.subframes)})

    usable = [c for c in chans if c["eph"] and c["anchor"]]
    print(f"\n{len(usable)} satellites with ephemeris + timing anchor")
    if len(usable) < 4:
        print("fewer than 4: no fix (a longer capture, or more sky, is needed)")
        return 2
    iono, iono_src = None, "not available"
    for c in usable:
        if "iono_a" in c["eph"]:
            iono, iono_src = {"a": c["eph"]["iono_a"], "b": c["eph"]["iono_b"]}, "decoded (subframe 4 page 18)"
    if iono is None:
        # page 18 comes once per 12.5 min; a short capture usually misses it. The terms vary slowly
        # (a day's are fine), so use an archived broadcast if there is one: numpy-gps keeps its last
        # decoded set in lab_local/iono_terms.json, and this receiver writes its own the same way.
        for cand in (os.path.join(HERE, "..", "lab_local", "iono_terms.json"),
                     os.path.join(os.environ.get("NUMPY_GPS_DIR", ""), "lab_local", "iono_terms.json")):
            if os.path.exists(cand):
                d = json.load(open(cand))
                iono, iono_src = {"a": d["iono_a"], "b": d["iono_b"]}, f"archived ({d.get('src', cand)})"
                break
    else:
        with open(os.path.join(HERE, "..", "lab_local", "iono_terms.json"), "w") as fh:
            json.dump({"iono_a": iono["a"], "iono_b": iono["b"], "src": os.path.basename(a.capture)}, fh)
    fx = pvt.fix_from_channels(usable, fs, iono=iono)
    # several instants: each channel kept its observables every second (below); average the fixes
    # the way numpy-gps does over epochs, and report the scatter as the honest precision number
    fixes = []
    n_ep = min(len(c["obs_hist"]) for c in usable)
    for k in range(max(0, n_ep - 15), n_ep):
        ent = [dict(c, obs=c["obs_hist"][k]) for c in usable]
        f = pvt.fix_from_channels(ent, fs, iono=iono)
        if f and f["altitude_plausible"]:
            fixes.append(f)
    if len(fixes) >= 3:
        P = np.array([f["ecef"] for f in fixes])
        mean = P.mean(axis=0)
        scatter = float(np.sqrt(np.mean(np.sum((P - mean) ** 2, axis=1))))
        fx["ecef_mean"] = mean.tolist()
        fx["llh_mean"] = pvt.ecef_to_llh(mean)
        fx["epochs"] = len(fixes)
        fx["scatter_m"] = scatter
        fx["rms_mean_m"] = float(np.mean([f["rms_m"] for f in fixes]))
        print(f"{len(fixes)} epochs: mean rms {fx['rms_mean_m']:.1f} m, epoch scatter {scatter:.1f} m")
    print(f"FIX: {fx['n']} satellites {fx['prns']}, residual rms {fx['rms_m']:.1f} m, PDOP {fx['pdop']:.1f}, "
          f"altitude {'plausible' if fx['altitude_plausible'] else 'NOT plausible'}, "
          f"ionosphere {iono_src}; {time.time() - t0:.0f} s total")
    print("residuals (m): " + " ".join(f"{r:+.1f}" for r in fx["residuals_m"]))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump({"capture": os.path.abspath(a.capture), "fix": fx, "channels": [
            {k: v for k, v in c.items() if k not in ("eph", "obs_hist")} for c in chans]}, fh, indent=1, default=float)
    print(f"position written to {a.out} (not printed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
