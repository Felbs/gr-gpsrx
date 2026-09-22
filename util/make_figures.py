#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The walkthrough's figures, every one from SYNTHETIC satellites (no capture, no radio, no
position): docs/img/walk_*.png. Re-run after changing the engines."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))
from gpsrx import acquire, cacode, nav, navgen, synth  # noqa: E402
from gpsrx.track import Channel  # noqa: E402

OUT = os.path.join(ROOT, "docs", "img")
os.makedirs(OUT, exist_ok=True)
FS = 2.048e6
plt.rcParams.update({"figure.dpi": 110, "font.size": 9, "axes.grid": True, "grid.alpha": 0.3})


def save(name):
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, name))
    plt.close()
    print("wrote", name)


# ---- 1. the C/A code and its autocorrelation: why a correlator can find a signal under the noise
c = cacode.ca_code(7)
ac = np.array([np.dot(c, np.roll(c, k)) for k in range(1023)])
fig, ax = plt.subplots(1, 2, figsize=(9, 2.8))
ax[0].step(np.arange(60), c[:60], where="post")
ax[0].set(title="PRN 7's C/A code, first 60 of 1023 chips", xlabel="chip", ylabel="+1 / -1")
ax[1].plot(np.arange(-511, 512), np.roll(ac, 511))
ax[1].set(title="autocorrelation: 1023 at zero lag, |<= 65| elsewhere", xlabel="lag (chips)")
save("walk_1_code.png")

# ---- 2. acquisition: the Doppler x code-phase map of one satellite in a synthetic sky
birds = [dict(prn=p, doppler_hz=d, code_phase_samples=cp, cn0_dbhz=cn) for p, d, cp, cn in
         ((7, -2510.0, 1500, 44), (3, 1234.0, 700, 45), (12, 90.0, 10, 40), (25, 4020.0, 2000, 38))]
x, truths = synth.sky(FS, 0.3, birds)
xd = x.astype(np.complex128)
dops = np.arange(-7000, 7001, 250.0)
m = acquire.search(xd, FS, prns=[7], dopplers=dops, n_noncoh=20)[7]
n1 = 2048
blocks = xd[:n1 * 20].reshape(20, n1)
t = np.arange(n1) / FS
code_f = np.conj(np.fft.fft(cacode.sampled_code(7, FS, n1)))
grid = np.empty((len(dops), n1))
for i, fd in enumerate(dops):
    BF = np.fft.fft(blocks * np.exp(-2j * np.pi * fd * t)[None, :], axis=1)
    grid[i] = (np.abs(np.fft.ifft(BF * code_f[None, :], axis=1)) ** 2).sum(axis=0)
fig = plt.figure(figsize=(9, 3.4))
ax = fig.add_subplot(1, 2, 1)
ax.imshow(grid, aspect="auto", origin="lower", extent=[0, n1, dops[0], dops[-1]], cmap="magma")
ax.set(title=f"PRN 7 search: peak at {m['code_phase']:.0f} samples, {m['doppler_hz']:+.0f} Hz (metric {m['metric']:.1f})",
       xlabel="code phase (samples)", ylabel="Doppler (Hz)")
ax.grid(False)
ax = fig.add_subplot(1, 2, 2)
ax.plot(grid[int(np.argmax(grid.max(axis=1)))])
ax.set(title="the Doppler row through the peak", xlabel="code phase (samples)", ylabel="correlation power")
save("walk_2_acquisition.png")

# ---- 3. tracking: prompt I/Q through pull-in, then the data bits; both stages
bits = np.where(navgen.frame_bits(navgen.EXAMPLE_EPH, tow_count0=50400, n_subframes=3) > 0, 1.0, -1.0)
x1, tr = synth.satellite(7, FS, 2.5, doppler_hz=987.0, code_phase_samples=333, cn0_dbhz=42, bits=bits, seed=2)
x1 = x1.astype(np.complex128)
ch = Channel(7, FS, 987.0 + 30.0, 333.4)
pos = 0
I, Q, F, L, sync_k = [], [], [], [], None
while pos + ch.samples_needed() <= len(x1):
    n = ch.samples_needed()
    ip, qp = ch.step(x1[pos:pos + n])
    pos += n
    I.append(ip)
    Q.append(qp)
    F.append(ch.s.carrier_hz - tr["doppler_hz"])
    L.append(ch.s.lock)
    if sync_k is None and ch.bit_offset is not None:
        sync_k = len(I)
fig, ax = plt.subplots(3, 1, figsize=(9, 6), sharex=True)
tt = np.arange(len(I)) * 1e-3
ax[0].plot(tt, I, lw=0.6, label="I (prompt)")
ax[0].plot(tt, Q, lw=0.6, label="Q (prompt)")
ax[0].set(title="the prompt correlator: energy moves into I as the PLL locks; then I carries the 50 bit/s data", ylabel="correlation")
ax[0].legend(loc="upper right")
ax[1].plot(tt, F)
ax[1].set(ylabel="Doppler error (Hz)", title="carrier loop: 30 Hz off at handover, pulled in; narrow stage after bit sync")
ax[2].plot(tt, L)
ax[2].set(ylabel="lock (I^2-Q^2)/(I^2+Q^2)", xlabel="time (s)")
for a in ax:
    if sync_k:
        a.axvline(sync_k * 1e-3, color="k", ls="--", lw=0.8)
ax[2].text(sync_k * 1e-3 + 0.02, 0.1, "bit sync -> narrow loops, 20 ms integration", fontsize=8)
save("walk_3_tracking.png")

# ---- 4. bit sync: the flip histogram
bs = nav.BitSync()
for k, ip in enumerate(I[300:]):
    bs.feed(ip, k)
fig, ax = plt.subplots(figsize=(6, 2.6))
ax.bar(range(20), bs.hist)
ax.set(title="bit sync: at which period (mod 20) does the prompt's sign flip?", xlabel="period mod 20", ylabel="flips")
save("walk_4_bitsync.png")

# ---- 5. the discriminators: E/P/L over code error, atan(Q/I) vs atan2 over phase error
err = np.linspace(-1.5, 1.5, 301)
tri = np.clip(1 - np.abs(err), 0, None)
E, P, Lc = np.clip(1 - np.abs(err - 0.5), 0, None), tri, np.clip(1 - np.abs(err + 0.5), 0, None)
fig, ax = plt.subplots(1, 2, figsize=(9, 2.8))
ax[0].plot(err, E, label="early (+0.5 chip)")
ax[0].plot(err, P, label="prompt")
ax[0].plot(err, Lc, label="late (-0.5 chip)")
ax[0].plot(err, 0.5 * (E - Lc) / np.maximum(E + Lc, 1e-9), "k--", label="DLL discriminator (E-L)/(E+L)/2")
ax[0].set(title="code loop: three correlators half a chip apart", xlabel="code error (chips)")
ax[0].legend(fontsize=7)
ph = np.linspace(-np.pi, np.pi, 361)
ax[1].plot(np.degrees(ph), np.degrees(np.arctan(np.tan(ph))), label="atan(Q/I): blind to the data flip")
ax[1].plot(np.degrees(ph), np.degrees(ph), "--", label="atan2(Q, I): sees a flip as 180 deg")
ax[1].set(title="carrier loop: the Costas discriminator", xlabel="phase error (deg)", ylabel="reported error (deg)")
ax[1].legend(fontsize=7)
save("walk_5_discriminators.png")

# ---- 6. the navigation message: parity-clean subframes and what the anchor is
dec = nav.NavDecoder(7)
for k, ip in enumerate(I):
    dec.feed(ip, k)
fig, ax = plt.subplots(figsize=(9, 2.4))
ax.plot(tt, np.sign(I), lw=0.5)
for i, (sf, tow, bit_i) in enumerate(dec.subframes):
    per = dec.anchors[i][0]
    ax.axvline(per * 1e-3, color="r")
    ax.text(per * 1e-3 + 0.01, 0.6, f"subframe {sf}\nTOW {tow:.0f}\nperiod {per}", fontsize=7, color="r")
ax.set(title="the prompt's sign is the bit stream; a parity-clean subframe ANCHORS a code period to a TOW - the receiver's clock",
       xlabel="time (s)", ylabel="sign(I)", ylim=(-1.4, 1.4))
save("walk_6_nav.png")
