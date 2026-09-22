# SPDX-License-Identifier: GPL-3.0-or-later
"""Acquisition: which satellites are in the sky, and where each one's code and carrier are.

Parallel code-phase search (Borre et al. ch. 6; the same search numpy-gps uses): for each
Doppler bin, wipe the carrier off 1 ms blocks, FFT them, multiply by the conjugate FFT of the
sampled C/A code, inverse FFT -> the circular correlation at every code phase in one go. Take
the power, sum it over `n_noncoh` blocks (non-coherent: the data bit may flip between blocks),
and the peak of the Doppler x code-phase map is the satellite. The detection metric is peak
over the second peak (excluding +-1 chip of the first at every Doppler): a satellite that is
there gives one clear peak; noise gives many of the same height.

The coarse Doppler bin (250 Hz) is wider than the tracking loop's pull-in, so `refine_doppler`
runs 100 ms of open-loop prompts at the coarse estimate and reads the residual tone's frequency
off an FFT of the squared prompts (squaring removes the data bits; halve the result).

The engine is pure NumPy; the Acquisition block runs it on a snapshot in its own thread."""
import numpy as np

from .cacode import CODE_RATE, sampled_code
from .track import track_array

ALL_PRNS = list(range(1, 33))


def search(x, fs, prns=ALL_PRNS, dopplers=None, n_noncoh=40):
    """x: complex baseband, at least n_noncoh ms. Returns {prn: {metric, doppler_hz, code_phase
    (samples, fractional: where the code starts within the first ms), peak_over_floor}}."""
    if dopplers is None:
        dopplers = np.arange(-7000.0, 7001.0, 250.0)
    n1 = int(round(fs * 1e-3))
    blocks = x[: n1 * n_noncoh].reshape(n_noncoh, n1)
    t = np.arange(n1) / fs
    excl = int(round(fs / CODE_RATE))                     # +-1 chip in samples
    code_f = {p: np.conj(np.fft.fft(sampled_code(p, fs, n1))) for p in prns}
    maps = {p: np.empty((len(dopplers), n1)) for p in prns}
    for di, fd in enumerate(dopplers):
        BF = np.fft.fft(blocks * np.exp(-2j * np.pi * fd * t)[None, :], axis=1)
        for p in prns:
            corr = np.fft.ifft(BF * code_f[p][None, :], axis=1)
            maps[p][di] = (corr.real ** 2 + corr.imag ** 2).sum(axis=0)
    out = {}
    for p in prns:
        m = maps[p]
        di, ci = np.unravel_index(int(np.argmax(m)), m.shape)
        mask = np.ones(n1, bool)
        mask[[(ci + off) % n1 for off in range(-excl, excl + 1)]] = False
        row = m[di]
        y1, y2, y3 = row[(ci - 1) % n1], row[ci], row[(ci + 1) % n1]
        den = y1 - 2.0 * y2 + y3
        frac = float(np.clip(0.5 * (y1 - y3) / den, -0.5, 0.5)) if den else 0.0     # parabolic sub-sample peak
        out[p] = dict(metric=float(m[di, ci] / m[:, mask].max()), doppler_hz=float(dopplers[di]),
                      code_phase=float(ci + frac), peak_over_floor=float(m[di, ci] / np.median(m)))
    return out


def refine_doppler(x, fs, prn, doppler_hz, code_phase, n_ms=100):
    """Doppler to ~1/(2 n_ms) kHz: open-loop prompts at the coarse estimate, FFT of their square."""
    pr, _ = track_array(x, fs, prn, doppler_hz, code_phase, n_ms, open_loop=True)
    nfft = 4096
    spec = np.abs(np.fft.fftshift(np.fft.fft(pr ** 2, nfft)))
    fax = np.fft.fftshift(np.fft.fftfreq(nfft, 1e-3))
    return float(doppler_hz + fax[int(np.argmax(spec))] / 2.0)


def lo_offset(x, fs, span_hz=50000.0, step_hz=500.0, n_noncoh=10, prns=ALL_PRNS, threshold=2.5):
    """Where is the radio's local oscillator? A cheap crystal (an RTL-SDR: ~30 ppm) puts every
    satellite up to +-47 kHz from where it should be, far outside a +-7 kHz search. One wide,
    coarse pass (500 Hz steps: a 1 ms coherent sum resolves ~1 kHz) over a few blocks finds the
    strongest satellites wherever they are; the median of their Dopplers is the LO's offset to
    within the constellation's own +-5 kHz spread. Returns (offset_hz, n_found)."""
    dops = np.arange(-span_hz, span_hz + 1, step_hz)
    found = [r["doppler_hz"] for r in search(x, fs, prns, dopplers=dops, n_noncoh=n_noncoh).values()
             if r["metric"] > threshold]
    if not found:
        return 0.0, 0
    return float(np.median(found)), len(found)


def sky(x, fs, threshold=2.5, n_noncoh=40, prns=ALL_PRNS, refine=True, centre_hz=0.0, doppler_max=7000.0):
    """The satellites present, refined, best first: [{prn, metric, doppler_hz, code_phase}].
    centre_hz: the LO offset from lo_offset(), so the +-doppler_max window sits on the sky."""
    dops = np.arange(centre_hz - doppler_max, centre_hz + doppler_max + 1, 250.0)
    found = [dict(prn=p, **r) for p, r in search(x, fs, prns, dopplers=dops, n_noncoh=n_noncoh).items()
             if r["metric"] > threshold]
    found.sort(key=lambda r: -r["metric"])
    if refine:
        n_ref = int(0.11 * fs)
        for r in found:
            if len(x) >= n_ref:
                r["doppler_hz"] = refine_doppler(x[:n_ref], fs, r["prn"], r["doppler_hz"], r["code_phase"])
    return found
