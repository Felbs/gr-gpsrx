# SPDX-License-Identifier: GPL-3.0-or-later
"""Galileo E1: the code table, the BOC replica's autocorrelation, the I/NAV decoder on 12 s of
real E1-B symbols (PRN 29, tracked from a wideband attic capture; numpy-gps decoded the same
satellite 387/387 pages CRC-clean), and a synthetic I/NAV page round trip."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import gale1, nav_gal  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def test_e1_codes_and_boc_autocorrelation():
    c = gale1.e1_code(29)
    assert c.shape == (4092,) and set(np.unique(c)) == {-1.0, 1.0}
    assert not np.array_equal(c, gale1.e1_code(14))
    # the BOC(1,1) replica at 2 samples per chip: its autocorrelation is 1 at zero lag and
    # NEGATIVE half a chip away (the subcarrier's signature), which is why E/L sit a quarter chip out
    r = gale1.sampled_code(29, 2.046e6, 8184)
    ac = np.array([np.dot(r, np.roll(r, k)) for k in range(-3, 4)]) / len(r)
    assert ac[3] == 1.0 and ac[2] < -0.3 and ac[4] < -0.3, ac


def test_inav_decoder_on_real_symbols():
    syms = np.load(os.path.join(HERE, "data", "e1b_prn29_symbols.npy"))
    dec = nav_gal.INavDecoder(29)
    words = []
    for k, v in enumerate(syms):
        words += dec.feed(float(v), k)
    assert dec.grid is not None and dec.n_crc >= 4 and dec.n_crc_fail == 0, (dec.grid, dec.n_crc, dec.n_crc_fail)
    assert dec.wn == 1405
    tows = [t for _, t in dec.anchors]
    per = [p for p, _ in dec.anchors]
    assert all(abs((b - a) - (q - p) * 4e-3) < 1e-9 for (p, a), (q, b) in zip(dec.anchors, dec.anchors[1:])), (per, tows)
    assert {w for w, _ in words} >= {1, 3, 5}


def test_inav_page_round_trip():
    """A word-1 page built with the ICD's encoder chain (CRC-24Q, rate-1/2 K=7 with the second branch
    inverted, 30x8 interleave, sync pattern) decodes back to its fields."""
    rng = np.random.default_rng(5)
    even = np.zeros(114, np.int8)
    odd = np.zeros(114, np.int8)
    even[0], odd[0] = 0, 1                                 # even/odd
    payload = rng.integers(0, 2, 128).astype(np.int8)
    payload[:6] = [0, 0, 0, 0, 0, 1]                        # word type 1
    even[2:114] = payload[:112]
    odd[2:18] = payload[112:128]
    page = np.concatenate([even, odd])
    crc = nav_gal.crc24q(page[:196])
    page[196:220] = [(crc >> (23 - i)) & 1 for i in range(24)]

    def part(bits114):
        bits = np.concatenate([bits114, np.zeros(6, np.int8)])          # tail
        sym = nav_gal.interleave(nav_gal.conv_encode(bits).astype(float))
        return np.concatenate([nav_gal.PRE, 2.0 * sym - 1.0])                # bit 1 -> +1, as the sync pattern
    # on air the parts are back to back, one per second: a random prefix, then the page four times
    stream = np.concatenate([rng.choice([-1.0, 1.0], 137)] + [part(page[:114]), part(page[114:])] * 4)
    dec = nav_gal.INavDecoder(1, need=600)
    got = []
    for k, v in enumerate(stream):
        got += dec.feed(v, k)
    assert dec.n_crc >= 1 and dec.n_crc_fail == 0, (dec.n_crc, dec.n_crc_fail)
    wt, f = got[0]
    assert wt == 1 and "sqrtA" in f


def test_pilot_secondary_sync_and_pure_pll():
    """A synthetic E1 satellite (E1-B x symbols + E1-C x CS25) through the pilot-aided channel:
    the secondary-code sync must name the TRUE offset (the channel's first prompt is the partial
    period before the code start, so its epoch k is the satellite's period k-1: offset 1), and
    after it the four-quadrant PLL on the wiped pilot holds phase to a few degrees. An off-by-one
    in the sync's epoch bookkeeping once put the wipe one chip late on 12 of 25 positions: the
    Costas loop was blind to it (it only lost coherent gain), the pure PLL lost lock in a second."""
    from gpsrx import synth, track
    fs = 4.096e6
    x, tr = synth.satellite(11, fs, 2.6, doppler_hz=-1810.0, code_phase_samples=900, cn0_dbhz=45.0, system="GAL", seed=3)
    ch = track.Channel(11, fs, -1790.0, 900, pll_bw=12.0, dll_bw=1.0, pll_bw_narrow=15.0, dll_bw_narrow=0.5,
                       coherent_ms=20, pll_order=3, signal=gale1.SIGNAL_E1)
    pos, errs = 0, []
    while pos + ch.samples_needed() <= len(x):
        n = ch.samples_needed()
        ch.step(x[pos:pos + n])
        pos += n
        if ch.bit_offset is not None and (ch.s.epochs - ch.bit_offset) % ch.coh == 0:
            errs.append(ch.s.pll_e_prev)
    assert ch.bit_offset == 1, ch.bit_offset
    assert len(errs) > 30
    assert np.max(np.abs(errs[5:])) < 0.06, np.max(np.abs(errs[5:]))       # cycles: < 22 degrees, every window
    assert ch.s.lock > 0.95 and abs(ch.s.carrier_hz - tr["doppler_hz"]) < 2.0, (ch.s.lock, ch.s.carrier_hz)
