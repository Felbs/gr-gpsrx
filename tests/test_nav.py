# SPDX-License-Identifier: GPL-3.0-or-later
"""The navigation decoder: parity, framing, the ephemeris round trip, and the timing anchor -
first from clean bits, then through the tracker from a synthetic satellite."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import nav, navgen, synth  # noqa: E402
from gpsrx.track import Channel  # noqa: E402

FS = 2.048e6
EPH = navgen.EXAMPLE_EPH


def test_parity_round_trip_and_detects_a_flipped_bit():
    d = np.random.default_rng(3).integers(0, 2, 24).astype(np.int8)
    tx = navgen.word_with_parity(d, 1, 0)
    ok, back = nav.parity_ok(tx, 1, 0)
    assert ok and np.array_equal(back, d)
    bad = tx.copy()
    bad[7] ^= 1
    assert not nav.parity_ok(bad, 1, 0)[0]


def test_ephemeris_round_trip_from_clean_bits():
    bits = navgen.frame_bits(EPH, tow_count0=50400, n_subframes=8)
    dec = nav.NavDecoder(7)
    dec.bits = nav.Bits(0)                                  # skip bit sync: feed bits as 20-period sums
    for k, b in enumerate(bits):
        for p in range(20):
            dec.bits.feed(1.0 if b else -1.0, k * 20 + p)
    new = dec.framer.scan(dec.bits.bits)
    assert [s[0] for s in new] == [1, 2, 3, 4, 5, 1, 2]        # the 8th is still accumulating its last bit
    for sf, tow_next, words, i, pol in new:
        nav.parse_subframe(dec.eph, sf, words)
    assert nav.complete(dec.eph)
    # every field comes back within ONE LSB of its broadcast scale factor (IS-GPS-200 Table 20-I/III)
    lsb = dict(WN=1, af0=2 ** -31, af1=2 ** -43, af2=2 ** -55, toc=16, Crs=2 ** -5, dn=2 ** -43 * nav.PI,
               M0=2 ** -31 * nav.PI, Cuc=2 ** -29, e=2 ** -33, Cus=2 ** -29, sqrtA=2 ** -19, toe=16, Cic=2 ** -29,
               Omega0=2 ** -31 * nav.PI, Cis=2 ** -29, i0=2 ** -31 * nav.PI, Crc=2 ** -5, omega=2 ** -31 * nav.PI,
               OmegaDot=2 ** -43 * nav.PI, IDOT=2 ** -43 * nav.PI)
    for k in nav.EPH_KEYS:
        assert abs(dec.eph[k] - EPH[k]) <= lsb[k] * 1.01, (k, dec.eph[k], EPH[k])
    # the first subframe's HOW says TOW count 50401 -> next subframe starts at 302406 s; this one at 302400
    assert new[0][1] == 50401 * 6.0


def test_decoder_from_a_tracked_synthetic_satellite_gives_ephemeris_and_anchor():
    bits = navgen.frame_bits(EPH, tow_count0=50400, n_subframes=12)          # 72 s
    pm = np.where(bits > 0, 1.0, -1.0)
    x, tr = synth.satellite(7, FS, 45.0, doppler_hz=987.0, code_phase_samples=333, cn0_dbhz=46.0, bits=pm)
    ch = Channel(7, FS, 987.0 + 25.0, 333.0)
    dec = nav.NavDecoder(7)
    xd = x.astype(np.complex128)
    pos, k = 0, 0
    while True:
        n = ch.samples_needed()
        if pos + n > len(xd):
            break
        ip, qp = ch.step(xd[pos:pos + n])
        pos += n
        dec.feed(ip, k)
        k += 1
    assert dec.bits is not None, "bit sync never found the boundary"
    assert len(dec.subframes) >= 5, dec.subframes
    assert dec.complete, [k for k in nav.EPH_KEYS if k not in dec.eph]
    assert abs(dec.eph["sqrtA"] - EPH["sqrtA"]) < 1e-3
    assert abs(dec.eph["e"] - EPH["e"]) < 1e-9
    # timing anchor: subframe k starts at TOW 302400 + 6k, and the anchor's period index must
    # advance by exactly 6000 periods between consecutive subframes
    ps = [a[0] for a in dec.anchors]
    tows = [a[1] for a in dec.anchors]
    assert all(np.diff(ps) == 6000), np.diff(ps)
    assert all(abs(np.diff(tows) - 6.0) < 1e-9)
    # ...and the first bit of the first found subframe began on a period that is a multiple of 20
    # from the bit grid the synthetic satellite used (bits are code-period aligned, 20 per bit)
    assert (ps[0] - dec.bits.offset) % 20 == 0
