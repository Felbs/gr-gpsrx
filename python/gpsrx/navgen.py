#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""A navigation-message ENCODER for tests: an ephemeris -> the 50 bit/s bit stream a satellite
would transmit (subframes 1-3 with real parity, TLM and HOW words, subframes 4/5 as filler).
The decoder must get the ephemeris back exactly. Also the inverse of the field parser, so a
mistake in either scale factor shows up as a round-trip failure."""
import numpy as np

from .nav import PAR, PI, PREAMBLE


def _bits(v, n):
    v = int(v) & ((1 << n) - 1)
    return [(v >> (n - 1 - k)) & 1 for k in range(n)]


def _u(x, n, scale):
    return _bits(int(round(x / scale)), n)


def _s(x, n, scale):
    return _bits(int(round(x / scale)), n)          # two's complement via the mask in _bits


def word_with_parity(data24, d29s, d30s):
    """Transmitted 30 bits from 24 data bits (source bits, before the complement convention)."""
    d = np.asarray(data24, np.int8)
    D = list((d ^ d30s).tolist())                   # transmitted data bits = d XOR D30*
    for star, mask in PAR:
        acc = d29s if star == 29 else d30s
        for m in mask:
            acc ^= int(d[m - 1])
        D.append(acc)
    return np.array(D, np.int8)


def subframe_words(sfid, tow_count, eph):
    """The ten 24-bit data words of subframe sfid (words 0..9 = TLM, HOW, 3..10)."""
    e = eph
    W = [None] * 10
    W[0] = PREAMBLE.tolist() + _bits(0, 14) + [0, 0]                       # TLM: preamble, message, reserved
    W[1] = _bits(tow_count, 17) + [0, 0] + _bits(sfid, 3) + [0, 0]          # HOW: TOW, alert, A/S, subframe id, parity fill
    if sfid == 1:
        W[2] = _bits(e["WN"], 10) + _bits(0, 2) + _bits(0, 4) + _bits(0, 6) + _bits(0, 2)
        W[3] = _bits(0, 24)
        W[4] = _bits(0, 24)
        W[5] = _bits(0, 24)
        W[6] = _bits(0, 24)
        W[7] = _bits(0, 8) + _u(e["toc"], 16, 16)                          # word 8: IODC lsbs, toc
        W[8] = _s(e["af2"], 8, 2 ** -55) + _s(e["af1"], 16, 2 ** -43)      # word 9
        W[9] = _s(e["af0"], 22, 2 ** -31) + [0, 0]                         # word 10
    elif sfid == 2:
        M0 = _s(e["M0"], 32, 2 ** -31 * PI)
        ecc = _u(e["e"], 32, 2 ** -33)
        sqrtA = _u(e["sqrtA"], 32, 2 ** -19)
        W[2] = _bits(e["IODE2"], 8) + _s(e["Crs"], 16, 2 ** -5)
        W[3] = _s(e["dn"], 16, 2 ** -43 * PI) + M0[:8]
        W[4] = M0[8:]
        W[5] = _s(e["Cuc"], 16, 2 ** -29) + ecc[:8]
        W[6] = ecc[8:]
        W[7] = _s(e["Cus"], 16, 2 ** -29) + sqrtA[:8]
        W[8] = sqrtA[8:]
        W[9] = _u(e["toe"], 16, 16) + _bits(0, 8)
    elif sfid == 3:
        Om0 = _s(e["Omega0"], 32, 2 ** -31 * PI)
        i0 = _s(e["i0"], 32, 2 ** -31 * PI)
        om = _s(e["omega"], 32, 2 ** -31 * PI)
        W[2] = _s(e["Cic"], 16, 2 ** -29) + Om0[:8]
        W[3] = Om0[8:]
        W[4] = _s(e["Cis"], 16, 2 ** -29) + i0[:8]
        W[5] = i0[8:]
        W[6] = _s(e["Crc"], 16, 2 ** -5) + om[:8]
        W[7] = om[8:]
        W[8] = _s(e["OmegaDot"], 24, 2 ** -43 * PI)
        W[9] = _bits(e["IODE3"], 8) + _s(e["IDOT"], 14, 2 ** -43 * PI) + [0, 0]
    else:
        for w in range(2, 10):
            W[w] = _bits(0x5A5A5A, 24)
    return W


def frame_bits(eph, tow_count0, n_subframes=15):
    """A stream of n_subframes consecutive subframes starting at subframe 1 with HOW TOW count
    tow_count0 (the count of the NEXT subframe = this one's start + 6 s, per IS-GPS-200)."""
    out = []
    d29s = d30s = 0
    for k in range(n_subframes):
        sfid = k % 5 + 1
        words = subframe_words(sfid, tow_count0 + k + 1, eph)
        for w in range(10):
            tx = word_with_parity(words[w], d29s, d30s)
            out.extend(tx.tolist())
            d29s, d30s = int(tx[28]), int(tx[29])
    return np.array(out, np.int8)


EXAMPLE_EPH = dict(prn=7, WN=345, toc=302400, af0=1.234e-4, af1=-2.5e-12, af2=0.0, IODE2=77, IODE3=77,
                   Crs=-23.5, dn=4.5e-9, M0=1.234, Cuc=-1.2e-6, e=0.0123, Cus=8.7e-6, sqrtA=5153.7,
                   toe=302400, Cic=3.0e-8, Omega0=-2.1, Cis=-1.5e-7, i0=0.96, Crc=250.0, omega=0.7,
                   OmegaDot=-8.1e-9, IDOT=2.0e-10)
