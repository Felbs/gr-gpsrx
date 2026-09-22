#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""The GPS L1 C/A code (IS-GPS-200 §3.3.2.3) and the constants every stage shares.

The generator is numpy-gps's (MIT, Felbs), kept bit for bit: two 10-stage LFSRs, G1 = x^10+x^3+1,
G2 = x^10+x^9+x^8+x^6+x^3+x^2+1, the PRN's two G2 taps XORed with G1's output. The published
first-ten-chip octal values are the self-check."""
import numpy as np

CODE_RATE = 1.023e6            # chips / s
CODE_LEN = 1023                # chips per code period = 1 ms
L1_HZ = 1575.42e6
C_LIGHT = 299792458.0

# IS-GPS-200 Table 3-Ia: G2 phase-select taps for PRN 1..32
G2_TAPS = {1: (2, 6), 2: (3, 7), 3: (4, 8), 4: (5, 9), 5: (1, 9), 6: (2, 10), 7: (1, 8), 8: (2, 9),
           9: (3, 10), 10: (2, 3), 11: (3, 4), 12: (5, 6), 13: (6, 7), 14: (7, 8), 15: (8, 9), 16: (9, 10),
           17: (1, 4), 18: (2, 5), 19: (3, 6), 20: (4, 7), 21: (5, 8), 22: (6, 9), 23: (1, 3), 24: (4, 6),
           25: (5, 7), 26: (6, 8), 27: (7, 9), 28: (8, 10), 29: (1, 6), 30: (2, 7), 31: (3, 8), 32: (4, 9)}
# IS-GPS-200 Table 3-Ia: first 10 chips, octal
PUBLISHED_OCTAL = {1: "1440", 2: "1620", 3: "1710", 4: "1744", 5: "1133", 6: "1455", 7: "1131", 8: "1454",
                   9: "1626", 10: "1504", 11: "1642", 12: "1750", 13: "1764", 14: "1772", 15: "1775",
                   16: "1776", 17: "1156", 18: "1467", 19: "1633", 20: "1715", 21: "1746", 22: "1763",
                   23: "1063", 24: "1706", 25: "1743", 26: "1761", 27: "1770", 28: "1774", 29: "1127",
                   30: "1453", 31: "1625", 32: "1712"}

_CA = {}
_SAMPLED = {}


def ca_code(prn):
    """The 1023-chip C/A code as +1/-1 float64, read-only, cached."""
    c = _CA.get(prn)
    if c is not None:
        return c
    s1, s2 = G2_TAPS[prn]
    g1 = np.ones(10, dtype=int)
    g2 = np.ones(10, dtype=int)
    out = np.empty(CODE_LEN, dtype=int)
    for i in range(CODE_LEN):
        out[i] = g1[9] ^ g2[s1 - 1] ^ g2[s2 - 1]
        fb1 = g1[2] ^ g1[9]
        fb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g1 = np.concatenate(([fb1], g1[:9]))
        g2 = np.concatenate(([fb2], g2[:9]))
    c = 1.0 - 2.0 * out
    c.setflags(write=False)
    _CA[prn] = c
    return c


def selfcheck():
    """Every PRN's first ten chips against the published octal. Raises on a mismatch."""
    for prn, want in PUBLISHED_OCTAL.items():
        chips = (1.0 - ca_code(prn)[:10]) / 2.0
        got = format(int("".join(str(int(b)) for b in chips), 2), "04o")
        if got != want:
            raise RuntimeError(f"C/A generator: PRN {prn} first chips {got}, published {want}")
    return True


def sampled_code(prn, fs, n_samp):
    """The code sampled at fs over n_samp samples at zero phase, nominal chip rate (no Doppler on
    the code). Read-only, cached. numpy-gps's acquisition/prompt reference uses exactly this."""
    key = (prn, float(fs), int(n_samp))
    c = _SAMPLED.get(key)
    if c is None:
        idx = (np.arange(n_samp) * CODE_RATE / fs).astype(np.int64) % CODE_LEN
        c = ca_code(prn)[idx]
        c.setflags(write=False)
        _SAMPLED[key] = c
    return c


def code_at(prn, phase_chips):
    """The code value at fractional chip phases (any shape): a nearest-chip lookup, the way a
    hardware correlator's code generator works. phase in chips, any real value."""
    idx = np.floor(phase_chips).astype(np.int64) % CODE_LEN
    return ca_code(prn)[idx]
