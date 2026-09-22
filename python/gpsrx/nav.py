#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""The navigation message, streaming: 1 ms prompts in, ephemeris and timing anchors out.

IS-GPS-200 §20.3: 50 bit/s, 20 code periods per bit; a 30-bit word carries 24 data bits and
6 parity bits (the Hamming (32,26) code with the two trailing bits of the previous word); ten
words make a 300-bit, 6 s subframe; five subframes make a 30 s frame. Subframe 1 has the clock,
2 and 3 the orbit (the ephemeris), 4 and 5 almanac and, on page 18 of 4, the ionosphere model.

Stages, each a small class with one method, so the block can show each one on a panel:

  BitSync   - which of the 20 code periods a bit starts on: the histogram of sign changes.
  Bits      - 20-period sums -> bits, with a running index that ties every bit to the code
              period it began on (the Nav Decoder's contribution to the pseudorange).
  Framer    - preamble 10001011 every 300 bits, parity on the next words, in either polarity.
  Ephemeris - the field parser, numpy-gps's (relativity.parse_harvest), unchanged in substance.

Timing anchor: every subframe's word 1 carries the time of week of the NEXT subframe's start
(HOW: TOW count x 6 s). The bit index of the subframe's first bit and its code-period index are
therefore a (period, GPS time) pair: from then on, every code period the channel counts is
1 ms of satellite time. That pair is what the PVT solver needs from each channel."""
import numpy as np

PREAMBLE = np.array([1, 0, 0, 0, 1, 0, 1, 1], dtype=np.int8)
# IS-GPS-200 Table 20-XIV: which data bits (1-indexed) each parity bit XORs, and D29*/D30*
PAR = [
    (29, [1, 2, 3, 5, 6, 10, 11, 12, 13, 14, 17, 18, 20, 23]),
    (30, [2, 3, 4, 6, 7, 11, 12, 13, 14, 15, 18, 19, 21, 24]),
    (29, [1, 3, 4, 5, 7, 8, 12, 13, 14, 15, 16, 19, 20, 22]),
    (30, [2, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 20, 21, 23]),
    (30, [1, 3, 5, 6, 7, 9, 10, 14, 15, 16, 17, 18, 21, 22, 24]),
    (29, [3, 5, 6, 8, 9, 10, 11, 13, 15, 19, 22, 23, 24]),
]
PI = 3.1415926535898          # the GPS pi (IS-GPS-200 §20.3.3.4.3)


# ---- bits from prompts -------------------------------------------------------------------
class BitSync:
    """Find the bit boundary: over N prompts, count at which period (mod 20) the prompt sign
    flips. The true boundary collects the flips; the rest is noise."""

    def __init__(self, need=200):
        self.need = need
        self.hist = np.zeros(20, np.int64)
        self.n = 0
        self.last = 0.0
        self.offset = None

    def feed(self, i_prompt, period_index):
        if self.last != 0.0 and np.sign(i_prompt) != np.sign(self.last):
            self.hist[period_index % 20] += 1
        self.last = i_prompt
        self.n += 1
        if self.offset is None and self.n >= self.need and self.hist.max() >= 8:
            top = np.sort(self.hist)[::-1]
            if top[0] >= 3 * max(top[1], 1):                # the winner must stand clear
                self.offset = int(np.argmax(self.hist))
        return self.offset


class Bits:
    """20-period sums on the bit grid -> bits, each tagged with the code-period index it began on."""

    def __init__(self, offset):
        self.offset = offset
        self.acc = 0.0
        self.count = 0
        self.start_period = None
        self.bits = []                # (bit, start_period_index, confidence)

    def feed(self, i_prompt, period_index):
        if (period_index - self.offset) % 20 == 0:
            if self.count == 20:
                self.bits.append((1 if self.acc > 0 else 0, self.start_period, abs(self.acc)))
            self.acc, self.count, self.start_period = 0.0, 0, period_index
        if self.start_period is not None:
            self.acc += i_prompt
            self.count += 1
        return len(self.bits)


# ---- words and subframes -----------------------------------------------------------------
def parity_ok(word, d29s, d30s):
    """word: 30 transmitted bits. Returns (ok, 24 data bits) with the complement convention:
    d_i = D_i XOR D30* for the data bits (IS-GPS-200 §20.3.5.2)."""
    D = np.asarray(word, np.int8)
    d = D[:24] ^ d30s
    for k, (star, mask) in enumerate(PAR):
        acc = d29s if star == 29 else d30s
        for m in mask:
            acc ^= d[m - 1]
        if acc != D[24 + k]:
            return False, d
    return True, d


def ubits(d, a, b):
    v = 0
    for k in range(a - 1, b):
        v = (v << 1) | int(d[k])
    return v


def sbits(d, a, b):
    v = ubits(d, a, b)
    n = b - a + 1
    return v - (1 << n) if v >= (1 << (n - 1)) else v


def cat(hi, hi_n, lo, lo_n, signed):
    v = (hi << lo_n) | lo
    n = hi_n + lo_n
    if signed and v >= (1 << (n - 1)):
        v -= 1 << n
    return v


class Framer:
    """Slides over the bit stream looking for a subframe: preamble in either polarity, then all
    ten words passing parity. Yields (subframe_id, tow_s, words, first_bit_index, polarity)."""

    def __init__(self):
        self.pos = 0                  # next bit index to examine

    def scan(self, bits):
        """bits: the growing list of (bit, start_period, conf). Returns new subframes."""
        out = []
        b = np.array([x[0] for x in bits], np.int8)
        while self.pos + 302 <= len(b):
            i = self.pos
            found = None
            for pol in (0, 1):
                bb = b ^ pol
                if not np.array_equal(bb[i:i + 8], PREAMBLE):
                    continue
                d29s, d30s = (int(bb[i - 2]), int(bb[i - 1])) if i >= 2 else (0, 0)
                words = {}
                ok_all = True
                for w in range(10):
                    word = bb[i + w * 30:i + (w + 1) * 30]
                    ok, d = parity_ok(word, d29s, d30s)
                    if not ok:
                        ok_all = False
                        break
                    words[w] = d
                    d29s, d30s = int(word[28]), int(word[29])
                if ok_all:
                    sf = ubits(words[1], 20, 22)
                    if 1 <= sf <= 5:
                        tow = ubits(words[1], 1, 17) * 6.0        # HOW: TOW count of the NEXT subframe
                        found = (sf, tow, words, i, pol)
                        break
            if found:
                out.append(found)
                self.pos = i + 300
            else:
                self.pos += 1
        return out


# ---- the ephemeris ------------------------------------------------------------------------
def parse_subframe(eph, sfid, words):
    """Fill `eph` from one parity-clean subframe's words (numpy-gps's parse_harvest, per subframe)."""
    W = words

    def put(key, val):
        eph[key] = val

    if sfid == 1:
        if 2 in W:
            put("WN", ubits(W[2], 1, 10))
        if 8 in W:
            put("af2", sbits(W[8], 1, 8) * 2 ** -55)
            put("af1", sbits(W[8], 9, 24) * 2 ** -43)
        if 9 in W:
            put("af0", sbits(W[9], 1, 22) * 2 ** -31)
        if 7 in W:
            put("toc", ubits(W[7], 9, 24) * 16)
    elif sfid == 2:
        if 2 in W:
            put("IODE2", ubits(W[2], 1, 8))
            put("Crs", sbits(W[2], 9, 24) * 2 ** -5)
        if 3 in W:
            put("dn", sbits(W[3], 1, 16) * 2 ** -43 * PI)
        if 3 in W and 4 in W:
            put("M0", cat(ubits(W[3], 17, 24), 8, ubits(W[4], 1, 24), 24, True) * 2 ** -31 * PI)
        if 5 in W:
            put("Cuc", sbits(W[5], 1, 16) * 2 ** -29)
        if 5 in W and 6 in W:
            put("e", cat(ubits(W[5], 17, 24), 8, ubits(W[6], 1, 24), 24, False) * 2 ** -33)
        if 7 in W:
            put("Cus", sbits(W[7], 1, 16) * 2 ** -29)
        if 7 in W and 8 in W:
            put("sqrtA", cat(ubits(W[7], 17, 24), 8, ubits(W[8], 1, 24), 24, False) * 2 ** -19)
        if 9 in W:
            put("toe", ubits(W[9], 1, 16) * 16)
    elif sfid == 3:
        if 2 in W:
            put("Cic", sbits(W[2], 1, 16) * 2 ** -29)
        if 2 in W and 3 in W:
            put("Omega0", cat(ubits(W[2], 17, 24), 8, ubits(W[3], 1, 24), 24, True) * 2 ** -31 * PI)
        if 4 in W:
            put("Cis", sbits(W[4], 1, 16) * 2 ** -29)
        if 4 in W and 5 in W:
            put("i0", cat(ubits(W[4], 17, 24), 8, ubits(W[5], 1, 24), 24, True) * 2 ** -31 * PI)
        if 6 in W:
            put("Crc", sbits(W[6], 1, 16) * 2 ** -5)
        if 6 in W and 7 in W:
            put("omega", cat(ubits(W[6], 17, 24), 8, ubits(W[7], 1, 24), 24, True) * 2 ** -31 * PI)
        if 8 in W:
            put("OmegaDot", sbits(W[8], 1, 24) * 2 ** -43 * PI)
        if 9 in W:
            put("IODE3", ubits(W[9], 1, 8))
            put("IDOT", sbits(W[9], 9, 22) * 2 ** -43 * PI)
    elif sfid == 4:
        if 2 in W and 3 in W and 4 in W and ubits(W[2], 3, 8) == 56:      # page 18: Klobuchar
            put("iono_a", [sbits(W[2], 9, 16) * 2 ** -30, sbits(W[2], 17, 24) * 2 ** -27,
                           sbits(W[3], 1, 8) * 2 ** -24, sbits(W[3], 9, 16) * 2 ** -24])
            put("iono_b", [sbits(W[3], 17, 24) * 2 ** 11, sbits(W[4], 1, 8) * 2 ** 14,
                           sbits(W[4], 9, 16) * 2 ** 16, sbits(W[4], 17, 24) * 2 ** 16])
    return eph


EPH_KEYS = ("WN", "af0", "af1", "af2", "toc", "Crs", "dn", "M0", "Cuc", "e", "Cus", "sqrtA", "toe",
            "Cic", "Omega0", "Cis", "i0", "Crc", "omega", "OmegaDot", "IDOT")


def complete(eph):
    return all(k in eph for k in EPH_KEYS)


# ---- the whole decoder for one channel ----------------------------------------------------
class NavDecoder:
    """Feed it (I_prompt, period_index) per code period; read `eph`, `anchors`, `subframes`."""

    def __init__(self, prn):
        self.prn = prn
        self.sync = BitSync()
        self.bits = None
        self.framer = Framer()
        self.eph = {"prn": prn}
        self.subframes = []           # (sfid, tow, first_bit_index)
        self.anchors = []             # (period_index of the subframe's first bit, tow of THAT subframe start)
        self.n_periods = 0

    def feed(self, i_prompt, period_index):
        self.n_periods += 1
        if self.bits is None:
            off = self.sync.feed(i_prompt, period_index)
            if off is not None:
                self.bits = Bits(off)
            return []
        self.bits.feed(i_prompt, period_index)
        if len(self.bits.bits) < self.framer.pos + 302:
            return []
        new = self.framer.scan(self.bits.bits)
        for sf, tow_next, words, i, pol in new:
            parse_subframe(self.eph, sf, words)
            # the HOW's TOW is the start of the NEXT subframe; THIS subframe began 6 s earlier,
            # at bit i, which began on code period self.bits.bits[i][1]
            tow_this = tow_next - 6.0
            self.subframes.append((sf, tow_this, i))
            self.anchors.append((int(self.bits.bits[i][1]), float(tow_this)))
        return new

    @property
    def complete(self):
        return complete(self.eph)
