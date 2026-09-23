# SPDX-License-Identifier: GPL-3.0-or-later
"""Galileo I/NAV, streamed: E1-B symbols in (one per 4 ms code period), ephemeris and timing
anchors out. The FEC machinery, page framing, CRC and word parser are numpy-gps's (MIT; its
1161/1161 CRC-clean pages on a 13-minute capture are this file's oracle), rearranged so each
250-symbol page part is decoded as soon as it completes.

The message (OS SIS ICD 4.3): every second a page PART of 250 symbols - a 10-symbol sync
pattern, then 240 symbols that are 120 bits of a rate-1/2, K=7 convolutional code (the second
branch inverted), block-interleaved 30 x 8. An EVEN part and the following ODD part make one
page of 2 x 114 bits, protected by CRC-24Q; a page carries one WORD, and word types 1-4 are the
ephemeris (IODnav ties the four together), 5 the BGD/health/GST week and TOW, 0 and 6 more TOW,
10 the GST-GPS time offset. The broadcast TOW marks the first chip of the first symbol of its
EVEN part: that symbol's code period is the timing anchor, exactly as a GPS subframe's first
bit is (nav.py)."""
import numpy as np

SYNC = "0101100000"
G1_OCT, G2_OCT = 0o171, 0o133
CRC24_POLY = 0x1864CFB
SPP = 250                                  # symbols per page part (1 s)
PRE = np.array([1.0 if c == "1" else -1.0 for c in SYNC])
MU_GAL = 3.986004418e14
PI = 3.1415926535898


def _trellis():
    sp = np.arange(64)
    p0 = 2 * (sp & 31)
    p1 = p0 + 1
    b = (sp >> 5) & 1
    v0 = (b << 6) | p0
    v1 = (b << 6) | p1
    pop = np.array([bin(v).count("1") & 1 for v in range(128)])
    s1 = 2.0 * pop[np.arange(128) & G1_OCT] - 1.0
    s2 = 2.0 * pop[np.arange(128) & G2_OCT] - 1.0
    return p0, p1, v0, v1, s1, s2


TP0, TP1, TV0, TV1, TS1, TS2 = _trellis()


def viterbi120(soft240):
    """Soft-decision Viterbi over one page part (120 bits incl. 6 tail bits). Returns (bits, metric)."""
    r = soft240.reshape(120, 2)
    pm = np.full(64, -1e12)
    pm[0] = 0.0
    dec = np.zeros((120, 64), dtype=bool)
    for k in range(120):
        m0 = pm[TP0] + TS1[TV0] * r[k, 0] + TS2[TV0] * r[k, 1]
        m1 = pm[TP1] + TS1[TV1] * r[k, 0] + TS2[TV1] * r[k, 1]
        take1 = m1 > m0
        pm = np.where(take1, m1, m0)
        dec[k] = take1
    s = 0
    bits = np.zeros(120, dtype=np.int8)
    for k in range(119, -1, -1):
        bits[k] = (s >> 5) & 1
        s = TP1[s] if dec[k, s] else TP0[s]
    return bits, float(pm[0])


def conv_encode(bits):
    """Transmit side (tests): 2N symbols with the second branch inverted."""
    out = np.zeros(2 * len(bits), dtype=np.int8)
    s = 0
    pop = [bin(v).count("1") & 1 for v in range(128)]
    for k, b in enumerate(bits):
        v = (int(b) << 6) | s
        out[2 * k] = pop[v & G1_OCT]
        out[2 * k + 1] = pop[v & G2_OCT] ^ 1
        s = (int(b) << 5) | (s >> 1)
    return out


def deinterleave(sym240):
    return sym240.reshape(8, 30).T.ravel()


def interleave(sym240):
    return sym240.reshape(30, 8).T.ravel()


def crc24q(bits):
    r = 0
    for b in bits:
        r = ((r << 1) ^ (CRC24_POLY if ((r >> 23) ^ int(b)) & 1 else 0)) & 0xFFFFFF
    return r


def decode_part(sym250):
    """250 symbols from the sync pattern, polarity already fixed -> (bits120, efficiency, tail ok)."""
    d = deinterleave(np.asarray(sym250[10:250], float).copy())
    d[1::2] *= -1.0
    bits, pm = viterbi120(d)
    eff = pm / (np.abs(d).sum() + 1e-12)
    return bits, eff, bool((bits[114:] == 0).all())


def _u(b, i, n):
    v = 0
    for k in range(i, i + n):
        v = (v << 1) | int(b[k])
    return v


def _s(b, i, n):
    v = _u(b, i, n)
    return v - (1 << n) if v >= (1 << (n - 1)) else v


def parse_word(page):
    """228-bit page -> (word type, fields); numpy-gps's parser, ICD tables 39-45."""
    w = np.concatenate([page[2:114], page[116:132]])
    wt = _u(w, 0, 6)
    f = {}
    if wt == 0:
        f["time_flag"] = _u(w, 6, 2)
        if f["time_flag"] == 2:
            f["WN"] = _u(w, 96, 12)
            f["TOW"] = _u(w, 108, 20)
    elif wt == 1:
        f["IODnav"] = _u(w, 6, 10)
        f["toe"] = _u(w, 16, 14) * 60
        f["M0"] = _s(w, 30, 32) * 2 ** -31 * PI
        f["e"] = _u(w, 62, 32) * 2 ** -33
        f["sqrtA"] = _u(w, 94, 32) * 2 ** -19
    elif wt == 2:
        f["IODnav"] = _u(w, 6, 10)
        f["Omega0"] = _s(w, 16, 32) * 2 ** -31 * PI
        f["i0"] = _s(w, 48, 32) * 2 ** -31 * PI
        f["omega"] = _s(w, 80, 32) * 2 ** -31 * PI
        f["IDOT"] = _s(w, 112, 14) * 2 ** -43 * PI
    elif wt == 3:
        f["IODnav"] = _u(w, 6, 10)
        f["OmegaDot"] = _s(w, 16, 24) * 2 ** -43 * PI
        f["dn"] = _s(w, 40, 16) * 2 ** -43 * PI
        f["Cuc"] = _s(w, 56, 16) * 2 ** -29
        f["Cus"] = _s(w, 72, 16) * 2 ** -29
        f["Crc"] = _s(w, 88, 16) * 2 ** -5
        f["Crs"] = _s(w, 104, 16) * 2 ** -5
        f["SISA"] = _u(w, 120, 8)
    elif wt == 4:
        f["IODnav"] = _u(w, 6, 10)
        f["SVID"] = _u(w, 16, 6)
        f["Cic"] = _s(w, 22, 16) * 2 ** -29
        f["Cis"] = _s(w, 38, 16) * 2 ** -29
        f["toc"] = _u(w, 54, 14) * 60
        f["af0"] = _s(w, 68, 31) * 2 ** -34
        f["af1"] = _s(w, 99, 21) * 2 ** -46
        f["af2"] = _s(w, 120, 6) * 2 ** -59
    elif wt == 5:
        f["ai0"] = _u(w, 6, 11) * 2 ** -2
        f["ai1"] = _s(w, 17, 11) * 2 ** -8
        f["ai2"] = _s(w, 28, 14) * 2 ** -15
        f["BGD_E1E5a"] = _s(w, 47, 10) * 2 ** -32
        f["BGD_E1E5b"] = _s(w, 57, 10) * 2 ** -32
        f["E5bHS"] = _u(w, 67, 2)
        f["E1BHS"] = _u(w, 69, 2)
        f["E5bDVS"] = _u(w, 71, 1)
        f["E1BDVS"] = _u(w, 72, 1)
        f["WN"] = _u(w, 73, 12)
        f["TOW"] = _u(w, 85, 20)
    elif wt == 6:
        f["A0"] = _s(w, 6, 32) * 2 ** -30
        f["A1"] = _s(w, 38, 24) * 2 ** -50
        f["dtLS"] = _s(w, 62, 8)
        f["TOW"] = _u(w, 105, 20)
    elif wt == 10:
        f["ggto_valid"] = not (_u(w, 86, 16) == 0xFFFF and _u(w, 102, 12) == 0xFFF
                               and _u(w, 114, 8) == 0xFF and _u(w, 122, 6) == 0x3F)
        f["A0G"] = _s(w, 86, 16) * 2 ** -35
        f["A1G"] = _s(w, 102, 12) * 2 ** -51
        f["t0G"] = _u(w, 114, 8) * 3600
        f["WN0G"] = _u(w, 122, 6)
    return wt, f


EPH_WORDS = {1: ("toe", "M0", "e", "sqrtA"), 2: ("Omega0", "i0", "omega", "IDOT"),
             3: ("OmegaDot", "dn", "Cuc", "Cus", "Crc", "Crs", "SISA"), 4: ("toc", "af0", "af1", "af2", "Cic", "Cis")}


class INavDecoder:
    """Feed (I_prompt, period_index) per 4 ms E1-B period; read `eph` (complete when words 1-4 of
    one IODnav are in, with af0 folded by BGD(E1,E5b) for a single-frequency E1 user - ICD 5.1.5),
    `anchors` [(period of the even part's first symbol, GST TOW)], `ggto`, `wn`, `pages`.
    The preamble grid is found from the first 600 symbols (both polarities); afterwards every
    250-symbol part is Viterbi-decoded as it completes and paired even/odd through the CRC."""

    def __init__(self, prn, need=600):
        self.prn = prn
        self.need = need
        self.sym = []                 # (period_index, I); self.sym[0] is symbol number self.base
        self.base = 0
        self.grid = None              # (phase, polarity)
        self.words = {}               # IODnav -> {wt: fields}
        self.extras = {}
        self.eph = {"prn": prn, "sys": "GAL"}
        self.anchors = []
        self.pages = []
        self.ggto = None
        self.wn = None
        self.n_parts = self.n_crc = self.n_crc_fail = 0
        self._pending_even = None
        self._next_part = None

    def feed(self, i_prompt, period_index):
        self.sym.append((period_index, float(i_prompt)))
        out = []
        if self.grid is None:
            if len(self.sym) >= self.need:
                self._find_grid()
                if self.grid is not None:
                    self._next_part = self.base + self.grid[0]
            return out
        # decode every complete part from _next_part on (absolute symbol indices; self.sym[0] is
        # symbol self.base)
        while self._next_part is not None and len(self.sym) + self.base >= self._next_part + SPP:
            i = self._next_part - self.base
            seg = np.array([v for _, v in self.sym[i:i + SPP]]) * self.grid[1]
            if np.dot(np.sign(seg[:10]), PRE) < 6:           # the sync pattern is not here: grid slipped
                self.grid = None
                self._next_part = None
                self._pending_even = None
                return out
            bits, eff, tail_ok = decode_part(seg)
            self.n_parts += 1
            part = dict(i=self._next_part, period=self.sym[i][0], bits=bits, eo=int(bits[0]))
            if part["eo"] == 0:
                self._pending_even = part
            elif self._pending_even is not None and part["i"] - self._pending_even["i"] == SPP:
                page = np.concatenate([self._pending_even["bits"][:114], bits[:114]])
                ok = crc24q(page[:196]) == _u(page, 196, 24)
                if ok:
                    self.n_crc += 1
                    wt, f = parse_word(page)
                    self._take(wt, f, self._pending_even["period"])
                    self.pages.append((wt, self._pending_even["period"]))
                    out.append((wt, f))
                else:
                    self.n_crc_fail += 1
                self._pending_even = None
            self._next_part += SPP
        # memory: keep only what the next part needs
        if self._next_part is not None:
            keep_from = self._next_part - self.base
            if keep_from > SPP:
                self.sym = self.sym[keep_from:]
                self.base += keep_from
        return out

    def _find_grid(self):
        s = np.sign(np.array([v for _, v in self.sym]))
        c = np.correlate(s, PRE, mode="valid")
        best = None
        for pol in (1.0, -1.0):
            for ph in range(SPP):
                hits = c[ph::SPP] * pol
                # strength, not a count: two chance 8-of-10 matches tied two perfect 10s once
                score = float(np.sum(np.where(hits >= 8, hits, 0.0)))
                if best is None or score > best[0]:
                    best = (score, ph, pol, int(np.sum(hits >= 8)))
        score, ph, pol, n = best
        if n >= 2 and score >= 18:
            self.grid = (ph, pol)

    def _take(self, wt, f, period):
        if wt in EPH_WORDS:
            iod = f["IODnav"]
            self.words.setdefault(iod, {})[wt] = f
            ws = self.words[iod]
            if all(w in ws for w in (1, 2, 3, 4)):
                e = {"prn": self.prn, "sys": "GAL", "mu": MU_GAL, "IODnav": iod}
                for w in (1, 2, 3, 4):
                    for k in EPH_WORDS[w]:
                        e[k] = float(ws[w][k])
                e.update({k: v for k, v in self.extras.items()})
                if "BGD_E1E5b" in e:
                    e["af0"] -= e["BGD_E1E5b"]                 # single-frequency E1 user clock
                    e["bgd_folded"] = True
                self.eph = e
        elif wt == 5:
            for k in ("BGD_E1E5a", "BGD_E1E5b", "E1BHS", "E1BDVS", "E5bHS", "E5bDVS", "ai0", "ai1", "ai2"):
                self.extras[k] = f[k]
            if "IODnav" in self.eph and "BGD_E1E5b" in f and not self.eph.get("bgd_folded"):
                self.eph["af0"] -= f["BGD_E1E5b"]
                self.eph["bgd_folded"] = True
            self.eph["BGD_E1E5b"] = f["BGD_E1E5b"]
        elif wt == 10 and f.get("ggto_valid"):
            self.ggto = {k: f[k] for k in ("A0G", "A1G", "t0G", "WN0G")}
        if "TOW" in f and f.get("time_flag", 2) == 2:
            self.anchors.append((int(period), float(f["TOW"])))
            if "WN" in f:
                self.wn = int(f["WN"])

    @property
    def complete(self):
        return "IODnav" in self.eph and bool(self.anchors)


def ggto_eval(ggto, wn_gst, tow_gst):
    """ICD 5.1.8: t_Galileo - t_GPS (s) at GST (wn, tow)."""
    dw = (wn_gst - ggto["WN0G"]) % 64
    if dw > 31:
        dw -= 64
    return ggto["A0G"] + ggto["A1G"] * (tow_gst - ggto["t0G"] + 604800.0 * dw)
