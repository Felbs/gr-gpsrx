# SPDX-License-Identifier: GPL-3.0-or-later
"""Galileo E1 open service: the codes, the BOC(1,1) replica, and what differs from GPS L1 C/A.

Galileo shares the L1 frequency (E1, 1575.42 MHz) and the 1.023 Mchip/s rate, and differs in
every other way that matters to a receiver:

  * the primary codes are 4092 chips - 4 ms - and are memory codes (tables in the OS SIS ICD,
    Annex C), not shift-register sequences. The table here is the one gnss-sdr carries
    (Galileo_E1.h, GPL-3, checked against the ICD excerpt by numpy-gps), vendored as
    data/galileo_e1_codes.npz: E1-B (data) and E1-C (pilot) for PRN 1..36;
  * every chip carries a BOC(1,1) subcarrier: +1 for its first half, -1 for its second, so the
    correlation peak is narrower than BPSK's and has side lobes half a chip out. Correlators
    are spaced a quarter chip, not a half;
  * the E1-B data symbol is 4 ms long - one per code period - at 250 symbols/s, so there is no
    bit sync and no integration beyond 4 ms on the data channel; E1-C carries a 25-chip
    secondary code instead of data (100 ms), for pilot tracking (not used here yet);
  * the navigation message (I/NAV) is convolutionally coded and interleaved: nav_gal.py.

The channel engine in track.py is the same for both systems; only the code function, the
period and the correlator spacing change (Channel(..., signal=gale1.SIGNAL_E1B))."""
import os

import numpy as np

CODE_RATE = 1.023e6
CODE_LEN = 4092
L1_HZ = 1575.42e6
_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "galileo_e1_codes.npz")
_CODES = {}


def e1_code(prn, component="b"):
    """The 4092-chip primary code as +1/-1, read-only, cached. component 'b' (data) or 'c' (pilot)."""
    key = (component, prn)
    c = _CODES.get(key)
    if c is None:
        d = np.load(_DATA)
        c = np.asarray(d[f"{component}{prn}"], np.float64)
        c.setflags(write=False)
        _CODES[key] = c
    return c


def code_at(prn, phase_chips, component="b"):
    """E1 code value with the sine-phased BOC(1,1) subcarrier at fractional chip phases: the chip
    at floor(phase), times +1 in its first half and -1 in its second (a nearest-half-chip lookup)."""
    ph = np.asarray(phase_chips, np.float64)
    idx = np.floor(ph).astype(np.int64) % CODE_LEN
    half = np.floor(2.0 * (ph - np.floor(ph))).astype(np.int64)          # 0 first half, 1 second
    return e1_code(prn, component)[idx] * (1.0 - 2.0 * half)


def sampled_code(prn, fs, n_samp, component="b"):
    """One code period (or n_samp samples of it) at fs, zero phase, nominal rate, with the BOC
    subcarrier - the acquisition replica. Band-limiting is left to the capture."""
    ph = np.arange(n_samp) * CODE_RATE / fs
    return code_at(prn, ph, component)


# the description a Channel needs to track this signal instead of GPS L1 C/A: E1-B alone (data
# channel: no integration beyond a period), or PILOT-AIDED - the loops on E1-C, whose 25-chip
# secondary code (ICD 3.8.2) is KNOWN: once the channel finds where it is in the sequence it
# wipes the chips off and integrates up to 100 ms, while a fourth correlator on E1-B delivers
# the data symbols for the I/NAV decoder.
CS25 = "0011100000001010110110010"
SECONDARY = np.array([1.0 - 2.0 * int(c) for c in CS25])


def code_at_c(prn, phase_chips):
    return code_at(prn, phase_chips, "c")


SIGNAL_E1B = dict(name="E1B", code_len=CODE_LEN, code_rate=CODE_RATE, code_at=code_at, spacing=0.25,
                  bit_periods=1, coherent_max=1, carrier_hz=L1_HZ)
SIGNAL_E1 = dict(name="E1", code_len=CODE_LEN, code_rate=CODE_RATE, code_at=code_at_c, data_code_at=code_at,
                 spacing=0.25, bit_periods=25, coherent_max=25, secondary=SECONDARY, carrier_hz=L1_HZ)
