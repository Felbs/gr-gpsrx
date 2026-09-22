#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""A synthetic GPS satellite: known PRN, Doppler, code phase, data bits, C/N0.

For QA and gates that must run with no capture and no radio. The signal is exactly what the
receiver models - a C/A code at 1.023 MHz x (1 + Doppler/L1), BPSK data at 50 bit/s, a carrier
at the Doppler offset, white Gaussian noise - so it tests the receiver against its own
assumptions. Real air (numpy-gps's captures) tests the assumptions."""
import numpy as np

from .cacode import CODE_LEN, CODE_RATE, L1_HZ, code_at


def satellite(prn, fs, secs, doppler_hz=1234.5, code_phase_samples=700, cn0_dbhz=45.0, bits=None,
              carrier_phase=0.3, seed=1, code_doppler=True, noise=True):
    """Baseband complex64 samples. Noise power 1 per sample (complex); signal power set from C/N0
    at this sample rate. `bits`: array of +-1 at 50 bit/s (random if None). Returns (x, truth)."""
    rng = np.random.default_rng(seed)
    n = int(round(fs * secs))
    t = np.arange(n) / fs
    rate = CODE_RATE * (1 + doppler_hz / L1_HZ) if code_doppler else CODE_RATE
    # code phase: the code STARTS at sample code_phase_samples
    ph = (t - code_phase_samples / fs) * rate
    code = code_at(prn, ph)
    n_bits = int(np.ceil(secs * 50)) + 1
    if bits is None:
        bits = rng.choice([-1.0, 1.0], n_bits)
    bits = np.asarray(bits, np.float64)
    bit_of_sample = np.floor(np.maximum(ph, 0) / (20 * CODE_LEN)).astype(np.int64)      # 20 code periods per bit
    data = bits[np.minimum(bit_of_sample, len(bits) - 1)]
    # amplitude from C/N0: complex noise of unit power over bandwidth fs -> N0 = 1/fs; C = A^2
    amp = np.sqrt(10 ** (cn0_dbhz / 10) / fs)
    sig = amp * code * data * np.exp(1j * (carrier_phase + 2 * np.pi * doppler_hz * t))
    if noise:
        sig = sig + (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    x = sig.astype(np.complex64) if noise else sig
    truth = {"prn": prn, "fs": fs, "doppler_hz": doppler_hz, "code_phase_samples": code_phase_samples,
             "code_rate": rate, "cn0_dbhz": cn0_dbhz, "bits": bits, "carrier_phase": carrier_phase, "amp": amp}
    return x, truth


def sky(fs, secs, birds, seed=0):
    """Several satellites over ONE noise floor: `birds` = list of dicts with prn, doppler_hz,
    code_phase_samples, cn0_dbhz (and optionally bits). Returns (x, truths). Summing single-
    satellite signals adds their noise floors too - this is the right way."""
    rng = np.random.default_rng(seed)
    n = int(round(fs * secs))
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    truths = []
    for b in birds:
        xi, tr = satellite(b["prn"], fs, secs, doppler_hz=b.get("doppler_hz", 0.0),
                           code_phase_samples=b.get("code_phase_samples", 0), cn0_dbhz=b.get("cn0_dbhz", 45.0),
                           bits=b.get("bits"), seed=int(rng.integers(1 << 30)), noise=False)
        x += xi
        truths.append(tr)
    return x.astype(np.complex64), truths
