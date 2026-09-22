#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 gr-gpsrx authors.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""One tracking channel for one satellite: the textbook receiver, one millisecond at a time.

No GNU Radio here - this is the engine the Channel block wraps, so it can be tested and timed
under plain pytest. Per code period (nominally 1 ms, 2048 samples at 2.048 MS/s):

  1. wipe the carrier off with the NCO (carrier phase + Doppler, both tracked);
  2. correlate with three copies of the C/A code: EARLY, PROMPT, LATE, half a chip apart;
  3. PLL (Costas discriminator on PROMPT: data bits flip the sign, atan(Q/I) does not care) ->
     carrier Doppler;
  4. DLL (normalised early-minus-late envelope) -> code rate, aided by the carrier
     (code Doppler = carrier Doppler / 1540);
  5. advance the code phase by the code rate x the samples consumed; the next code period
     begins at the sample where the code phase crosses 1023 chips.

The integer count of code periods since the channel started, plus the fractional code phase,
IS the pseudorange (modulo the navigation message's time of week, which the Nav Decoder
supplies). That count and phase are published per period as the observable.

Open-loop mode (loop bandwidths 0, integer sample phase) reproduces numpy-gps's prompts_ms()
exactly: that is gate 0."""
from dataclasses import dataclass, field

import numpy as np

from .cacode import CODE_LEN, CODE_RATE, L1_HZ, code_at

CARRIER_TO_CODE = CODE_RATE / L1_HZ     # 1/1540: code Doppler per Hz of carrier Doppler


def loop_gains(bw_hz, zeta=0.707, k=1.0):
    """Second-order loop filter constants (tau1, tau2) for a noise bandwidth, the SoftGNSS /
    Kaplan & Hegarty form (Borre et al. 2007, ch. 7). Per update with discriminator error e
    over interval T the FREQUENCY correction is
        f += (tau2 * (e - e_prev) + T * e) / tau1
    and the NCO's own phase accumulation closes the loop. k is the discriminator gain."""
    wn = bw_hz * 8 * zeta / (4 * zeta * zeta + 1)
    return k / (wn * wn), 2 * zeta / wn


@dataclass
class ChannelState:
    prn: int
    fs: float
    carrier_hz: float               # current carrier Doppler estimate
    code_phase: float               # chips, [0, 1023): phase at the START of the next period
    epochs: int = 0                 # code periods completed since start (the integer pseudorange)
    carrier_phase: float = 0.0      # radians, at the start of the next period
    code_rate: float = CODE_RATE    # chips/s, including code Doppler
    samples_in: int = 0             # samples consumed since start
    pll_e_prev: float = 0.0
    dll_e_prev: float = 0.0
    cn0_db: float = 0.0             # carrier-to-noise density, dB-Hz (moment estimator over 20 prompts)
    carrier_cycles: float = 0.0     # accumulated NCO carrier phase, cycles: the carrier-phase observable
    lock: float = 0.0               # 0..1 PLL lock indicator, smoothed
    prompts: list = field(default_factory=list)


class Channel:
    """Track one PRN through a stream, one code period per step."""
    RAMP_TOL = 5.0          # Hz: rebuild the NCO ramp when the Doppler estimate moves this much (5 Hz over 1 ms = 1.8 deg, under the loop noise)

    def __init__(self, prn, fs, doppler_hz, code_phase_samples, pll_bw=18.0, dll_bw=2.0, spacing=0.5,
                 open_loop=False, pll_bw_narrow=15.0, dll_bw_narrow=0.5, coherent_ms=20):
        # Two stages, as gnss-sdr does it: wide loops and 1 ms integration to pull in; then, once
        # the data-bit edges are known, NARROW loops and coherent integration over a whole bit
        # (coherent_ms, a divisor of 20). The correlators sum across the bit - 13 dB more
        # coherent gain - and the discriminators run once per bit on a 20x quieter number.
        # Measured before this (head-to-head, docs/GNSS_SDR_COMPARISON.md): 15-epoch scatter
        # 11.4 m here vs 7.6 m for gnss-sdr with these two things on. pll_bw_narrow=0 disables.
        self.prn, self.fs = prn, float(fs)
        self.pll_bw_narrow, self.dll_bw_narrow = float(pll_bw_narrow), float(dll_bw_narrow)
        self.coh = int(coherent_ms) if pll_bw_narrow else 1
        self.bit_offset = None                 # period index (mod 20) at which a data bit begins
        self._flips = np.zeros(20, np.int64)   # sign-flip histogram for bit sync
        self._last_ip = 0.0
        self._acc = np.zeros(3, complex)       # E, P, L accumulated over the coherent window
        self._acc_n = 0
        self._acc_dt = 0.0
        self._aligned = False
        self._f_avg = None                     # carrier Doppler averaged over ~100 periods (stage 1)
        self._m2 = self._m4 = 0.0              # C/N0: running second and fourth moments of |P|
        self._mn = 0
        # acquisition hands over a code phase in SAMPLES (the sample at which the code starts);
        # the code phase in chips at sample 0 is therefore -(that) * chips/sample, modulo 1023
        chips_per_sample = CODE_RATE / self.fs
        cp0 = (-float(code_phase_samples) * chips_per_sample) % CODE_LEN
        self.s = ChannelState(prn=prn, fs=self.fs, carrier_hz=float(doppler_hz), code_phase=cp0)
        self.spacing = float(spacing)
        self.open_loop = bool(open_loop)
        self.T = CODE_LEN / CODE_RATE                       # nominal period, 1 ms
        self.pll_t1, self.pll_t2 = loop_gains(pll_bw, k=0.25)   # atan/2pi is +-1/4 cycle full scale
        self.dll_t1, self.dll_t2 = loop_gains(dll_bw)
        self.doppler0 = float(doppler_hz)
        self.carr_corr = 0.0                                 # PLL frequency correction, Hz
        self.code_dop = 0.0                                  # DLL rate correction, chips/s
        self.n_nominal = int(round(self.fs * self.T))      # samples per period, nominal
        self._ramp_cache, self._ramp_hz = None, 0.0
        self._epl = np.array([self.spacing, 0.0, -self.spacing])

    def _ramp(self, n):
        """exp(-2 pi j f t) for the current Doppler over n samples. The complex exponential is
        the receiver's most expensive line (109 us of a 170 us period, measured), so the ramp is
        cached and only rebuilt when the tracked Doppler has moved more than RAMP_TOL Hz: a
        0.5 Hz error over one period is 0.18 degrees, far under the loop's own noise. The phase
        continuity across periods is carried by carrier_phase, which is exact."""
        f = self.s.carrier_hz
        if self._ramp_cache is None or abs(f - self._ramp_hz) > self.RAMP_TOL or n > len(self._ramp_cache):
            m = max(n, self.n_nominal + 64)
            self._ramp_cache = np.exp(-2j * np.pi * f * (np.arange(m) / self.fs))
            self._ramp_hz = f
        return self._ramp_cache[:n]

    # ---- how many samples the next period needs ----------------------------------------
    def samples_needed(self):
        """Samples until the code phase reaches the next 1023-chip boundary at the current rate."""
        if self.open_loop:
            return self.n_nominal
        chips_left = CODE_LEN - self.s.code_phase
        return max(1, int(np.ceil(chips_left * self.fs / self.s.code_rate)))

    # ---- one code period ---------------------------------------------------------------
    def step(self, x):
        """x: exactly samples_needed() complex samples. Returns (I_p, Q_p) for the period and
        updates the state. Never allocates more than a few arrays of len(x)."""
        s = self.s
        n = len(x)
        # 1. carrier wipe-off. The NCO phase ramp is built once per period as a real vector; the
        #    complex exponential is the single most expensive line in the receiver.
        carr = self._ramp(n) * np.exp(-1j * s.carrier_phase)
        xb = x * carr
        # 2. code replicas at E / P / L as ONE (3, n) matrix; the six correlations are then a
        #    single matrix product with [Re, Im] - one BLAS call instead of six dots.
        chips_per_sample = s.code_rate / self.fs
        ph = s.code_phase + np.arange(n) * chips_per_sample
        if self.open_loop:
            prompt = code_at(self.prn, ph)
            ip, qp = float(np.dot(xb.real, prompt)), float(np.dot(xb.imag, prompt))
            ie = qe = il = ql = 0.0
        else:
            codes = code_at(self.prn, ph[None, :] + self._epl[:, None])          # (3, n): E, P, L
            c = codes @ np.column_stack((xb.real, xb.imag))                     # (3, 2)
            (ie, qe), (ip, qp), (il, ql) = (float(c[0, 0]), float(c[0, 1])), (float(c[1, 0]), float(c[1, 1])),                 (float(c[2, 0]), float(c[2, 1]))
        # advance phases by what this period consumed
        dt = n / self.fs
        s.carrier_phase = (s.carrier_phase + 2 * np.pi * s.carrier_hz * dt) % (2 * np.pi)
        s.carrier_cycles += s.carrier_hz * dt              # the NCO's own count: -lambda x this = range change
        s.code_phase = (s.code_phase + n * chips_per_sample) % CODE_LEN
        s.samples_in += n
        s.epochs += 1
        if not self.open_loop:
            # PLL lock indicator: cos(2 x phase error) = (I^2 - Q^2) / (I^2 + Q^2), averaged over
            # ~50 periods; +1 = all the energy in I, 0 = random phase. Sign-blind, so data bits
            # do not disturb it. (Kaplan & Hegarty §5.11 phase-lock detector.)
            p = ip * ip + qp * qp
            li = (ip * ip - qp * qp) / p if p > 0 else 0.0
            s.lock = 0.98 * s.lock + 0.02 * li
            # C/N0 by the moment method (gnss-sdr's 'SNR estimation for complex signals', Pauluzzi &
            # Beaulieu): over N prompts M2 = mean|P|^2, M4 = mean|P|^4; signal power Pd = sqrt(2 M2^2 - M4),
            # noise Pn = M2 - Pd; C/N0 = Pd / Pn / T. One number per 20 prompts, smoothed.
            self._m2 += p
            self._m4 += p * p
            self._mn += 1
            if self._mn == 20:
                m2, m4 = self._m2 / 20, self._m4 / 20
                pd = np.sqrt(max(2 * m2 * m2 - m4, 0.0))
                pn = m2 - pd
                if pd > 0 and pn > 0:
                    cn0 = 10 * np.log10(pd / pn / self.T)
                    s.cn0_db = cn0 if s.cn0_db == 0.0 else 0.9 * s.cn0_db + 0.1 * cn0
                self._m2 = self._m4 = 0.0
                self._mn = 0
            # bit sync (stage 1 only): where, mod 20, does the prompt's sign flip?
            if self.bit_offset is None and self.coh > 1:
                self._f_avg = s.carrier_hz if self._f_avg is None else 0.99 * self._f_avg + 0.01 * s.carrier_hz
                self._bit_sync(ip)
            # coherent window: a whole data bit once the edges are known, one period before
            in_window = self.bit_offset is not None
            if in_window:
                at_edge = (s.epochs - self.bit_offset) % self.coh == 0
                if not self._aligned:
                    # between bit sync and the first bit edge: no loop update at all. A partial
                    # first window with the error history reset kicked the narrow loop 18 Hz
                    # off on a synthetic sky (measured); the NCO free-runs a few ms instead.
                    self._aligned = at_edge
                    return ip, qp
                self._acc += (ie + 1j * qe, ip + 1j * qp, il + 1j * ql)
                self._acc_n += 1
                self._acc_dt += dt
                # the window ends on the period BEFORE a bit edge: the next period starts a bit
                if not at_edge:
                    return ip, qp
                (ie, qe), (ip_, qp_), (il, ql) = [(c.real, c.imag) for c in self._acc]
                dt_loop = self._acc_dt
                self._acc[:] = 0
                self._acc_n, self._acc_dt = 0, 0.0
            else:
                ip_, qp_, dt_loop = ip, qp, dt
            # 3. PLL: Costas discriminator (atan(Q/I): a data-bit sign flip does not move it),
            #    error in cycles, through the loop filter into a frequency correction; the NCO's
            #    phase accumulation (above) is the loop's integrator.
            # atan(Q/I), NOT atan2: the single-argument form is blind to a 180-degree data flip
            # (atan2 read a bit transition as a 165-degree phase error and slewed the carrier 100 Hz)
            e_pll = np.arctan(qp_ / ip_) / (2 * np.pi) if ip_ != 0 else 0.0     # cycles, (-1/4, 1/4)
            self.carr_corr += (self.pll_t2 * (e_pll - s.pll_e_prev) + dt_loop * e_pll) / self.pll_t1
            s.pll_e_prev = e_pll
            s.carrier_hz = self.doppler0 + self.carr_corr
            # 4. DLL: normalised early-minus-late envelope, error in chips, the same filter,
            #    carrier-aided (the carrier loop already knows the Doppler; this trims the residual)
            E, L = np.hypot(ie, qe), np.hypot(il, ql)
            e_dll = 0.5 * (E - L) / (E + L) if (E + L) > 0 else 0.0            # chips
            self.code_dop += (self.dll_t2 * (e_dll - s.dll_e_prev) + dt_loop * e_dll) / self.dll_t1
            s.dll_e_prev = e_dll
            s.code_rate = CODE_RATE + s.carrier_hz * CARRIER_TO_CODE + self.code_dop
        return ip, qp

    def _bit_sync(self, ip):
        """Stage 1 -> 2: after the loops have settled, histogram the prompt's sign flips mod 20;
        the bit edge collects them. When a bin stands 3x clear of the runner-up with at least 8
        flips, the edges are known: narrow the loops and start integrating over whole bits."""
        s = self.s
        # s.epochs has already counted this period, so THIS period's index is epochs - 1
        if self._last_ip != 0.0 and np.sign(ip) != np.sign(self._last_ip) and s.lock > 0.5:
            self._flips[(s.epochs - 1) % 20] += 1
        self._last_ip = ip
        if s.epochs >= 300 and self._flips.max() >= 8 and s.lock > 0.8:
            top = np.sort(self._flips)[::-1]
            if top[0] >= 3 * max(top[1], 1):
                # a flip counted for period k means period k is the first of a new bit; a window
                # therefore ends after period k-1, i.e. when epochs == k (mod 20)
                self.bit_offset = int(np.argmax(self._flips))
                # k=1 here, not the 0.25 of the 1 ms stage: the proportional path of this filter
                # form corrects (2 zeta wn / k) * T cycles per cycle of error per update, 1.7 at
                # T=20 ms with k=0.25 - over unity, and the loop tore itself apart on every
                # bandwidth tried (measured); with k=1 it is 0.4 and 5-15 Hz all hold
                self.pll_t1, self.pll_t2 = loop_gains(self.pll_bw_narrow, k=1.0)
                self.dll_t1, self.dll_t2 = loop_gains(self.dll_bw_narrow)
                # start the narrow loop from the AVERAGED frequency: the 1 ms loop jitters
                # +-5-10 Hz, and a 20 ms window can only pull in ~10 Hz (a handover that
                # landed at +10 Hz drove one synthetic bird to a false lock 18 Hz off, measured)
                if self._f_avg is not None:
                    self.carr_corr = self._f_avg - self.doppler0
                    s.carrier_hz = self._f_avg
                s.pll_e_prev = s.dll_e_prev = 0.0
                self._acc[:] = 0
                self._acc_n, self._acc_dt = 0, 0.0

    # ---- the observable ------------------------------------------------------------------
    def observable(self):
        """Everything the solver needs to make this channel a pseudorange once the Nav Decoder has
        anchored an epoch to a TOW.

        After step() the code phase is a little PAST the 1023-chip boundary (the period consumed
        whole samples; the crossing fell between two of them). The crossing itself - the instant
        the satellite's code epoch `epochs` arrived - is `epoch_sample` samples into the stream,
        fractional: samples_in minus the overshoot converted at the current code rate. Reported
        to the solver as one number, so no sign convention is left for a caller to get wrong.
        (Measured before this: 0.05-0.33 chips of overshoot read as a wrong-signed range error of
        up to 100 m per satellite, and the boundary rounded to a whole sample was 146 m more.)"""
        s = self.s
        overshoot_samples = s.code_phase * self.fs / s.code_rate
        return {"prn": self.prn, "epochs": s.epochs, "code_phase": s.code_phase, "carrier_hz": s.carrier_hz,
                "samples_in": s.samples_in, "epoch_sample": s.samples_in - overshoot_samples, "lock": s.lock,
                "cn0_db": s.cn0_db, "carrier_cycles": s.carrier_cycles - s.carrier_hz * overshoot_samples / self.fs}


def track_array(x, fs, prn, doppler_hz, code_phase_samples, n_ms, **kw):
    """Convenience: run a channel over an in-memory array for n_ms periods. Returns the prompt
    I+jQ per period and the channel."""
    ch = Channel(prn, fs, doppler_hz, code_phase_samples, **kw)
    out = np.empty(n_ms, np.complex128)
    pos = 0
    for k in range(n_ms):
        n = ch.samples_needed()
        if pos + n > len(x):
            out = out[:k]
            break
        ip, qp = ch.step(x[pos:pos + n])
        out[k] = ip + 1j * qp
        pos += n
    return out, ch
