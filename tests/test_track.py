# SPDX-License-Identifier: GPL-3.0-or-later
"""The tracking engine, with no GNU Radio and no capture.

Gate 0 (correlator identity) is `test_open_loop_prompts_equal_numpy_gps`: it needs numpy-gps on
disk (NUMPY_GPS_DIR) and is skipped otherwise; everything else runs anywhere."""
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from gpsrx import cacode, synth  # noqa: E402
from gpsrx.track import Channel, track_array  # noqa: E402

FS = 2.048e6


def test_ca_code_matches_the_published_first_chips():
    assert cacode.selfcheck()
    c = cacode.ca_code(7)
    assert c.shape == (1023,) and set(np.unique(c)) == {-1.0, 1.0}
    assert abs(np.dot(c, np.roll(c, 1))) < 70              # a Gold code: tiny off-peak autocorrelation


@pytest.mark.skipif(not os.environ.get("NUMPY_GPS_DIR"), reason="gate 0 needs numpy-gps (NUMPY_GPS_DIR)")
def test_open_loop_prompts_equal_numpy_gps():
    """GATE 0: open-loop prompts == numpy-gps prompts_ms(), the same 1 ms sums, to float precision."""
    sys.path.insert(0, os.environ["NUMPY_GPS_DIR"])
    import measure as gt
    x, tr = synth.satellite(11, FS, 0.2, doppler_hz=-2500.0, code_phase_samples=1234, cn0_dbhz=44.0)
    xd = x.astype(np.complex128)
    ref = gt.prompts_ms(xd, FS, 11, tr["doppler_hz"], tr["code_phase_samples"], 150)
    ours, _ = track_array(xd, FS, 11, tr["doppler_hz"], tr["code_phase_samples"], 150, open_loop=True)
    # numpy-gps rolls a zero-phase sampled code by the code phase; we generate the code at a
    # fractional chip phase. Same nearest-chip model, so the sums must agree to rounding.
    assert np.allclose(ours, ref, rtol=1e-9, atol=1e-6), np.max(np.abs(ours - ref))
    # ...and the correct PRN correlates far above a wrong one: the synthetic bird is really there.
    # At C/N0 44 dB-Hz a 1 ms coherent sum sits ~5x above the noise floor (44 - 30 = 14 dB).
    wrong, _ = track_array(xd, FS, 12, tr["doppler_hz"], tr["code_phase_samples"], 150, open_loop=True)
    assert np.median(np.abs(ours)) > 3.5 * np.median(np.abs(wrong))


def test_closed_loop_pulls_in_and_recovers_the_data_bits():
    x, tr = synth.satellite(21, FS, 1.0, doppler_hz=1234.5, code_phase_samples=700, cn0_dbhz=45.0)
    # start 40 Hz and 0.3 chip off the truth, as acquisition would
    prompts, ch = track_array(x.astype(np.complex128), FS, 21, tr["doppler_hz"] + 40.0,
                              tr["code_phase_samples"] + 0.6, 1000)
    assert len(prompts) == 1000
    s = ch.s
    assert abs(s.carrier_hz - tr["doppler_hz"]) < 2.0, s.carrier_hz          # PLL converged
    assert abs(s.code_rate - tr["code_rate"]) < 0.5, s.code_rate               # DLL converged
    assert s.lock > 0.8                                     # PLL lock indicator, smoothed over ~50 ms
    # after pull-in the prompt I channel carries the bits, 20 periods each. The tracker's epoch
    # grid is offset from the signal's by however many periods acquisition's handover cost (here
    # one), so vote each 20-period group at every candidate offset and take the best: a bit-sync
    # block does exactly this. The PLL takes ~3 periods to re-settle after a 180-degree data
    # transition, so the vote must tolerate a few wrong periods per bit.
    I = prompts.real[400:]
    truth = tr["bits"]
    best = 0.0
    for off in range(20):
        seg = I[off:]
        m = len(seg) // 20
        votes = np.sign(seg[:m * 20].reshape(-1, 20).mean(axis=1))
        for k0 in range(18, 23):                             # the handover cost 0..2 periods
            want = truth[k0:k0 + m]
            if len(want) == m:
                best = max(best, abs(float(np.mean(votes * want))))
    assert best > 0.97, best


def test_epoch_count_and_phase_are_a_consistent_clock():
    """Every step consumes exactly one code period: after the first boundary, the samples
    consumed per epoch must equal 1023 chips at the tracked code rate, to within one sample."""
    x, tr = synth.satellite(3, FS, 0.5, doppler_hz=800.0, code_phase_samples=100)
    ch = Channel(3, FS, 800.0, 100)
    xd = x.astype(np.complex128)
    pos, marks = 0, []
    for _ in range(400):
        n = ch.samples_needed()
        ch.step(xd[pos:pos + n])
        pos += n
        marks.append((ch.s.epochs, ch.s.samples_in, ch.s.code_rate))
    e0, s0, _ = marks[100]
    e1, s1, rate = marks[399]
    t_chips = (e1 - e0) * cacode.CODE_LEN / rate
    t_samples = (s1 - s0) / FS
    assert abs(t_chips - t_samples) < 1.0 / FS * 2, (t_chips, t_samples)


BIRDS = [dict(prn=p, doppler_hz=d, code_phase_samples=c, cn0_dbhz=44.0) for p, d, c in zip(
    (1, 3, 7, 11, 14, 17, 22, 30), (-3210, 1234.5, 2800, -900, 410, -2200, 1750, -150),
    (100, 700, 1500, 2000, 300, 1200, 900, 50))]


def test_eight_satellites_track_at_once():
    """A real sky: eight birds over one noise floor, eight channels, every one locks."""
    x, truths = synth.sky(FS, 1.5, BIRDS)
    xd = x.astype(np.complex128)
    for b in BIRDS:
        ch = Channel(b["prn"], FS, b["doppler_hz"] + 30.0, b["code_phase_samples"] + 0.4)
        pos = 0
        for _ in range(1400):
            n = ch.samples_needed()
            ch.step(xd[pos:pos + n])
            pos += n
        assert ch.s.lock > 0.85, (b["prn"], ch.s.lock)


def test_throughput_eight_channels_real_time():
    """GATE 1: eight LOCKED channels on one second of samples. The bar for live use is 1.0x real
    time; this prints the number and only fails below 0.5x (the hard gate is in the report)."""
    x, truths = synth.sky(FS, 1.0, BIRDS)
    xd = x.astype(np.complex128)
    chans = [Channel(b["prn"], FS, b["doppler_hz"], b["code_phase_samples"]) for b in BIRDS]
    t0 = time.perf_counter()
    for ch in chans:
        pos = 0
        while True:
            n = ch.samples_needed()
            if pos + n > len(xd):
                break
            ch.step(xd[pos:pos + n])
            pos += n
    wall = time.perf_counter() - t0
    print(f"\n8 locked channels x 1 s of samples: {wall:.2f} s CPU = {1.0 / wall:.2f}x real time")
    assert wall < 2.0
