# gr-gpsrx against gnss-sdr: head-to-head, and what to borrow

2026-09-22. gnss-sdr 0.0.20 (Debian package, on a Raspberry Pi 5) and gr-gpsrx (C++ channels, on
the same IQ files). Offsets and scatter only; no coordinates.

## The head-to-head

**A 240 s static capture (attic antenna, 7 satellites).** Reference = gr-gpsrx's live fix from a
different constellation hours later (confirmed by the owner against the antenna's position).

| receiver / configuration | epochs | mean vs reference | scatter over the run | median 15-epoch scatter |
|---|---|---|---|---|
| gr-gpsrx (as is) | 180 | 7.3 m | 14.2 m | 11.4 m |
| gnss-sdr default (`DLL_PLL`, 1 ms, RTKLIB single) | 145 | 4.8 m | 13.0 m | 9.5 m |
| gnss-sdr, tracking tuned (20 ms coherent, PLL 5 Hz / DLL 0.5 Hz narrow) | 173 | 3.6 m | 15.5 m | **7.6 m** |
| gnss-sdr, tuned + Kalman PVT + RAIM | 173 | 4.6 m | 5.4 m | **1.6 m** |

gnss-sdr's mean and gr-gpsrx's mean on this capture are 3.7 m apart. All three configurations
and gr-gpsrx put the antenna in the same place; the differences are per-epoch noise (tracking)
and filtering (PVT).

**Three 90 s legs of a drive, from the car roof.** Both receivers drawn on the road at zoom 18:
gr-gpsrx and gnss-sdr lie on top of each other in the lane, 2-3 m apart laterally, on both moving
legs (~52 mph); at the stop they agree to 6 m. (The 20 m "median distance at matched times" on
the moving legs is one second of along-track timestamp alignment, not position.)

**Verdict.** On the same air, the two receivers give the same answer. gnss-sdr's per-epoch
measurements are ~30% quieter with its tracking tuned (7.6 m vs 11.4 m at 15 epochs), and its
Kalman PVT turns that into 1.6 m for a static user. gr-gpsrx is ~3k lines you can read; gnss-sdr
is ~200k. The gap is specific and measured, and closable.

## What gnss-sdr does at each stage, against ours

**Acquisition.** PCPS (parallel code phase search), same FFT search as ours, `max_dwells`
non-coherent dwells, threshold on peak/second-peak like ours, a channel manager that assigns
free channels on acquisition events (ours: the allocator inside the Acquisition block). Same.

**Tracking (`GPS_L1_CA_DLL_PLL_Tracking`).** Costas PLL (two-quadrant atan - ours too), a
2nd/3rd-order loop filter (`pll_filter_order`, default 3; ours: 2nd), a non-coherent early-minus-late
envelope DLL (ours too), 0.5-chip spacing (ours too). **Three things ours lacks:**
1. **Two-stage bandwidths.** Wide loops for pull-in (`pll_bw_hz` 50, `dll_bw_hz` 2), then after bit
   synchronisation **narrow** loops (`pll_bw_narrow_hz` 20 default, 5 in our tuned run;
   `dll_bw_narrow_hz`). Ours runs 18 Hz / 2 Hz forever.
2. **Extended coherent integration** (`extend_correlation_symbols`, up to 20 ms, must divide the
   20 ms bit) once the bit edges are known: the correlators sum over a whole data bit, 13 dB
   more coherent gain, and the discriminators run 20x less often with 20x less noise each.
3. **An FLL for pull-in** (`enable_fll_pull_in`, `fll_bw_hz` 35, `pull_in_time_s` 5): a frequency
   discriminator (cross/dot product) that pulls a 100 Hz-off carrier in where a PLL would cycle-slip;
   ours relies on the 100 ms Doppler refinement in acquisition instead (which works, but is not a loop).

Also: a **C/N0 estimator** (second/fourth moments of the prompt, `cn0_smoother_alpha`), a carrier
lock detector cos(2Δφ) > `carrier_lock_th` 0.85 (ours: the same statistic, threshold 0.2 over 500
periods), and loss-of-lock on `cn0_min` 25 dB-Hz / `max_lock_fail` 50. Ours has no C/N0 number.

**Telemetry.** Preamble + parity framing, TOW from the HOW, one symbol counter per channel.
Ours: the same, plus continuity rejection of parity-clean false frames.

**Observables (`Hybrid_Observables`).** Common receive time by choosing the channel with the most
recent TOW and sliding the others to it ("without waiting for a particular bit front on each
channel") - the same idea as our `refer_to_common_sample`. Transmit time from TOW plus the
symbol count - the same counted-epoch idea as ours (theirs counts symbols, ours code periods).
**Plus carrier phase and Doppler as observables**, and **carrier smoothing of the pseudorange**
(`enable_carrier_smoothing`, Hatch filter, `smoothing_factor` 200). Ours has neither carrier phase
nor smoothing; numpy-gps has a Hatch filter (`hatch.py`) that was measured at 11.8 m joint scatter.

**PVT (`RTKLIB_PVT`).** Weighted least squares with an elevation-dependent variance per satellite
(a + b / sin(el), plus ephemeris, iono, tropo and code-bias variances; optional C/N0 de-weighting
`error_factor_snr`); **RAIM/FDE** (`raim_fde`: chi-square test, drop one satellite at a time, needs
six); Klobuchar or SBAS ionosphere, Saastamoinen troposphere (ours: Klobuchar + a simple
mapping); `elevation_mask` 15 (ours: 0.02 rad ~1); and the **Kalman PVT** (`enable_pvt_kf`:
state [x y z vx vy vz], constant-velocity model, process/measurement noise as parameters). Ours:
unweighted least squares, no RAIM, a running mean of the last 15 fixes.

**What ours has that theirs does not** (and should keep): readable Python for every block with a
C++ twin that is port-for-port identical; the observable as one number (the epoch's fractional
arrival sample) with no millisecond ambiguity anywhere; every intermediate as a message on the
canvas; the position never printed, a screenshot-safe panel; a same-capture gate against an
independent oracle; warm start from a plain JSON file.

## The borrow list, ranked by measured gain per line of code

| # | borrow | where | expected gain | measured (theirs) |
|---|---|---|---|---|
| 1 | two-stage loop bandwidths after bit sync | Channel (py + C++) | most of the 11.4 -> 7.6 m | yes |
| 2 | 20 ms coherent integration after bit sync | Channel (py + C++) | the rest of it, and weak-signal margin | yes |
| 3 | C/N0 estimator | Channel | weighting, honest lock/loss decisions | - |
| 4 | elevation + C/N0 weighted LS | PVT | fewer bad-satellite pulls (the drive's 4-5 sv legs) | (RTKLIB) |
| 5 | RAIM/FDE | PVT | drop the one bad pseudorange that made rms 38 m on a drive leg | (RTKLIB) |
| 6 | carrier-phase observable + Hatch smoothing | Channel + PVT | pseudorange noise / 5-10 on a static user | (numpy-gps: 11.8 m) |
| 7 | Kalman PVT (constant velocity) | PVT | 7.6 -> 1.6 m static; a smooth car track | yes |
| 8 | FLL pull-in | Channel | robustness at large Doppler error; not needed today | - |

Both are GPL-3, so their code may flow here - but every item above is a few dozen lines when
written against our engine, and written that way it stays readable, which is the point of this
receiver. The order is the order to build.
