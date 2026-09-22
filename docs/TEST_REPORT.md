# gr-gpsrx test report

2026-09-22 · pure NumPy engines; the Channel block under GNU Radio 3.10.12 (radioconda)

## Tests (no capture, no radio): 13 pass

| file | proves |
|---|---|
| `test_track.py` (6) | C/A code vs the published chips; **gate 0:** open-loop prompts equal numpy-gps's `prompts_ms()` to 1e-9; closed loop pulls in from 40 Hz / 0.6 chip off and recovers the data bits; the epoch count is a consistent clock; eight satellites over one noise floor all lock (PLL lock > 0.85); throughput number |
| `test_nav.py` (3) | parity round trip and flipped-bit detection; ephemeris round trip to one LSB per field from clean bits; ephemeris + timing anchors 6000 periods apart from a tracked synthetic satellite carrying a real message |
| `test_pvt.py` (4) | six-satellite synthetic constellation round trip to < 0.5 m; transmit time from the epoch count and the common receive instant; the epoch's fractional arrival sample to < 0.15 chip; the size of the clock-order error |
| `qa_channel.py` (3, GNU Radio) | idle channel; eight channel blocks lock on a synthetic sky and publish observables; a lost signal is reported and the channel idles |

## Gate 1: throughput (decides Python vs C++ for the Channel)

Eight Python channel blocks in the GNU Radio scheduler: 1 channel 3.1x real time, 2: 1.1x,
4: 0.40x, **8: 0.17x**. Eight are slower than eight times one - the interpreter lock, not the
arithmetic. Batching periods per `work()` and caching the NCO ramp did not help; the engine
alone single-threaded is 0.75x. **The Channel goes to C++.**

## Real air: two captures, one attic antenna, one evening

| capture | satellites | subframes / channel | fix |
|---|---|---|---|
| 120 s | 8 tracked, 8 decoded | 18 | residual rms 3.6 m, PDOP 2.5, 15 epochs scatter 6.5 m, archived ionosphere |
| 240 s, 2 h later | 7 tracked, 7 decoded | 38 | **residual rms 0.3 m**, PDOP 2.4, 15 epochs scatter 9.3 m, ionosphere decoded live (page 18) |

Against numpy-gps on the SAME captures, same antenna, same spot (offsets only; no coordinates):

| | ours | numpy-gps |
|---|---|---|
| repeatability, 240 s fix vs 120 s fix | **9.9 m** (up +8.5) | 102 m (up -98: its two altitudes differ by 100 m) |
| residual rms | 0.3 / 3.6 m | 12 / 26 m |
| the two receivers, same capture | 150 m apart (E +90, N -44, U +112) | |

Cross-check that separates solver from measurements: numpy-gps's own cached pseudoranges through
THIS solver land 9 m from numpy-gps's fix. The solvers agree; the difference is in the
pseudoranges, and the receiver whose fixes repeat to 10 m and whose residuals are 0.3 m is the
one with the better pseudoranges. The remaining question - which is *true* - needs an
independent position and is left to the owner.

## Defects found by testing (all fixed)

1. Costas discriminator written as `atan2(Q, I)`: a 180-degree data flip read as a 165-degree phase error, the carrier slewed 100 Hz, every bit transition glitched. Must be `atan(Q/I)`. (Two hours.)
2. Test synthetic sky built by summing single-satellite signals summed their noise floors too; `synth.sky` puts them over one.
3. Channel block: with a 1 kHz output the scheduler called `work()` once per period; `set_output_multiple` demands whole batches or nothing; `set_min_input_buffer` is not exposed to Python.
4. Nav encoder (test tool) had subframe 1's clock words one word early; the round trip caught it.
5. Acquisition's 250 Hz Doppler bins are wider than the PLL's pull-in: two satellites failed to lock until a 100 ms FFT refinement to ~10 Hz was added.
6. The observable reported the code phase at the period's END (0.05-0.33 chip past the epoch) and the boundary on a whole sample: up to 250 m of range error, fixes wandering 190 m between epochs. Now: the epoch's arrival as a fractional sample. (190 m -> 6.5 m scatter.)
7. Klobuchar evaluated at a nonsense elevation while the solve was still far from Earth returned NaN and the least squares "did not converge". Atmosphere now applied only near the surface, above the horizon, if finite.
8. A parity-clean false frame (ten words by chance) gave one satellite a TOW past the end of the week and a residual of 10^13 m. Frames now need TOW < 604800 and continuity (6 s per 6000 periods, ids stepping 1-5) with the previous one.
