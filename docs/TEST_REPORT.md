# gr-gpsrx test report

2026-09-22 · pure NumPy engines; the Channel block under GNU Radio 3.10.12 (radioconda)

## Tests (no capture, no radio): 14 engine tests + 4 flowgraph QA, all pass

| file | proves |
|---|---|
| `test_track.py` (6) | C/A code vs the published chips; **gate 0:** open-loop prompts equal numpy-gps's `prompts_ms()` to 1e-9; closed loop pulls in from 40 Hz / 0.6 chip off and recovers the data bits; the epoch count is a consistent clock; eight satellites over one noise floor all lock (PLL lock > 0.85); throughput number |
| `test_nav.py` (3) | parity round trip and flipped-bit detection; ephemeris round trip to one LSB per field from clean bits; ephemeris + timing anchors 6000 periods apart from a tracked synthetic satellite carrying a real message |
| `test_pvt.py` (4) | six-satellite synthetic constellation round trip to < 0.5 m; transmit time from the epoch count and the common receive instant; the epoch's fractional arrival sample to < 0.15 chip; the size of the clock-order error |
| `test_acquire.py` (1) | four synthetic satellites found at their code phase (< 1 sample) and Doppler (< 15 Hz after refinement), no absent PRN reported, the metric separates present from absent by > 1.5x |
| `qa_receiver.py` (1, GNU Radio) | the whole flowgraph on a synthetic sky: Acquisition finds all four and assigns them, every Channel locks and reports an absolute observable, every Nav Decoder anchors on the 6 s grid with 6000 periods between subframes |
| `qa_channel_cc.py` (2, GNU Radio) | the C++ Channel against the Python one on the same synthetic sky: identical epoch counts, code phase within 0.02 chip, Doppler within 2 Hz, data bits agree > 99%; **gate 1 in C++: eight channels x 4 s in 0.25 s = 15.7x real time** |
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

## The flowgraph on real air (same 240 s capture)

| run | channels | result |
|---|---|---|
| `apps/gpsrx_replay.py` (Receiver hier block), 60 s | 8 | first fix at 7 s search + 36 s of stream; 7 satellites, rms 2.5 m, PDOP 2.4, 15-fix scatter 10.8 m, iono decoded from the stream; **9.9 m from the offline engine's fix** (E +3.9, N +4.2, U -8.1) |
| `build/grc/gpsrx_canvas.py` (grcc output, unedited), 240 s | 4 | 203 fixes; 4 satellites exactly -> rms 0 by construction, PDOP 17 (the four strongest cluster in the sky), scatter 45-77 m: geometry, so the canvas example now carries six channels |
| `build/grc/gpsrx_replay_qt.py` + Sky Panel, 330 s wall | 8 | 7 satellites, rms 2.2 m, PDOP 2.4, scatter 13.9 m over 11; the panel draws the sky from the fix's az/el and shows no coordinate |

Acquisition metric on air, 100 ms non-coherent: the seven real satellites 3.7-22.7, the best absent
PRN 1.7; threshold 2.5. Search 6.7-7.0 s on one core.

Replay lesson: a file source runs as fast as its readers, and idle channels read fast - 30 s of
capture streamed past during the first 7 s search and were never tracked. The Acquisition block's
`hold` stops the stream while a search runs (replay only; live, the radio paces it).

## The flowgraph with C++ channels, and LIVE

| run | result |
|---|---|
| Receiver (C++ channels), whole 240 s capture | 90 s wall including four 7 s search holds; 179 fixes (one per stream second); rms 2.5 m, PDOP 2.4, 15-fix scatter 9.5 m; **1.9 m from the offline engine's fix** on the same file |
| **`build/grc/gpsrx_radio_qt.py` (grcc output, unedited) on an SDRplay RSPdx, active antenna, 200 s** | 8 satellites found in a 7.4 s search, all eight tracked (lock 0.80-0.96), 31 subframes each; **fix from 8 satellites, rms 2.6 m, PDOP 2.9, 15-fix scatter 9.8 m, ionosphere decoded live** |
| live fix vs the offline fix on the earlier recording | **7.2 m** (E +6.0, N -1.3, U -3.7): a different constellation, hours apart, the same antenna |

That last number is the strongest evidence on the 150 m question above: a receiver whose
pseudoranges carried a systematic error would land somewhere else with a different geometry.
Two skies agreeing to 7 m says the error is small. numpy-gps's two fixes on the earlier
captures differ from each other by 102 m. The owner's known coordinates remain the final word.

## The 150 m question, answered - and numpy-gps fixed

The owner confirmed gr-gpsrx's position. Comparing numpy-gps's cached per-satellite code phases
with this receiver's tracking channels at the same instant of the same capture: every satellite's
range was off by -0.15 s x its range rate (-114 m at +4.5 kHz of Doppler, +20 m at -1 kHz; slope
0.93 to that prediction, 6.8 m residual). numpy-gps measures the code phase from a 300 ms
non-coherent snapshot, through which the code drifts up to a chip, so the summed peak sits at the
phase of the window's CENTRE while the transmit time was evaluated at its START. One line
(attribute the phase to the centre) and numpy-gps lands **1.7 m** from this receiver's live fix
(was 152 m); residual rms 12 -> 6.6 m. Its 26 m scatter against this receiver's 10 m is the
snapshot-versus-tracking-loop difference, not a defect.

The same fix on a waterside field capture from two weeks earlier (five satellites low in the sky,
PDOP 3-4): numpy-gps moved from 68 m to **13.5 m** from this receiver's mean of 229 fixes, its
residual rms 41 -> 5 m, and the "creek bias" in its height (-58 m; sea level is about -33 m on
the ellipsoid there) was the same bug. This receiver: median rms 3.3 m, scatter 27 m at that site.
A 90 s capture from the next day gives no fix in either receiver (never four satellites with
ephemeris and timing at once).

That field capture found a defect here too: a channel that lost its satellite at 134 s kept its
last observable in the solver, and every later fix carried it (rms 3 -> 50 m). PVT now drops a
channel's observable on 'lost'/'idle' and ignores any observable more than 2.5 s older than the
newest; the fix file keeps every fix's quality and ECEF as a history.

## A drive: 27 x 90 s from a car roof, up to 24 m/s (54 mph), both receivers on the same IQ

Cold (nothing carried between captures): gr-gpsrx fixed 23 of 27. The misses had seven or eight
satellites acquired strongly but only three ephemerides finished - at speed, fades break subframes
and 90 s is not long enough to decode everything from scratch. numpy-gps fixed those because it
carries orbits between captures, so PVT now keeps every decoded ephemeris in a file and lends it to
a channel that has timing but has not finished its own decode (toe within 2 h). **Warm: 27 of 27**,
37-81 fixes per capture at 1 Hz, residual rms 1-3 m, 14 captures using a borrowed orbit.

Against numpy-gps on the same captures (offsets between the two receivers' tracks at matching
instants, no coordinates): numpy-gps AFTER its window-centre fix sits a median **28 m** from this
receiver's track and agrees on the car's speed to ~1 m/s on every moving leg; numpy-gps as run ON
THE DRIVE NIGHT (before the fix) was a median **217 m** off, up to 1 km, with one 2,500 km "fix" -
the window-centre error scales with Doppler, and a car roof sees the full +-5 kHz. numpy-gps gave no
fix on three captures this receiver fixed, and its two four-satellite solutions were 100-700 m out
(no redundancy). The owner's eyes on a private map of both tracks are the truth check for the road.

## Defects found by testing (all fixed)

1. Costas discriminator written as `atan2(Q, I)`: a 180-degree data flip read as a 165-degree phase error, the carrier slewed 100 Hz, every bit transition glitched. Must be `atan(Q/I)`. (Two hours.)
2. Test synthetic sky built by summing single-satellite signals summed their noise floors too; `synth.sky` puts them over one.
3. Channel block: with a 1 kHz output the scheduler called `work()` once per period; `set_output_multiple` demands whole batches or nothing; `set_min_input_buffer` is not exposed to Python.
4. Nav encoder (test tool) had subframe 1's clock words one word early; the round trip caught it.
5. Acquisition's 250 Hz Doppler bins are wider than the PLL's pull-in: two satellites failed to lock until a 100 ms FFT refinement to ~10 Hz was added.
6. The observable reported the code phase at the period's END (0.05-0.33 chip past the epoch) and the boundary on a whole sample: up to 250 m of range error, fixes wandering 190 m between epochs. Now: the epoch's arrival as a fractional sample. (190 m -> 6.5 m scatter.)
7. Klobuchar evaluated at a nonsense elevation while the solve was still far from Earth returned NaN and the least squares "did not converge". Atmosphere now applied only near the surface, above the horizon, if finite.
8. A parity-clean false frame (ten words by chance) gave one satellite a TOW past the end of the week and a residual of 10^13 m. Frames now need TOW < 604800 and continuity (6 s per 6000 periods, ids stepping 1-5) with the previous one.
9. Replay raced past the receiver during the search (above): `hold`.
10. A code phase found seconds earlier was wrapped forward with the NOMINAL code period; at 5 kHz of Doppler the code runs 3 chips/s fast, 22 chips over a 7 s search. Wrapped with the Doppler-shifted period.
11. The Channel block's `epoch_sample` was relative to its own engine start, not the flowgraph's sample clock: fine for one channel, wrong for a solve across eight. Made absolute.
12. The GNU Radio Python gateway looks message handlers up by NAME: a lambda handler raised `no attribute '<lambda>'` on a scheduler thread and the test saw an empty sink.
13. An idle C++ channel consumed two seconds of test samples faster than the test's sleep could post the assignment: post before start (messages queue). Also true of any fast block.
14. PVT throttled on the wall clock: one fix per wall second, five per second of capture at 5x replay. The receiver's clock is the sample counter; throttle on that.
15. gr-soapy's per-channel `settings` cannot carry SDRplay's bias-T (`biasT_ctrl` is a device-level setting): `ValueError: Unsupported setting`. Applied with the block's device-level `write_setting`.
16. The extension module copied beside the package made plain CPython's `import gpsrx` raise `ImportError` (GNU Radio's DLLs not loadable there), which the `ModuleNotFoundError` guard let through. Guard on `ImportError`.
17. A lost channel's last observable stayed in the solver (above).
