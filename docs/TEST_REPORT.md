# gr-gpsrx test report

2026-09-22 · pure NumPy engines; the Channel block under GNU Radio 3.10.12 (radioconda)

## Tests (no capture, no radio): 23 engine tests + 7 flowgraph QA, all pass

(The table below is the original set; since then: `test_track.py` 7 - a Doppler-rate satellite
with the third-order loop; `test_pvt.py` 8 - RAIM isolation, Hatch smoothing, the timing product to
< 5 ns, a two-system solve recovering a 25 ns inter-system bias; `test_acquire.py` 2 - the LO
offset search; `test_gal.py` 3 - Galileo codes and BOC autocorrelation, the I/NAV decoder on real
PRN 29 symbols, a page round trip; `qa_channel_cc.py` 3 - the C++ Galileo E1 channel against the
Python one on two synthetic Galileo satellites: identical epochs, same secondary-code sync, E1-B
symbols out equal to the symbols that went in.)

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

## The plan of 22 September, run in order

**1. A live hour** (`apps/gpsrx_live.py`, RSPdx, attic antenna, everything on): 3601 s, **3525 fixes,
99.8% valid**, a median of 7 satellites, rms 1.9 m, 15-epoch scatter 0.6 m raw / 0.4 m filtered,
4.5 m spread over the whole hour, the hour's mean 5.8 m from the previous day's live fix; three
satellites set during the hour and one rose and was picked up; no thread or memory trouble.

**2. Cheap crystals.** An RTL-SDR-class LO can sit 47 kHz off at L1, outside a +-7 kHz search.
One wide coarse pass (500 Hz steps, 10 ms) finds the strongest satellites wherever they are; the
median of their Dopplers is the LO offset, and every later search is centred on it with a
widened window. Synthetic test: four satellites 33 kHz off, invisible to the normal search,
all found after the offset. Hardware test pending an RTL-SDR on the bench.

**3. Linux.** First build on a Raspberry Pi 5 (Debian, GNU Radio 3.10.12, g++ 14): clean; the 17
engine tests pass on aarch64/Python 3.13/numpy 2.2; the C++ QA passes and **eight C++ channels run
at 13.3x real time on the Pi**. GitHub Actions CI: the engine tests on ubuntu-latest, and a
from-scratch Ubuntu 24.04 build against the distribution's GNU Radio with the three QA files.

**4. Deterministic replay.** Three sources of run-to-run difference on the same file, removed:
the acquisition snapshot now starts on an exact sample (not a scheduler chunk boundary); each
channel starts on the exact sample the assignment names; and the PVT solves once per stream
second using, from every live channel, its observable at or before that second (not whichever
message happened to arrive first). Two runs of a drive leg: 78 fixes at identical instants,
identical satellite sets, **position difference 0.000000 m** - and the residual rms fell from
0.5 to 0.11 m as a side effect, because every observable in a solve now belongs to one instant.

**5. Sensitivity** (synthetic, one satellite, 40 s): tracking and decoding hold to **34 dB-Hz**
with the two-stage loops (single stage: decodes 4 subframes at 34, Doppler 13 Hz off); 31 dB-Hz
fails in the first stage (the 1 ms loops never lock). Found on the way: the two-stage handover
lost the code at 40 dB-Hz because the DLL's rate correction, +-1 chip/s of jitter at 1 ms, was
handed over unaveraged - now averaged like the carrier; and the switch was gated on lock > 0.8,
which weak signals never reach - now 0.5. An FLL for pull-in was tried (cross/dot discriminator)
and measured: no gain at any C/N0 or initial error; kept as an off-by-default option.

**6. The walkthrough** - `docs/WALKTHROUGH.md`, the receiver one block at a time, six figures from
synthetic satellites (`util/make_figures.py`). **7.** An announcement draft, not posted.

## The "not done" list, done (22 September, evening)

**Third-order PLL** for the narrow stage (Kaplan & Hegarty form), both engines, default on. Synthetic
satellite accelerating at 26 Hz/s (a car at 5 m/s^2): phase error 12 -> 7 degrees, lock 0.86 -> 0.96;
at 50 Hz/s 22 -> 7 degrees. Real air: the static attic unchanged to the decimal (nothing
accelerates); the 52 mph leg keeps five or more satellites in 42 fixes against 36. Test.

**The timing product.** Every fix carries GPS TOW at its receive sample, the week (mod 1024), the
sample clock's offset and drift (a 30 s fit) and the sample index of the next whole GPS second -
a 1PPS on the sample clock. The RSPdx's drift reads **-796 ppb; numpy-gps measured this TCXO at
+796.7 ppb in July** (opposite sign convention): the same number to half a ppb from two
receivers two months apart. Test: the solve returns GPS time to < 5 ns.

**Galileo E1.** Codes (gnss-sdr's table, GPL-3, vendored), a BOC(1,1) code function with the
subcarrier, quarter-chip correlators, 4 ms periods, one data symbol per period; the same Channel
engine with a `signal` description; the same Acquisition block with `system=GAL` (125 Hz steps);
a streaming I/NAV decoder (numpy-gps's Viterbi/deinterleaver/CRC/parser, incremental, page parts
decoded as they complete, anchors on the even part's first symbol); the solver with a fifth unknown
for the inter-system clock bias and each ephemeris's own mu. On the 4.096 MS/s wideband capture:
acquisition finds six Galileo satellites (the same three strongest numpy-gps decoded, plus three);
E1-B tracking locks at 0.97-0.98 with C/N0 42-44 dB-Hz; 11 of 11 page parts in 12 s, 5 pages
CRC-clean, 0 failures, anchors on the GST second grid, WN 1405; and **a joint fix from 11 satellites
(7 GPS + 4 Galileo), rms 4.1 m, PDOP 1.9, GST-GPS as this receiver sees it -111 ns**. Honestly:
on this capture the joint solve has lower residuals (6.5 -> 4.8 m median) and the same position as
GPS-only to 4 m, but MORE run scatter (8.9 vs 3.6 m) - the Galileo channels are first-stage only
(no bit-length integration exists on E1-B), joined the solve late, and the extra unknown costs a
degree of freedom. The gain will come from E1-C pilot tracking (100 ms coherent) - next. Galileo
needs >= 4 MS/s (BOC main lobes at +-1.023 MHz) and Python channels so far (the C++ twin is L1 C/A).
Tests: codes and BOC autocorrelation, the decoder on 12 s of real PRN 29 symbols, a page round trip
through the encoder chain, a two-system synthetic solve recovering a 25 ns bias to a nanosecond.

Found on the way: at 4.096 MS/s a code period (4096 samples; 16384 for Galileo) exceeds GNU
Radio's default 8191-item buffer, so the feeding block's output buffer is enlarged by the apps.

## Galileo E1-C pilot tracking (the same evening)

The loops now run on E1-C, whose 25-chip secondary code is known: the channel correlates its
prompt signs with the sequence at every offset (90% agreement names it), wipes the chips, and
integrates 20 ms with narrow loops - the same two-stage machinery as GPS - while a fourth correlator
on E1-B feeds the I/NAV decoder. Three real satellites: lock 0.97-0.98, every page CRC-clean, quieter
Doppler. 100 ms windows do not lock (a few Hz of post-handover error turns half a cycle inside the
window; a 20 -> 100 ms ramp would fix it - a static-antenna luxury, not done). A parameter lesson: the
window is milliseconds and is now converted to periods (a "20" that meant 80 ms at 4 ms periods
put the loop over unity gain).

The joint solve, settled (last 35 of 90 s, all 11 satellites in): **scatter 1.7-1.8 m, rms 4.6 m -
against GPS-only 3.8 m and 5.9 m** on the same samples. The whole-run scatter (7.6 m) is the transition
while Galileo joins; the fix's mean is 2.7 m from GPS-only's.

## The C++ Galileo channel (22 September, night)

The C++ Channel took the same `signal` generalisation as the Python engine: a code vector of any
length with an optional BOC(1,1) subcarrier and a settable spacing, an optional data code (a fourth
correlator, whose prompt goes to the decoder), an optional known secondary code (sync by correlating
the prompt's signs against it at every offset, then chip wipe-off), the coherent window in periods.
The E1-B and E1-C tables are compiled in (`lib/galileo_e1_codes.h`, generated from the vendored
table by `util/make_e1_codes_h.py`), so a GRC canvas or the live app runs Galileo with no Python
channel. `channel_cc(..., signal="E1")`; the Receiver uses it for its Galileo slots whenever the
engine is C++, which is what `apps/gpsrx_live.py --rate 4.096e6 --galileo 4` now does.

The same 90 s of the wideband capture, 8 GPS + 4 Galileo slots, C++ against Python channels:
**73 fixes each at identical instants, median position difference 0.07 m (max 1.68 m, during the
handover), settled scatter 1.78 vs 1.81 m, settled means 0.08 m apart, GST-GPS -115 vs -114 ns**;
E29 lock 0.97 at C/N0 45 dB-Hz, every page CRC-clean, all four with ephemerides. Wall time: 66 s
for the 90 s of 4.096 MS/s samples including the search holds (the Python channels: 587 s).
QA: two synthetic Galileo satellites (`synth.satellite(system="GAL")`: E1-B x symbols plus E1-C x
the secondary code, half the power each) through both engines - identical epoch counts, Doppler
within 2 Hz, code phase within 0.02 chip, the C++ symbol stream equal to the symbols that went in.

**Galileo LIVE** (the same evening, RSPdx, Antenna B, 4.096 MS/s, 8 GPS + 4 Galileo C++ channels,
8 minutes): four Galileo satellites tracked pilot-aided from the air (C/N0 33-37 dB-Hz, lock
0.80-0.92, every one decoding I/NAV to a complete ephemeris - 163-187 pages each); **426 fixes,
100% valid, 355 of them with Galileo in the solve; median 10 satellites (7 GPS + 3 Galileo), median
rms 2.2 m, PDOP 2.1, 15-fix scatter 0.6-0.7 m**, 2.6 m scatter over the Galileo era, the run's mean
3 m from the previous live fix and 7 m from the live hour's. The first Galileo satellite joined
at 4 minutes (the search runs at 125 Hz steps over 36 PRNs and the I/NAV ephemeris takes ~30 s
of clean pages). One number to keep watching: the inter-system bias read **-834 ns live against
-115 ns on the July wideband capture** from the same radio - the receiver-side part of GST-GPS
depends on the analog filter (the BOC lobes sit at +-1 MHz, where a 5 MHz IF filter's group delay
differs from the centre), and the two runs had different bandwidth settings. It is solved per
fix, so the position does not depend on it; the number itself is not a GGTO measurement.

## The pilot PLL, and the bug a Costas loop cannot see (22 September, late)

A reflection over the lab's earlier GPS notes found one thing numpy-gps did that this receiver did
not: a **pure four-quadrant PLL on the dataless pilot** (once the secondary code is wiped, E1-C
carries no sign flips, so `atan2(Q, I)` is legitimate: +-half a cycle of pull-in, no half-cycle
ambiguity). Switching to it made things WORSE - every Galileo channel lost lock within a second of
the handover, on the synthetic satellite and on the wideband capture alike. The trace showed the
discriminator at 0 for two windows and then +0.47 cycles, and again every few windows: the wipe was
applying the wrong chip on some periods. **The secondary-code sync had an off-by-one** (the last
entry of the sign history is the current period, whose count is `epochs`, not `epochs - 1`), so it
named an offset one short of the truth, and CS25 shifted by one disagrees on 12 of its 25 chips.
The Costas loop, sign-blind, had never noticed: it only lost coherent gain on those periods, and
tracked. The pure PLL is not blind, and failed loudly - which is the whole reason to prefer the
loop that cannot hide a sign error. Both engines fixed; a test asserts the offset the sync names
against the truth (`test_pilot_secondary_sync_and_pure_pll`).

With the wipe right and the pure PLL on: Galileo lock 0.62 -> 0.91 (E21, the weak one),
0.93 -> 0.99 (E29); the joint position unchanged (settled scatter 1.76 vs 1.78 m - it is code-
and GPS-limited). And the **100 ms pilot window now locks** - once the narrow PLL's bandwidth is
right for it: 15 Hz x 100 ms is a bandwidth-times-integration product of 1.5 (unstable; every
satellite lost in a second), 3 Hz x 100 ms is 0.3 (locks 0.93-0.99, settled scatter 1.71 m).
Both engines now cap the narrow bandwidth at 0.3 / window, so `--coherent-ms 100` works with the
defaults (result identical to the hand-set 3 Hz to 0.04 m); the GPS default, 15 Hz x 20 ms, is
exactly at the cap and is untouched.

**Found on the way - the Qt live flowgraph overruns at 4.096 MS/s.** Two attempts at a live
Galileo screenshot with the `gpsrx_radio_qt` flowgraph (spectrum display + Sky Panel) at
4.096 MS/s printed `O` (SoapySDR overflow) bursts, and each burst took every channel down: the
observable is a COUNT of code periods against the sample stream, so a dropped block silently
shifts every channel's clock, the solve goes wrong and the channels re-acquire. The headless
app at the same rate had no overflow in 8 minutes; the difference is the Qt frequency sink
running on the full 4 MS/s stream. Two things follow: the radio flowgraph's spectrum display
needs to be decimated (or removed) for Galileo rates, and the receiver should treat an
overflow as a resynchronisation event (drop every anchor) rather than carry on counting -
gr-soapy reports the overflow on the console, not as a stream tag, so this needs a watchdog
(a channel whose epoch count and `sample_abs` disagree with the sample rate by more than a
period) - done the same night, below. The Galileo screenshot in the README is from the replay,
where no sample is ever dropped.

## The stream watchdog: what a missing millisecond looks like (22 September, night)

The plan was the obvious one - compare each fix's sample-clock offset from GPS time with the last
fix's, flag a jump - and it is in (`pvt.samples_lost`, unit-tested). Then a capture with **exactly
one code period (2048 samples at 2.048 MS/s) cut out at 50 s** was replayed to see it fire, and it
did not: 49 valid fixes, rms 1.9 m, the offset on a perfectly straight -0.8 us/s line through the
gap. The lesson is worth the section. A gap of a whole number of code periods is invisible to the
loops (the code repeats) AND to the solver (the sample counter and every channel's epoch count skip
the same millisecond, so the stream stays self-consistent), and the nav decoders then RE-ANCHOR ON
THE SHIFTED BIT GRID - consistently, and 1 ms wrong in absolute time, for ever. The position is
right; the timing product is a millisecond off and nothing in the receiver can know. A gap that is
NOT a whole number of periods is the easy case: the DLL loses the code and the channels report
`lost` (the 2016-sample version of the same test: every channel lost at 50 s, re-acquired later).

The one thing a whole-period gap does move is where the data bits flip: one period earlier than
the bit grid says. So the channel's flip histogram now keeps running in stage 2 (`_monitor_grid`,
both engines; the pilot re-scores its secondary-code offset instead), and when another bin wins
clearly the channel moves its window, reports **`slip`** on its status port and puts a `gpsrx_slip`
tag on the prompt stream. PVT treats a slip like a new assignment (anchor and Hatch filter dropped,
timing history cleared); the Nav Decoder re-finds its bit grid on the tag without restarting its
period count (the I/NAV decoder already re-finds its page grid when the preamble moves). Found on
the way: the decoder's continuity check floored the period difference (`11999 // 6000 == 1`) and
refused every subframe after a slip - now rounded.

Result on the one-period-gap capture: all seven channels report the slip within 0.3-0.8 s of the
gap, the fixes pause 11 s (bit re-sync plus one subframe), and at 61.9 s they return **with the
clock offset stepped by exactly -1000 us** - the receiver's GPS time is right again - at rms
1-2 m. A clean 90 s replay and a live run show no false slips. Engine test: a synthetic satellite
with one period removed, the channel names the slip and the new grid one period earlier.

The live app also keeps the law's own accounting - the radio's item count against the wall clock,
printed every report. The count advances in scheduler chunks (+-100 ms of jitter, measured), so a
loss is read as a step in the deficit's floor, not from one reading.

## Scintillation indices (22 September, night)

Every channel now reports, per 60 s block, **S4** (the amplitude scintillation index: std/mean of
the coherent-window prompt power) and **sigma_phi** (the standard deviation of the PLL's residual
phase error, in radians - the phase jitter the loop did not track, i.e. high-passed at the loop
bandwidth; a proxy for the ICD's detrended sigma_phi, not the same thing). They travel in `obs`
(`s4`, `sigma_phi`) and print in the status table and the Sky Panel. The first version computed
S4 on 1 ms prompts and read 0.25 on a clean synthetic satellite at 45 dB-Hz: thermal noise alone
gives S4 = sqrt(1 / (C/N0 x T)), which is why the index is defined on 20 ms and why the standard
practice (Van Dierendonck 1993) subtracts the noise part from the measured C/N0; both are done,
in both engines, and the same synthetic satellite now reads below 0.1 while one whose amplitude
is modulated 1 + 0.5 sin(2 pi 0.7 t) reads its true S4 (0.35) to 0.06 (test).

Real air, the July attic capture (130 s, one full block per channel): **S4 0.07-0.24, sigma_phi
0.04-0.11 rad** across seven GPS and four Galileo satellites - a quiet mid-latitude ionosphere
seen through an attic roof, with the strongest satellite (G27, 48 dB-Hz) the most faded at 0.24:
multipath, not the sky. numpy-gps's July numbers on the same kind of capture (0.55-0.78) were 1 ms
prompts without the noise correction; these are the honest ones. What the index is for: the
space-weather probe - a disturbed ionosphere reads S4 above 0.3 on every high satellite at once,
which multipath never does.

## The decoded ephemerides against the IGS broadcast file (22 September, night)

The one check a synthetic round trip cannot give: the synthetic satellites carry messages made by
our own encoder, so a parser that mis-scales a field would agree with itself. The RINEX nav file
this receiver wrote from the July wideband capture was compared field by field with **BKG's IGS
broadcast ephemeris file for that day** (`BRDC00WRD_S_20262090000_01D_MN.rnx`, assembled from the
real-time streams), matched by PRN and time of ephemeris, both read by georinex
(`lab_local/eph_vs_broadcast.py`). **11 ephemerides, 7 GPS + 4 Galileo, 275 fields: every orbit,
clock, harmonic and rate term identical to the printed 12 digits**; the Galileo af0, unfolded from
the single-frequency BGD correction, equal to the broadcast value to 0.00 ns. The only two
differences are representation: one URA index that BKG prints as 2.8 m and we as 2^1.5 = 2.83 m,
and one Galileo health word where the IGS file carries an E5a flag that is not in I/NAV. Found on
the way and fixed before the comparison came clean: URA as metres by RTKLIB's convention, the
Galileo health bits packed the RINEX way, SISA carried into the ephemeris.

## An elevation mask, and what it cost (22 September, night)

`mask_deg` on the PVT (`--mask-deg` on the apps): after a first solve gives the elevations,
satellites below the mask are left out. The first version kept as few as four. On the attic
capture every satellite sits above 21 degrees and a 10 degree mask changed nothing (0.15 m). On
the Smith Point capture - five satellites, one at 8 degrees - **the mask took the fix from five
satellites to four: residual rms 0 by construction, the 15-fix scatter 0.9 -> 4.6 m, and the
mean moved 26 m.** The low satellite carried more redundancy than error. So the mask now keeps
at least five (six with two systems) - RAIM's spare - and on that sky it leaves the low satellite
in. The weighting already demotes low satellites; the mask is for skies with satellites to spare.

## RINEX 3 output (22 September, night)

`--rinex-file run.obs` on either app (or `rinex_file` on the PVT / Receiver blocks) writes a
RINEX 3.04 observation file - C1C, L1C, D1C, S1C for GPS (G) and Galileo (E) at every valid fix,
epochs in GPS time - and `run.nav`, a mixed navigation file with every decoded ephemeris (GPS LNAV
with IODE/IODC/TGD/URA/health, now parsed for it; Galileo I/NAV with the BGDs and the af0 unfolded
to the broadcast value). The observables are the raw ones: the code pseudorange with SV clock,
ionosphere and troposphere still in it, referred to the solve's receive instant; the carrier phase
in cycles increasing with range; Doppler positive closing. No position is written (APPROX POSITION
is zero on purpose).

Checked two ways. Unit: GPS time to calendar at the 1980, 1999 and 2019 epochs; a write/read round
trip of the synthetic constellation to the millimetre; both files read by **georinex** (an
independent parser) with every satellite and field present. Real air: the 90 s wideband replay
wrote 61 epochs x 11 satellites and 11 ephemerides; georinex read both; our solver, fed **the
file's numbers alone** (raw code, no smoothing, no ionosphere), solved 61 of 61 epochs with a mean
1.6 m from the receiver's own fixes. The check that matters - RTKLIB's `rnx2rtkp` on our files -
waits for the Pi (Debian's `rtklib` package); neither Linux machine was on tonight.

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
