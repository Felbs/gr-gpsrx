# gr-gpsrx design: a GPS L1 C/A receiver made of GNU Radio blocks

> **Status 2026-09-22:** built as designed, in Python: all six blocks plus a Receiver hier block,
> four flowgraphs, the tests in the gates section; gate 0 passed, gate 1 failed (0.17x for eight
> Python channels -> the Channel goes to C++), gate 2 superseded by a finding (this receiver's
> pseudoranges repeat better than the reference's; see TEST_REPORT), gate 3 passed on air
> (rms 2.2-3.6 m), gate 4 (live) not yet run. Differences from the text below: the allocator
> lives inside the Acquisition block (no separate Channel Bank block; the Receiver hier block
> plays that role), and the observable carries the epoch's fractional arrival sample rather
> than a whole-sample boundary plus code phase.

Status: built (see the note above). Public: https://github.com/Felbs/gr-gpsrx

## Purpose

Educational first, working second, and the two are not in tension: a receiver whose every
stage is a block on the canvas, whose every intermediate is a message you can plot, and which
produces a real position from a real antenna. The reference and oracle is numpy-gps, the
offline pure-NumPy receiver: every algorithm here already exists there and has been proven on
air. This project is not a port of numpy-gps into blocks; it is the same receiver re-shaped
from batch (whole file, one satellite at a time) into streaming (all satellites, sample by
sample), which is where the education is.

Scope: GPS L1 C/A only, one antenna, one radio. No Galileo (numpy-gps has it; later). No
assistance, no almanac, no RTK. Never anything encrypted.

## The receiver, as blocks

```
                     â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
 Soapy Source â”€â”€â”€â”€â”€â”€â–ºâ”‚ gpsrx Acquisition                                                  â”‚
 2.048 MS/s, cf32    â”‚  every N s: snapshot 4 ms, FFT search 32 PRNs x Doppler bins       â”‚â”€â”€â–º 'sky' (msg: PRN, Doppler, code phase, metric)
 bias-T on           â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
        â”‚
        â”‚  the same stream, to every channel
        â”œâ”€â”€â–º gpsrx Channel (PRN a) â”€â”€â–º 'prompt' (1 kHz I/Q), 'obs' (code phase + carrier + epoch counter, 1 Hz)
        â”œâ”€â”€â–º gpsrx Channel (PRN b) â”€â”€â–º ...
        â””â”€â”€â–º gpsrx Channel (PRN n)
                       â”‚ 'prompt'                      â”‚ 'obs'
                       â–¼                               â–¼
              gpsrx Nav Decoder (per channel)   gpsrx PVT Solver
              bits -> preamble -> parity ->     pseudoranges from >=4 channels at a common
              subframes -> ephemeris â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â–º receive epoch -> least squares -> position,
              'eph' (msg)                       clock bias, DOP â”€â”€â–º 'fix' (msg)
                                                                    â”‚
                                                          gpsrx Sky Panel (Qt)
                                                          sky plot, C/N0 bars, eye of the
                                                          nav bits, fix quality (not coordinates)
```

Six blocks, five in Python:

| block | in | out | numpy-gps source | notes |
|---|---|---|---|---|
| **Acquisition** | stream | `sky` msg per detection | `measure.acquire()` | Runs on a 4 ms snapshot every few seconds on its own thread (the ATSC 3.0 lesson: never do heavy work in `work()`). Emits only NEW satellites, and re-emits a lost one after the channel reports loss. |
| **Channel** | stream | `prompt` (msg or a 1 kHz float stream), `obs` msg | `measure.track_sv()`, `prompts_ms()` | One instance per satellite, created by the flowgraph from a channel *bank* (see below). NCO, C/A code with fractional phase, early/prompt/late correlators, DLL (early-late envelope) and PLL (Costas). Works on whole-millisecond chunks; carries an integer code-epoch counter from acquisition onward - **that counter is the pseudorange**. |
| **Nav Decoder** | `prompt` | `eph` msg | `fix.decode_eph()` | Bit sync (20 ms histogram), preamble + parity (IS-GPS-200 Â§20.3.5), subframes 1-3 -> ephemeris, 4/18 -> Klobuchar. Publishes the TOW of each subframe with the channel's epoch counter at that instant: the **timing anchor** numpy-gps calls it. |
| **PVT Solver** | `obs` x N, `eph` x N | `fix` msg | `fix.solve_snapshot()`, `sat_ecef`, `clock_corr`, `klobuchar`, `tropo_delay` | At each solve epoch: for every channel with an ephemeris, transmit time = anchored TOW + (epochs since anchor) x 1 ms - code phase (numpy-gps rule 1: phase is to the NEXT epoch), in SV time; SV clock correction after (rule 2); common receive epoch = max transmit time + ~70 ms, then the same relative-integer search numpy-gps uses (rule 3). Least squares with elevation weights, tropo + iono. |
| **Sky Panel** (Qt) | `sky`, `obs`, `fix` | - | `dash.py` | Az/el sky plot from the ephemerides, C/N0 per channel, bit eye, the fix as **count, rms, DOP, plausible altitude** - coordinates go to a file the user names, never to the panel, never to a screenshot. |
| **Channel Bank** | - | - | - | Not a block: a Python helper (and later a hier block) that owns N Channel instances and a Nav Decoder each, connects them, and maps acquisition messages onto free slots. GNU Radio cannot add blocks to a running flowgraph, so the bank pre-allocates (12) and idles the unused ones. This is what gnss-sdr's channel manager does. |

Ports carry what the textbook draws: the acquisition grid, the three correlator magnitudes,
the prompt I/Q constellation (two blobs = tracking), the nav bits, the ephemeris fields. A
Message Debug block on any port is a lesson.

## Where it will be hard, and the plan for each

**1. Tracking throughput in Python.** 8-12 channels x 3 correlators x 2048 samples/ms = 50-75 M
multiply-adds a second, plus code generation and NCO. Feasible in NumPy *if* each channel
handles whole milliseconds per `work()` call, generates its code and carrier for the chunk
vectorised (numpy-gps's `sampled_code` + a phase ramp), and does the three correlations as one
matrix product. annappo/GPS-SDR-Receiver tracks 12 birds in real time in Python this way. If
it does not keep up, the Channel is the one block that goes to C++ (a 200-line sync block; the
canonical gr_modtool exercise), built on the Ubuntu rig - this PC has no compiler. **Gate:**
8 channels at 1.0x real time on this PC before anything else is built on top.

**2. Time.** A GNU Radio stream has no clock. Each Channel counts code epochs from the sample
where acquisition handed it its first code phase; the Nav Decoder ties one epoch count to a
TOW; from then on every epoch is 1 ms of SV time. The PVT block needs every channel's epoch
count *at the same receive sample*, so channels tag their `obs` with the absolute sample index
(`nitems_read`) and the solver interpolates each channel's code phase to a common sample.
numpy-gps does this offline by construction (one file, one clock); here it is bookkeeping,
and it is exactly where the July bugs lived. **Gate:** pseudoranges from the blocks equal
numpy-gps's from the same capture to better than 1 m (a few ns) before a fix is attempted.

**3. Acquisition on a live stream.** 32 PRNs x 57 Doppler bins x 4 ms FFTs is ~2 s of CPU;
on its own thread with a snapshot, as gr-atsc3rx's Frame Sync learned to. Re-acquisition of
a lost channel is the same search restricted to one PRN.

**4. Static flowgraphs vs dynamic satellites.** The bank of 12 idle channels. A channel with
nothing assigned consumes and discards its input cheaply (returns immediately).

## Gate 1 result (2026-09-22): the Channel goes to C++

Measured, eight synthetic satellites, `gpsrx.channel` Python blocks in the GNU Radio scheduler,
one thread per block, on a desktop:

| channels | wall for 1 s of samples | real-time factor |
|---|---|---|
| 1 | 0.32 s | 3.1x |
| 2 | 0.89 s | 1.1x |
| 4 | 2.49 s | 0.40x |
| 8 | 5.90 s | 0.17x |

One channel is fast; eight are slower than eight times one. The per-period NumPy operations
(a 2048-point complex multiply, three correlations) are ~10 us each, so the interpreter lock
changes hands thousands of times a second between eight threads and the handoff costs more
than the arithmetic. Batching periods per `work()` call did not help (the upstream buffer is
8191 items, four periods, and is not settable from Python); caching the NCO ramp did not help
(the loop moves the Doppler every period while pulling in). The engine alone, single-threaded,
is 0.75x. **So: the tracking channel is the one C++ block**, built on the Ubuntu rig; the
Python engine in `track.py` stays as the readable reference and the block's test oracle
(gate 0 for the C++ block = its prompts equal the Python engine's on the synthetic sky).
Everything else - acquisition, nav decode, PVT, panel - is low-rate and stays Python.

## Gates (in order)

0. **Correlator identity:** the Channel's prompt output on a capture equals numpy-gps's
   `prompts_ms()` for the same PRN/Doppler/phase, bit for bit (both are NumPy; float64).
1. **Tracking throughput:** 8 channels at >= 1.0x real time, Python, this PC.
2. **Pseudorange identity:** blocks vs numpy-gps on the same capture, < 1 m.
3. **Fix identity:** the PVT block's fix from those pseudoranges vs numpy-gps's `--resolve`
   from its cache, same answer to the metre (same solver code, so it must).
4. **Live:** a fix from the attic antenna inside a GRC window, gr-rxtune holding the gain
   (the loop showed 9/21 the LNA state is the only knob that matters; the flowgraph should
   show that too).
5. **Educational:** a stranger can open `gpsrx_live.grc`, read the six blocks, and see each
   stage's output on a panel. That gate is a README walkthrough with screenshots.

## Repository

`gr-gpsrx`, GPL-3.0-or-later (GNU Radio convention). numpy-gps is MIT, so its algorithms can be
carried in with attribution; code flows from numpy-gps into here, never back. gr_modtool layout,
Python blocks under `python/gpsrx/`, GRC definitions, `examples/` with a replay flowgraph and
a live Qt flowgraph written from one generator (as the other modules), `util/gate_*.py` for
gates 0-3, QA against a synthetic signal (numpy-gps's `sampled_code` can *make* a fake
satellite: known PRN, Doppler, phase - so CI needs no capture and no radio).

Oracle captures (numpy-gps `lab_local/sky_capture_*.cs16`) never enter the repository; QA
uses synthesised signals; the README's screenshots show sky plots and C/N0, not coordinates.

## Order of work

1. Scaffold + Channel block + gate 0 (correlator identity on a synthetic satellite). One evening.
2. Gate 1 (throughput). Decides Python vs C++ for the Channel before anything depends on it.
3. Acquisition block + Channel Bank + a replay flowgraph that tracks every visible bird from a
   capture, with the Sky Panel showing C/N0. First thing worth a screenshot.
4. Nav Decoder + gate 2 (pseudoranges). The hard one.
5. PVT Solver + gate 3. Then gate 4, live, with gr-rxtune.
6. Walkthrough README.

What it is NOT: a replacement for numpy-gps (which stays the offline oracle and the war-drive
tool) or a competitor to gnss-sdr (which is the receiver to use). It is the one you can read.
