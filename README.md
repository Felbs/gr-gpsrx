# gr-gpsrx

**A GPS L1 C/A receiver made of GNU Radio blocks - the one you can read.**

> **Status (2026-09-22): the receiver works, offline.** Acquisition, tracking channels, the
> navigation decoder and the position solver exist as pure-NumPy engines with 13 tests, and
> `apps/gpsrx_offline.py` runs the whole chain on a capture: two 2-4 minute recordings from an
> attic antenna gave fixes from 8 and 7 satellites with **residual rms 0.3-3.7 m** and
> **10 m repeatability** between them, two hours apart. Only the tracking Channel is a GNU
> Radio block so far, and it is going to C++ (below). `docs/DESIGN.md` has the whole receiver;
> `docs/PRIOR_ART.md` the field; `docs/TEST_REPORT.md` what was and was not run.

The point is educational: every stage of a GPS receiver as a block on the canvas, every
intermediate as a message you can plot, and a real position at the end. The reference and
oracle is [numpy-gps](https://github.com/Felbs/numpy-gps), the offline pure-NumPy receiver; the
receiver to *use* is [gnss-sdr](https://github.com/gnss-sdr/gnss-sdr). This is the textbook, in
the tool people use.

## The receiver

```
 samples ─► Acquisition ─► Channel x N (NCO, C/A code, early/prompt/late, PLL, DLL) ─► prompts (1 kHz)
                                │ observables: code-epoch count + the epoch's arrival sample        │
                                ▼                                                                   ▼
                          PVT Solver  ◄──── ephemeris + timing anchor (period <-> TOW) ◄──── Nav Decoder
```

| stage | file | what it is |
|---|---|---|
| C/A code | `cacode.py` | IS-GPS-200 generator, checked against the published first chips |
| Channel | `track.py`, block in `channel.py` | one satellite: wipe-off, three correlators, Costas PLL, early-late DLL, one code period per step; **counts code epochs** |
| Nav decoder | `nav.py` | bit sync, preamble/parity framing in either polarity, subframes 1-4 -> ephemeris and Klobuchar terms, timing anchors; rejects false frames by TOW range and subframe continuity |
| PVT | `pvt.py` | orbit, SV clock, Sagnac, troposphere, ionosphere, least squares (numpy-gps's math) |
| synthetic sky | `synth.py`, `navgen.py` | any satellites, Dopplers, C/N0s, real nav messages with parity - so every test runs with no capture and no radio |

**What the streaming design buys.** numpy-gps decodes offline and must extrapolate each
satellite's time from a fitted anchor, which leaves a whole-millisecond ambiguity per satellite
that it searches over. A channel that runs continuously *counts*: every code period after the
anchored subframe is exactly 1 ms of that satellite's clock, and the epoch's arrival is reported
as a fractional input sample. No ambiguity, no search - and on the same two captures this
receiver repeated to 9.9 m where numpy-gps repeated to 102 m.

## Run it

```
pip install numpy pytest
NUMPY_GPS_DIR=/path/to/numpy-gps pytest -q tests                 # 13 tests, ~25 s
NUMPY_GPS_DIR=/path/to/numpy-gps python apps/gpsrx_offline.py capture.cs16 --secs 120
```

The capture is interleaved int16 I/Q at 2.048 MS/s, centred on 1575.42 MHz (numpy-gps's
`locate.py` records one from any SoapySDR radio). Acquisition is borrowed from numpy-gps in
this phase. The console prints satellites, lock, subframes, residual rms, PDOP and whether the
altitude is plausible; **the position goes only to a file in a gitignored directory and is
never printed.**

## What is not done

- **Only the Channel is a GNU Radio block.** Eight Python channels in the scheduler run at 0.17x
  real time (the interpreter lock changes hands thousands of times a second): the Channel goes
  to C++ (design gate 1). Acquisition, Nav Decoder, PVT and the Sky Panel blocks are next, all
  low-rate Python.
- No live run yet; no GRC flowgraph yet; gate 2 (pseudoranges against numpy-gps to the metre) is
  superseded by the finding above and needs a truth position to settle.

## Licence

GPL-3.0-or-later. The C/A generator, the orbit/clock/atmosphere math and the field parser come
from numpy-gps (MIT, Felbs). No captures, no positions, no coordinates are in this repository.
