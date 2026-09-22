# gr-gpsrx

**A GPS L1 C/A receiver made of GNU Radio blocks - the one you can read.**

> **Status (2026-09-22): it fixes, as a flowgraph.** Every stage is a GNU Radio block
> (Acquisition, Channel, Nav Decoder, PVT, Status, Sky Panel, plus a one-block Receiver), the
> four example flowgraphs validate and compile in GRC, and the grcc-generated flowgraph run
> unedited on a 4-minute attic capture gives a fix from 7 satellites with **residual rms
> 2.2-2.5 m, PDOP 2.4, 15-fix scatter 11-14 m**, 10 m from the offline engine's fix on the same
> file. The tracking Channel is Python and runs eight satellites at 0.17x real time; it goes
> to C++ next. No live run yet. `docs/DESIGN.md` has the receiver; `docs/TEST_REPORT.md` the
> evidence and the defects.

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
| Acquisition | `acquire.py`, block in `acquisition.py` | parallel code-phase search, 1 ms coherent x N non-coherent, Doppler refined on 100 ms; assigns satellites to channel slots |
| Channel | `track.py`, block in `channel.py` | one satellite: wipe-off, three correlators, Costas PLL, early-late DLL, one code period per step; **counts code epochs** |
| Nav decoder | `nav.py`, block in `nav_decoder.py` | bit sync, preamble/parity framing in either polarity, subframes 1-4 -> ephemeris and Klobuchar terms, timing anchors; rejects false frames by TOW range and subframe continuity |
| PVT | `pvt.py`, block in `pvt_solver.py` | orbit, SV clock, Sagnac, troposphere, ionosphere, least squares (numpy-gps's math) |
| Status, Sky Panel | `status_sink.py`, `sky_panel.py` (Qt) | the receiver's state, quality only - never the position |
| Receiver | `receiver.py` | the whole chain as one hier block, every inner message port exposed |
| synthetic sky | `synth.py`, `navgen.py` | any satellites, Dopplers, C/N0s, real nav messages with parity - so every test runs with no capture and no radio |

**What the streaming design buys.** numpy-gps decodes offline and must extrapolate each
satellite's time from a fitted anchor, which leaves a whole-millisecond ambiguity per satellite
that it searches over. A channel that runs continuously *counts*: every code period after the
anchored subframe is exactly 1 ms of that satellite's clock, and the epoch's arrival is reported
as a fractional input sample. No ambiguity, no search - and on the same two captures this
receiver repeated to 9.9 m where numpy-gps repeated to 102 m.

## Run it

```
pip install numpy pytest pyyaml
NUMPY_GPS_DIR=/path/to/numpy-gps pytest -q tests                  # 14 engine tests, ~30 s, no GNU Radio
python python/gpsrx/qa_channel.py ; python python/gpsrx/qa_receiver.py    # the blocks, under GNU Radio

python apps/gpsrx_replay.py capture.cs16 --fix-file /somewhere/private/fix.json     # the flowgraph, on a file
grcc -o build/grc examples/gpsrx_canvas.grc && python build/grc/gpsrx_canvas.py --capture capture.cs16 --fix-file ...
gnuradio-companion examples/gpsrx_replay_qt.grc                   # with the Sky Panel
```

From the source tree, `examples/_devpath.py` grafts `python/` onto the `gnuradio` package (the
apps and QA do this themselves); for GRC set `GRC_BLOCKS_PATH` to this repo's `grc/`. The
capture is interleaved int16 I/Q at 2.048 MS/s centred on 1575.42 MHz (numpy-gps's `locate.py`
records one from any SoapySDR radio, or use `examples/gpsrx_radio_qt.grc` live). The receiver
prints satellites, lock, subframes, residual rms, PDOP and scatter; **the position goes only
where you point `--fix-file` and is never printed.** Keep that file out of any repository.

![the receiver on the canvas](docs/img/grc_canvas.png)

## What is not done

- **The Channel is Python** and eight of them run at 0.17x real time in the scheduler (the
  interpreter lock, not the arithmetic - design gate 1). It goes to C++; every other block is
  low-rate and stays Python. Until then: replay only, or fewer satellites live.
- No live run yet (`examples/gpsrx_radio_qt.grc` is written and compiles, not yet run on a
  radio). Acquisition is snapshot-and-search; re-acquiring a lost satellite waits for the next
  interval.
- The two receivers on the same captures disagree by ~150 m and only an independent
  position can say which is right (`docs/TEST_REPORT.md`).

## Licence

GPL-3.0-or-later. The C/A generator, the orbit/clock/atmosphere math and the field parser come
from numpy-gps (MIT, Felbs). No captures, no positions, no coordinates are in this repository.
