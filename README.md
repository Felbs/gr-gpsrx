# gr-gpsrx

**A GPS L1 C/A receiver made of GNU Radio blocks - the one you can read.**

> **Status (2026-09-22): it fixes, live.** Every stage is a GNU Radio block (Acquisition,
> Channel in Python *and* C++, Nav Decoder, PVT, Status, Sky Panel, and a one-block Receiver);
> the four example flowgraphs validate and compile in GRC; the grcc-generated radio flowgraph
> run unedited on an SDRplay RSPdx gave a fix from **8 satellites, residual rms 2.6 m, PDOP
> 2.9, 15-fix scatter 9.8 m** within 200 s of start, 7 m from the fix the same receiver gets
> on a recording from a different constellation hours earlier. Eight C++ channels run at 15x
> real time. Built and run on Windows (radioconda + MSVC); Linux is the same CMake but has
> not been exercised yet. `docs/DESIGN.md` has the receiver; `docs/TEST_REPORT.md` the
> evidence and the defects found on the way.

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
| Channel | `track.py`, block in `channel.py`; C++ twin `lib/channel_cc_impl.cc` | one satellite: wipe-off, three correlators, Costas PLL, early-late DLL, one code period per step; **counts code epochs**. Same ports, tag and messages in both languages; the Python one is for reading, the C++ one for the radio |
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

util\build_win.cmd                    # Windows: the C++ Channel (radioconda + VS Build Tools); Linux: cmake as any OOT
python python/gpsrx/qa_channel_cc.py  # C++ vs Python on the same sky, and gate 1

python apps/gpsrx_replay.py capture.cs16 --fix-file /somewhere/private/fix.json          # a recording
grcc -o build/grc examples/*.grc
python util/run_qt_shot.py build/grc/gpsrx_radio_qt.py --set "antenna=Antenna B" \
       --call "src.write_setting('biasT_ctrl','true')" --set fix_file=/somewhere/private/live.json --seconds 200
gnuradio-companion examples/gpsrx_canvas.grc                      # the one to read
```

From the source tree, `examples/_devpath.py` grafts `python/` onto the `gnuradio` package and
picks up the built C++ module (the apps and QA do this themselves); for GRC set
`GRC_BLOCKS_PATH` to this repo's `grc/`. A capture is interleaved int16 I/Q at 2.048 MS/s
centred on 1575.42 MHz. The receiver prints satellites, lock, subframes, residual rms, PDOP and
scatter; **the position goes only where you point the fix file and is never printed.** Keep
that file out of any repository.

![the receiver on the canvas](docs/img/grc_canvas.png)

## What is not done

- Linux build and CI (the CMake is gr_modtool's; nothing platform-specific in the C++).
- Acquisition is snapshot-and-search every N seconds; a lost satellite waits for the next
  interval, and a satellite rising mid-run is only picked up then.
- Carrier-phase observables, a Kalman filter, other constellations: numpy-gps has some of
  these offline; they are not here.
- The walkthrough (each block's page with its plots) is the next document.

## Licence

GPL-3.0-or-later. The C/A generator, the orbit/clock/atmosphere math and the field parser come
from numpy-gps (MIT, Felbs). No captures, no positions, no coordinates are in this repository.
