# gr-gpsrx handoff — 2026-09-22 (tracking engine built + tested; LOCAL git only, no GitHub repo)

A GPS L1 C/A receiver made of GNU Radio blocks, for the education of it, held against
numpy-gps (Z:\src\numpy-gps, the offline pure-NumPy receiver; formerly GPSTuna) as its oracle.
Read docs/PRIOR_ART.md (name check: gr-gpsrx is free; gr-gps is taken twice) then docs/DESIGN.md
(six blocks, five gates, order of work). Owner's direction 9/21: "build a GPS out of GNU Radio
blocks — it will be educational"; a numpy-gps *bridge* was proposed and declined.

Rules carried over from the sibling projects: no repo creation or push without the owner's
word; dox gate before every push; NO coordinates, captures or fix results in the tree
(numpy-gps keeps fixes in gitignored lab_local — same here); GPU optional; radio via
RXTUNE_LOCK. First gate to build: correlator identity on a SYNTHETIC satellite (numpy-gps's
sampled_code makes one), so QA never needs a capture.

## 9/22 — what exists
`python/gpsrx/{cacode,track,synth}.py` + `tests/test_track.py` (6 pass). gate 0 PASSES (open-loop prompts ==
numpy-gps prompts_ms to 1e-9). Closed loop: SoftGNSS filter form (tau1/tau2), PLL 18 Hz / DLL 2 Hz.
★ BUG FOUND THE HARD WAY: the Costas discriminator must be `atan(Q/I)` (single-argument), NOT `atan2(Q, I)` -
atan2 sees a 180-degree data flip as a 165-degree phase error and slews the carrier 100 Hz, so every bit
transition produced a 3-period glitch and the sign flipped back. Cost 2 hours. Also: sum synthetic satellites
over ONE noise floor (`synth.sky`), never sum single-bird signals (that sums their noise floors).
Throughput: 8 locked channels 0.75x real time (exp() is 109 us of a 170 us period; the NCO ramp is now cached
and rebuilt only when Doppler moves > 0.5 Hz, which did not help because the dots/multiplies are the rest).
Next: decide Python vs C++ Channel (gate 1), then Acquisition block + Channel Bank + replay flowgraph.
- `python/gpsrx/channel.py` = the Channel as a Python block (assign/obs/status msg ports, lost detection,
  absolute-sample stamping) + `qa_channel.py` (3 tests pass). ★ GATE 1 MEASURED: 1 ch 3.1x, 8 ch 0.17x -
  Python threads contend on the GIL; batching and ramp caching did not help. DECISION: Channel -> C++ on
  the Ubuntu rig (this PC has no compiler). The Python block stays as the reference for the C++ one.
  GR gateway traps: `set_min_input_buffer` is not exposed to Python; `set_output_multiple(n)` means a
  short batch is an ERROR (pad or return 0); a Python source must be kept referenced (`tb.keep`).

## 9/22 night — IT FIXES. Offline, both attic captures. LOCAL git @de48f21 + docs; NO GitHub repo.
`apps/gpsrx_offline.py CAPTURE --secs N` (needs NUMPY_GPS_DIR for acquisition) -> lab_local/fix.json (GITIGNORED;
coordinates never printed). 120 s: 8 birds rms 3.6 m; 240 s: 7 birds rms 0.3 m; repeat 9.9 m (numpy-gps 102 m on the
same two). The two receivers disagree by ~150 m (E+90 N-44 U+112) on the same capture; solvers agree (their prs
through ours = 9 m), so it is the pseudoranges - ours repeat, theirs wander in altitude. WHICH IS TRUE needs the
owner's known coordinates (compare lab_local/fix_240.json llh_mean vs numpy-gps lab_local/fix_result.json).
Fixed tonight: epoch arrival as fractional sample (the big one), false-frame rejection, atmosphere guard, Doppler
refinement before tracking. Next: Acquisition + Nav + PVT + Sky Panel as blocks (Python, low-rate), C++ Channel on
Ubuntu, replay GRC flowgraph, then live.

## 9/22 later — THE FLOWGRAPH FIXES. LOCAL @464ff65+. Blocks: acquisition, channel, nav_decoder, pvt_solver, status_sink,
sky_panel (Qt), receiver (hier). grc/*.block.yml <- util/make_blocks.py; examples/*.grc <- util/make_grc.py (replay,
canvas[6 ch], replay_qt, radio_qt); check: util/grc_screenshot.py (needs radioconda Library/bin on PATH), grcc -o build/grc.
Run: apps/gpsrx_replay.py CAPTURE --fix-file lab_local/x.json  |  util/run_qt_shot.py build/grc/gpsrx_replay_qt.py --set capture=..
--set fix_file=.. --seconds 330 --png lab_local/shot.png. Results: 7 birds rms 2.2-2.5 m PDOP 2.4 scatter 11-14 m, 9.9 m from
the offline fix. Private: lab_local/{fix_replay,fix_canvas,fix_qt}.json, sky_panel_shot.png (a sky plot at a known time is a
coarse location hint - user decides if it is ever shown). Next: C++ Channel on Ubuntu (gate 1), live run (radio_qt.grc, RSPdx
ANT-B bias-T), then the walkthrough.

## 9/22 dawn — LIVE FIX + C++ CHANNEL. LOCAL @a50487b+ (still no GitHub). This PC CAN build C++: VS 2022 Build Tools
(cl 14.44) + VS's cmake/ninja + radioconda headers -> util/build_win.cmd -> build/cpp/{lib/gnuradio-gpsrx.dll,
python/gpsrx/bindings/gpsrx_python.cp312-win_amd64.pyd}; examples/_devpath.py copies the pyd beside the package (gitignored)
and add_dll_directory's the lib. Binding is hand-kept (GEN_AUTOMATIC 0, HASH 0: castxml has no -fPIC on MSVC).
Gate 1 C++: 15.7x. Live: util/run_qt_shot.py build/grc/gpsrx_radio_qt.py --set "antenna=Antenna B" --call
"src.write_setting('biasT_ctrl','true')" --set fix_file=lab_local/fix_live.json --seconds 200 --png lab_local/live_shot.png
-> 8 sats rms 2.6 m PDOP 2.9; 7.2 m from the offline fix (different constellation). Private: lab_local/fix_live.json,
live_shot.png. Ubuntu/Pi were unreachable tonight (Linux build untested). Next: walkthrough doc, Linux build, GitHub when told.
