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
