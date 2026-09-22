# gr-gpsrx

**A GPS L1 C/A receiver made of GNU Radio blocks - the one you can read.**

> **Status: the tracking engine exists and is tested; there are no GNU Radio blocks yet.**
> `python/gpsrx/track.py` is one satellite's channel - NCO, C/A code, early/prompt/late
> correlators, PLL and DLL - one code period at a time, in NumPy, with nothing from GNU Radio
> in it so that it can be tested and timed under plain pytest. It holds the design's gate 0
> (its open-loop correlator equals [numpy-gps](https://github.com/Felbs/numpy-gps)'s to float
> precision) and locks eight synthetic satellites at once. `docs/DESIGN.md` has the whole
> receiver: six blocks, five gates, the order of work; `docs/PRIOR_ART.md` the field.

The point is educational: every stage of a GPS receiver as a block on the canvas, every
intermediate as a message you can plot, and a real position at the end. The reference and
oracle is numpy-gps, the offline pure-NumPy receiver; the receiver to *use* is
[gnss-sdr](https://github.com/gnss-sdr/gnss-sdr). This is the textbook, in the tool people use.

## What works today

```
pip install numpy pytest
NUMPY_GPS_DIR=/path/to/numpy-gps pytest -s tests      # 6 tests; gate 0 needs numpy-gps, the rest need nothing
```

`python/gpsrx/synth.py` makes a synthetic sky - any satellites, Dopplers, code phases, C/N0s,
data bits, one noise floor - so every test runs with no capture and no radio. Measured on the
synthetic sky: eight channels lock (PLL lock indicator > 0.94) and track at **0.75x real time**
in Python on a desktop, which is under the design's bar of 1.0x. That decision (Python vs a
C++ channel block) is next.

## Licence

GPL-3.0-or-later. The C/A code generator and the algorithms come from numpy-gps (MIT, Felbs).
No captures, no positions, no coordinates are in this repository, ever.
