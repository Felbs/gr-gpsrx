# Prior art: GPS receivers in and around GNU Radio (checked 2026-09-21)

## Name check

| candidate | GitHub (exact name) | PyPI | CGRAN | verdict |
|---|---|---|---|---|
| `gr-gps` | **taken twice**: slaaja/gr-gps (2017, 3 blocks, abandoned), hahnpv/gr-gps (2019, "not a complete receiver, not in active development") | free | - | no: two dead repos own the name |
| `gr-gnss` | gnss-sdr's own satellite repos use the prefix (gr-gn3s, gr-dbfcttc) | free | - | no: reads as part of gnss-sdr |
| **`gr-gpsrx`** | free | free | free | **yes**: says receiver, matches gr-atsc3rx |
| `gr-numpygps` | free | free | free | ties the name to an implementation detail |
| `gr-gnssrx`, `gr-gps-rx`, `gr-l1ca` | free | free | free | fine, less clear |

CGRAN lists nothing GPS/GNSS-shaped at all (6 incidental mentions of "gps" on the index page).

## What exists

**gnss-sdr** (gnss-sdr/gnss-sdr, C++, GPL-3, 2278 stars, active the day of this check). The
real thing: a full multi-constellation receiver whose signal path IS a GNU Radio flowgraph
(signal source -> conditioner -> acquisition -> N channels of tracking -> telemetry decoder ->
observables -> PVT). Its architecture is the correct one and the block split below borrows it
deliberately. It is ~200k lines and a serious build; it is the receiver to *use*, not to learn
from by reading. Its tracking is `dll_pll_veml_tracking` (very-early/early/prompt/late/very-late
correlators, DLL + PLL, optional Kalman), one block instance per channel, channels created by a
channel manager on acquisition events. **This is the pattern to copy: one tracking block per
satellite, spawned by acquisition messages.**

**slaaja/gr-gps** (2017): three blocks - code generator, despreader, navdata - in the old XML
GRC era. Never reached a fix. **hahnpv/gr-gps** (2019): acquisition + a tracking block, the
author's own README sends readers to gnss-sdr and SoftGNSS. Both confirm the shape and that
nobody finished a small, readable one in GNU Radio.

**Python receivers, not GNU Radio.** perrysou/SoftGNSS-python (the classic textbook receiver
ported from MATLAB; offline). justkow/python-gps-receiver (offline, educational).
**annappo/GPS-SDR-Receiver** (16 stars): *real-time* tracking of up to 12 satellites in pure
Python from an RTL-SDR, a fix every 32 ms, 1-5 m standard deviation - the existence proof that
Python tracking keeps up on a PC when the loop is vectorised by the millisecond. And
**numpy-gps** (ours): offline, pure NumPy, GPS + Galileo, with the pseudorange-assembly rules
that took July to get right and a war-drive record. numpy-gps is this project's oracle.

## What is missing, that this project would be

A GPS L1 C/A receiver **made of GNU Radio blocks**, each stage a block you can watch, small
enough to read in an afternoon, in Python where Python keeps up, held bit-for-bit (pseudoranges,
then the fix) against numpy-gps on the same capture, and able to run live on a radio inside a
GRC window with the tuning loop (gr-rxtune) in charge of the gain. Nothing on the list above is
that. gnss-sdr is the goal state of the field; this is the textbook, in the tool people use.
