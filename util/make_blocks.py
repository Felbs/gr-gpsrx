#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes the GRC block definitions (grc/*.block.yml) from one description, so the canvas and
the Python blocks cannot drift apart. Re-run after changing a block's parameters."""
import os

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RATE = {"id": "samp_rate", "label": "Sample rate", "dtype": "real", "default": "samp_rate"}
SLOT = {"id": "slot", "label": "Slot", "dtype": "int", "default": "0"}

QT_MAKE = """<%
    win = '_%s_win' % id
%>    ${win} = gpsrx.BLOCK(label=${label}, parent=self, private=${private})
self.${id} = ${win}
${gui_hint() % win}"""

BLOCKS = [
    dict(id="gpsrx_acquisition", label="GPS Acquisition", doc="""\
Complex baseband in (L1 C/A at any rate; 2.048 MS/s is the tested one). Finds the satellites.
Every Interval seconds it copies a snapshot and searches it on its own thread: parallel
code-phase search, 1 ms coherent x N non-coherent, over +-7 kHz of Doppler in 250 Hz bins,
then refines each Doppler to ~10 Hz. Publishes what it found on 'sky' and assigns each new
satellite to a free Channel slot on 'assign' (connect it to every Channel's 'assign'; a
Channel only acts on its own slot). Connect every Channel's 'status' back to 'status' so
the slot table knows which are busy.
Threshold: peak over second peak; 2.5 with 100 ms separates real satellites (>3) from
noise (<1.8), measured on air. Hold: for file replay only - stops the stream while a
search runs, otherwise a file source races past an idle receiver.""",
         make="gpsrx.acquisition(samp_rate=${samp_rate}, n_slots=${n_slots}, snapshot_ms=${snapshot_ms}, "
              "interval_s=${interval_s}, threshold=${threshold}, n_noncoh=${n_noncoh}, hold=${hold}, settle_s=${settle_s})",
         params=[RATE, {"id": "n_slots", "label": "Channel slots", "dtype": "int", "default": "8"},
                 {"id": "interval_s", "label": "Interval (s)", "dtype": "real", "default": "20.0"},
                 {"id": "threshold", "label": "Threshold", "dtype": "real", "default": "2.5"},
                 {"id": "n_noncoh", "label": "Non-coherent ms", "dtype": "int", "default": "100", "hide": "part"},
                 {"id": "snapshot_ms", "label": "Snapshot (ms)", "dtype": "int", "default": "110", "hide": "part"},
                 {"id": "settle_s", "label": "First search after (s)", "dtype": "real", "default": "0.5", "hide": "part"},
                 {"id": "hold", "label": "Hold the stream while searching", "dtype": "enum", "default": "False",
                  "options": ["False", "True"], "option_labels": ["No (radio)", "Yes (file replay)"]}],
         inputs=[{"domain": "stream", "dtype": "complex"}, {"domain": "message", "id": "status", "optional": True}],
         outputs=[{"domain": "message", "id": "sky", "optional": True}, {"domain": "message", "id": "assign"}]),
    dict(id="gpsrx_channel", label="GPS Channel", doc="""\
One satellite's tracking loop. Complex baseband in; the prompt correlation out, one item per
code period (1 kHz complex: the sign of the real part is the 50 bit/s navigation data).
Carrier wipe-off (NCO), early/prompt/late C/A code correlators, Costas PLL (atan(Q/I),
blind to the data flips), early-late DLL, carrier-aided code rate. Consumes exactly one
code period per step, so its period count IS the satellite's clock; 'obs' reports once a
second the period count and the code epoch's arrival as a fractional input sample - the
two numbers a pseudorange is made of.
Two stages: wide loops and 1 ms integration to pull in, then - once the data-bit edges
are found - NARROW loops and coherent integration over a whole bit (as gnss-sdr does).
Static antenna: 5 Hz narrow PLL; moving: 8-15 Hz. 0 disables the second stage.
Idle until 'assign' names its slot (from the Acquisition block). 'status' says tracking /
lost / idle. Slot: this channel's number; Acquisition assigns by it.""",
         make="gpsrx.channel(samp_rate=${samp_rate}, slot=${slot}, pll_bw=${pll_bw}, dll_bw=${dll_bw}, "
              "obs_every_ms=${obs_every_ms}, pll_bw_narrow=${pll_bw_narrow}, dll_bw_narrow=${dll_bw_narrow}, coherent_ms=${coherent_ms})",
         params=[RATE, SLOT, {"id": "pll_bw", "label": "PLL bandwidth (Hz)", "dtype": "real", "default": "18.0"},
                 {"id": "dll_bw", "label": "DLL bandwidth (Hz)", "dtype": "real", "default": "2.0"},
                 {"id": "obs_every_ms", "label": "Observable every (ms)", "dtype": "int", "default": "1000", "hide": "part"},
                 {"id": "pll_bw_narrow", "label": "PLL bandwidth after bit sync (Hz)", "dtype": "real", "default": "15.0"},
                 {"id": "dll_bw_narrow", "label": "DLL bandwidth after bit sync (Hz)", "dtype": "real", "default": "0.5"},
                 {"id": "coherent_ms", "label": "Coherent integration after bit sync (ms)", "dtype": "int", "default": "20"}],
         inputs=[{"domain": "stream", "dtype": "complex"}, {"domain": "message", "id": "assign"}],
         outputs=[{"domain": "stream", "dtype": "complex"}, {"domain": "message", "id": "obs"},
                  {"domain": "message", "id": "status", "optional": True}]),
    dict(id="gpsrx_channel_cc", label="GPS Channel (C++)", doc="""The same Channel as GPS Channel, in C++: eight of them run at 15x real time (measured; the
Python ones at 0.17x). Same ports, same tag, same messages, the same arithmetic line for
line - a Nav Decoder or PVT cannot tell which it is wired to. Use this one live.""",
         make="gpsrx.channel_cc(${samp_rate}, ${slot}, ${pll_bw}, ${dll_bw}, ${obs_every_ms}, ${pll_bw_narrow}, ${dll_bw_narrow}, ${coherent_ms})",
         params=[RATE, SLOT, {"id": "pll_bw", "label": "PLL bandwidth (Hz)", "dtype": "real", "default": "18.0"},
                 {"id": "dll_bw", "label": "DLL bandwidth (Hz)", "dtype": "real", "default": "2.0"},
                 {"id": "obs_every_ms", "label": "Observable every (ms)", "dtype": "int", "default": "1000", "hide": "part"},
                 {"id": "pll_bw_narrow", "label": "PLL bandwidth after bit sync (Hz)", "dtype": "real", "default": "15.0"},
                 {"id": "dll_bw_narrow", "label": "DLL bandwidth after bit sync (Hz)", "dtype": "real", "default": "0.5"},
                 {"id": "coherent_ms", "label": "Coherent integration after bit sync (ms)", "dtype": "int", "default": "20"}],
         inputs=[{"domain": "stream", "dtype": "complex"}, {"domain": "message", "id": "assign"}],
         outputs=[{"domain": "stream", "dtype": "complex"}, {"domain": "message", "id": "obs"},
                  {"domain": "message", "id": "status", "optional": True}]),
    dict(id="gpsrx_nav_decoder", label="GPS Nav Decoder", doc="""\
A Channel's prompt stream in (1 kHz); the navigation message out. Bit sync (where the 20
code periods of each bit begin), words, parity in either polarity, subframes 1-5,
ephemeris, Klobuchar terms from subframe 4 page 18 - and the timing anchor: which code
period a subframe began on and the TOW it carries. Publishes on 'nav' after every accepted
subframe. Slot: must match the Channel it is wired to.""",
         make="gpsrx.nav_decoder(slot=${slot})", params=[SLOT],
         inputs=[{"domain": "stream", "dtype": "complex"}], outputs=[{"domain": "message", "id": "nav"}]),
    dict(id="gpsrx_pvt", label="GPS PVT", doc="""\
Every Channel's 'obs' and every Nav Decoder's 'nav' in; the position out on 'fix', at most
once a second, when four or more channels have an observable, an ephemeris and a timing
anchor. Transmit time = anchor TOW + period count (no millisecond ambiguity); all channels
referred to one receive sample; SV clock, Sagnac, troposphere and ionosphere; least squares.
The 'fix' message carries n, prns, rms_m, pdop, altitude_plausible, ecef, llh and a running
mean/scatter over the last Average fixes. Fix file: optional path to keep the latest fix
(keep it out of any repository). Iono file: an archived broadcast to use until a channel
decodes page 18. Connect every Channel's 'status' too: a lost channel's last
observable is dropped, and any observable more than 2.5 s older than the newest is ignored.""",
         make="gpsrx.pvt_solver(samp_rate=${samp_rate}, average=${average}, iono_file=${iono_file}, fix_file=${fix_file}, eph_file=${eph_file})",
         params=[RATE, {"id": "average", "label": "Average fixes", "dtype": "int", "default": "15"},
                 {"id": "iono_file", "label": "Iono file", "dtype": "file_open", "default": "''"},
                 {"id": "fix_file", "label": "Fix file", "dtype": "file_save", "default": "''"},
                 {"id": "eph_file", "label": "Ephemeris file (warm start)", "dtype": "file_save", "default": "''"}],
         inputs=[{"domain": "message", "id": "obs"}, {"domain": "message", "id": "nav"},
                 {"domain": "message", "id": "status", "optional": True}],
         outputs=[{"domain": "message", "id": "fix"}]),
    dict(id="gpsrx_status", label="GPS Status", doc="""\
The receiver's state as text, one table every Every seconds: satellites per slot, lock,
Doppler, seconds tracked, subframes, and the fix's quality (satellites, residual rms, PDOP,
scatter). It never prints the position. Connect whichever of 'sky', 'status', 'obs', 'nav',
'fix' you have. Log file: append the tables there.""",
         make="gpsrx.status_sink(every_s=${every_s}, log_file=${log_file})",
         params=[{"id": "every_s", "label": "Every (s)", "dtype": "real", "default": "2.0"},
                 {"id": "log_file", "label": "Log file", "dtype": "file_save", "default": "''"}],
         inputs=[{"domain": "message", "id": p, "optional": True} for p in ("sky", "status", "obs", "nav", "fix")],
         outputs=[]),
    dict(id="gpsrx_receiver", label="GPS Receiver", doc="""\
The whole receiver in one block: Acquisition -> N x (Channel -> Nav Decoder) -> PVT.
Complex baseband in; 'fix' out, plus every inner message port ('sky', 'status', 'obs',
'nav') for a Status block or a Message Debug. The same wiring laid out block by block is
examples/gpsrx_canvas.grc. Hold: file replay only. Channel engine: the C++
Channel keeps up with a radio; the Python one is for replay and reading.""",
         make="gpsrx.receiver(samp_rate=${samp_rate}, n_channels=${n_channels}, interval_s=${interval_s}, "
              "threshold=${threshold}, iono_file=${iono_file}, fix_file=${fix_file}, hold=${hold}, engine=${engine}, eph_file=${eph_file})",
         params=[RATE, {"id": "n_channels", "label": "Channels", "dtype": "int", "default": "8"},
                 {"id": "engine", "label": "Channel engine", "dtype": "enum", "default": "'auto'",
                  "options": ["'auto'", "'cpp'", "'python'"], "option_labels": ["C++ if built", "C++ (live)", "Python (replay only)"]},
                 {"id": "interval_s", "label": "Search interval (s)", "dtype": "real", "default": "20.0"},
                 {"id": "threshold", "label": "Threshold", "dtype": "real", "default": "2.5"},
                 {"id": "iono_file", "label": "Iono file", "dtype": "file_open", "default": "''"},
                 {"id": "fix_file", "label": "Fix file", "dtype": "file_save", "default": "''"},
                 {"id": "hold", "label": "Hold the stream while searching", "dtype": "enum", "default": "False",
                  "options": ["False", "True"], "option_labels": ["No (radio)", "Yes (file replay)"]},
                 {"id": "eph_file", "label": "Ephemeris file (warm start)", "dtype": "file_save", "default": "''"}],
         inputs=[{"domain": "stream", "dtype": "complex"}],
         outputs=[{"domain": "message", "id": p, "optional": True} for p in ("sky", "status", "obs", "nav", "fix")]),
    dict(id="gpsrx_sky_panel", label="GPS Sky Panel", doc="""Qt panel: the sky as a polar plot with every tracked satellite where the fix puts it
(coloured by lock; on a Doppler ring before the first fix), one line per channel, and the
fix's quality - satellites, residual rms, PDOP, scatter. It never shows the position.
Private view: for screenshots - PRN numbers hidden and the sky turned by an undisclosed
angle, because a named constellation at a known time can be inverted to a rough position.
Optional: the receiver runs without it. Docks like any QT GUI widget (GUI Hint).""",
         make=None,
         params=[{"id": "label", "label": "Label", "dtype": "string", "default": "GPS receiver"},
                 {"id": "private", "label": "Private view (screenshots: no PRNs, sky turned)", "dtype": "bool", "default": "False"},
                 {"id": "gui_hint", "label": "GUI Hint", "dtype": "gui_hint", "hide": "part"}],
         inputs=[{"domain": "message", "id": p, "optional": True} for p in ("sky", "status", "obs", "nav", "fix")],
         outputs=[]),
]


def main():
    names = []
    for b in BLOCKS:
        doc = {"id": b["id"], "label": b["label"], "category": "[gpsrx]", "flags": ["python"],
               "documentation": b["doc"],
               "templates": {"imports": "from gnuradio import gpsrx", "make": b["make"] or QT_MAKE.replace("BLOCK", b["id"][6:])},
               "parameters": b["params"], "inputs": b["inputs"]}
        if b["outputs"]:
            doc["outputs"] = b["outputs"]
        doc["file_format"] = 1
        path = os.path.join(ROOT, "grc", b["id"] + ".block.yml")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=None, width=100)
        names.append(b["id"] + ".block.yml")
        print("wrote", path)
    cm = os.path.join(ROOT, "grc", "CMakeLists.txt")
    txt = open(cm, encoding="utf-8").read()
    head = txt[:txt.index("install(FILES")]
    with open(cm, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(head + "install(FILES\n" + "".join(f"    {n}\n" for n in names) + "    DESTINATION share/gnuradio/grc/blocks)\n")


if __name__ == "__main__":
    main()
