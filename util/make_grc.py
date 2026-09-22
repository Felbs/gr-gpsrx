#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""Writes the example flowgraphs from one description, so they cannot drift apart.

  examples/gpsrx_replay.grc   a capture file through the Receiver block, status as text
  examples/gpsrx_canvas.grc   the same receiver laid out block by block: Acquisition, six
                              Channels, six Nav Decoders, PVT, Status - the educational one
  examples/gpsrx_replay_qt.grc  the Receiver block with the Sky Panel, on a capture
  examples/gpsrx_radio_qt.grc   a SoapySDR radio through the Receiver, Sky Panel + spectrum

Nothing about any place is stored in either: the capture path and the fix file are
flowgraph Parameters. Re-run after editing, then check both with util/grc_check.py."""
import os

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRC_VERSION = "3.10.12.0"


def blk(name_, id_, x, y, **params):
    return {"name": name_, "id": id_, "parameters": {k: str(v) for k, v in params.items()},
            "states": {"bus_sink": False, "bus_source": False, "bus_structure": None,
                       "coordinate": [x, y], "rotation": 0, "state": "enabled"}}


def options(name, title, desc, qt=False):
    return {"parameters": {"id": name, "title": title, "author": "gr-gpsrx",
                           "generate_options": "qt_gui" if qt else "no_gui",
                           "output_language": "python", "category": "[GRC Hier Blocks]", "run": "True",
                           "run_options": "prompt" if qt else "run", "gen_cmake": "On", "description": desc},
            "states": {"bus_sink": False, "bus_source": False, "bus_structure": None,
                       "coordinate": [8, 8], "rotation": 0, "state": "enabled"}}


def source(y):
    return [
        blk("samp_rate", "variable", 230, 12, value="2.048e6", comment="2.048 MS/s: 2 samples per chip, the tested rate"),
        blk("capture", "parameter", 400, 12, label="Capture file (.cs16, interleaved int16 I/Q)", type="str",
            value='"capture.cs16"', short_id="c"),
        blk("fix_file", "parameter", 760, 12, label="Fix file (keep it out of any repository)", type="str",
            value='""', short_id="o"),
        blk("src", "blocks_file_source", 8, y, file="capture", type="short", repeat="False", offset=0, length=0,
            begin_tag="pmt.PMT_NIL", comment="int16 I,Q pairs"),
        blk("to_c", "blocks_interleaved_short_to_complex", 300, y + 12, vector_input="False", swap="False",
            scale_factor="32768.0"),
    ]


def build_replay():
    y = 220
    b = source(y) + [
        blk("rx", "gpsrx_receiver", 560, y - 40, samp_rate="samp_rate", n_channels=8, interval_s=100.0, threshold=2.5,
            iono_file='""', fix_file="fix_file", hold="True", engine="'auto'",
            comment="Acquisition -> 8 x (Channel -> Nav Decoder) -> PVT"),
        blk("status", "gpsrx_status", 900, y - 40, every_s=2.0, log_file='""',
            comment="quality only; the position goes to the fix file"),
    ]
    c = [["src", "0", "to_c", "0"], ["to_c", "0", "rx", "0"]] + [["rx", p, "status", p] for p in
                                                                   ("sky", "status", "obs", "nav", "fix")]
    return {"options": options("gpsrx_replay", "GPS receiver: capture file to a fix",
                               "File -> Receiver (one block) -> Status."),
            "blocks": b, "connections": c, "metadata": {"file_format": 1, "grc_version": GRC_VERSION}}, "gpsrx_replay"


def build_canvas(n=6):
    y = 300
    b = source(y) + [
        blk("acq", "gpsrx_acquisition", 560, y - 200, samp_rate="samp_rate", n_slots=n, interval_s=100.0,
            threshold=2.5, n_noncoh=100, snapshot_ms=110, settle_s=0.5, hold="True",
            comment="finds the satellites; assigns each to a Channel slot"),
        blk("pvt", "gpsrx_pvt", 1300, y + 60 * n, samp_rate="samp_rate", average=15, iono_file='""',
            fix_file="fix_file", comment="period counts + anchors -> pseudoranges -> least squares"),
        blk("status", "gpsrx_status", 1560, y + 60 * n, every_s=2.0, log_file='""'),
    ]
    c = [["src", "0", "to_c", "0"], ["to_c", "0", "acq", "0"], ["acq", "sky", "status", "sky"],
         ["pvt", "fix", "status", "fix"]]
    for s in range(n):
        yy = y + 120 * s
        b.append(blk(f"ch{s}", "gpsrx_channel_cc" if s % 2 else "gpsrx_channel", 720, yy, samp_rate="samp_rate", slot=s,
                     pll_bw=18.0, dll_bw=2.0, obs_every_ms=1000,
                     comment=f"slot {s}: NCO, E/P/L, PLL, DLL; counts code epochs" + (" (C++)" if s % 2 else " (Python)")))
        b.append(blk(f"nav{s}", "gpsrx_nav_decoder", 1020, yy, slot=s, comment="bits -> words -> subframes -> anchor"))
        c += [["to_c", "0", f"ch{s}", "0"], ["acq", "assign", f"ch{s}", "assign"], [f"ch{s}", "status", "acq", "status"],
              [f"ch{s}", "0", f"nav{s}", "0"], [f"ch{s}", "obs", "pvt", "obs"], [f"nav{s}", "nav", "pvt", "nav"],
              [f"ch{s}", "status", "status", "status"], [f"ch{s}", "obs", "status", "obs"], [f"nav{s}", "nav", "status", "nav"]]
    return {"options": options("gpsrx_canvas", "GPS receiver, block by block",
                               "Acquisition -> Channels -> Nav Decoders -> PVT -> Status, every wire visible."),
            "blocks": b, "connections": c, "metadata": {"file_format": 1, "grc_version": GRC_VERSION}}, "gpsrx_canvas"


def build_replay_qt():
    y = 220
    b = source(y) + [
        blk("rx", "gpsrx_receiver", 560, y - 40, samp_rate="samp_rate", n_channels=8, interval_s=100.0, threshold=2.5,
            iono_file='""', fix_file="fix_file", hold="True", engine="'auto'"),
        blk("panel", "gpsrx_sky_panel", 900, y - 40, label='"GPS receiver (replay)"', gui_hint="0,0,1,1"),
    ]
    c = [["src", "0", "to_c", "0"], ["to_c", "0", "rx", "0"]] + [["rx", p, "panel", p] for p in
                                                                   ("sky", "status", "obs", "nav", "fix")]
    return {"options": options("gpsrx_replay_qt", "GPS receiver: capture file, with the sky panel",
                               "File -> Receiver -> Sky Panel.", qt=True),
            "blocks": b, "connections": c, "metadata": {"file_format": 1, "grc_version": GRC_VERSION}}, "gpsrx_replay_qt"


def build_radio_qt():
    y = 300
    b = [
        blk("samp_rate", "variable", 230, 12, value="2.048e6"),
        blk("driver", "parameter", 400, 12, label="SoapySDR driver", type="str", value='"sdrplay"', short_id="d"),
        blk("antenna", "parameter", 600, 12, label="Antenna port", type="str", value='""', short_id="a"),
        blk("gain", "parameter", 780, 12, label="Overall gain (dB)", type="eng_float", value="40", short_id="g",
            comment="or tune the elements with gr-rxtune"),
        blk("settings", "parameter", 960, 12, label="Channel settings (key=value,...)", type="str",
            value='""', short_id="s",
            comment="bias-T for an active antenna is a DEVICE setting on SDRplay (biasT_ctrl), not a channel one: "
                    "util/run_qt_shot.py ... --call src.write_setting('biasT_ctrl','true')"),
        blk("fix_file", "parameter", 1240, 12, label="Fix file (keep it out of any repository)", type="str",
            value='""', short_id="o"),
        blk("src", "soapy_custom_source", 8, y - 20, driver="driver", type="fc32", nchan=1, dev_args='""',
            samp_rate="samp_rate", center_freq0="1575.42e6", bandwidth0="2.5e6", antenna0="antenna", gain0="gain",
            agc0=True, settings0="settings", minoutbuf=str(1 << 22),
            comment="L1 C/A: 1575.42 MHz. GPS is under the noise, so the AGC only sees noise - and that is fine"),
        blk("rx", "gpsrx_receiver", 560, y - 40, samp_rate="samp_rate", n_channels=8, interval_s=20.0, threshold=2.5,
            iono_file='""', fix_file="fix_file", hold="False", engine="'cpp'",
            comment="live needs the C++ Channel (util/build_win.cmd, or cmake on Linux)"),
        blk("spectrum", "qtgui_freq_sink_x", 330, y - 260, type="complex", name='"L1 baseband (the signal is below the noise)"',
            fftsize=1024, fc=0, bw="samp_rate", average=0.02, gui_hint="0,0,1,1"),
        blk("panel", "gpsrx_sky_panel", 900, y - 40, label='"GPS receiver"', gui_hint="1,0,1,1"),
    ]
    c = [["src", "0", "rx", "0"], ["src", "0", "spectrum", "0"]] + [["rx", p, "panel", p] for p in
                                                                    ("sky", "status", "obs", "nav", "fix")]
    return {"options": options("gpsrx_radio_qt", "GPS receiver: radio, with the sky panel",
                               "SoapySDR radio -> Receiver -> Sky Panel, spectrum.", qt=True),
            "blocks": b, "connections": c, "metadata": {"file_format": 1, "grc_version": GRC_VERSION}}, "gpsrx_radio_qt"


def main():
    for build in (build_replay, build_canvas, build_replay_qt, build_radio_qt):
        doc, name = build()
        path = os.path.join(ROOT, "examples", name + ".grc")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
