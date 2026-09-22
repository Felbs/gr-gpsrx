# SPDX-License-Identifier: GPL-3.0-or-later
'''
gr-gpsrx: a GPS L1 C/A receiver made of GNU Radio blocks. The one you can read.
'''
try:
    from .gpsrx_python import *          # no compiled bindings yet (Python phase)
except ModuleNotFoundError:
    pass

from . import cacode, synth, track     # the engine: importable without GNU Radio

try:
    from .channel import channel
except ImportError:                    # gnuradio not installed: the engine still works
    pass
