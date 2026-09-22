# SPDX-License-Identifier: GPL-3.0-or-later
'''
gr-gpsrx: a GPS L1 C/A receiver made of GNU Radio blocks. The one you can read.
'''
try:
    from .gpsrx_python import *          # the C++ blocks, when built (channel_cc)
except ImportError:                      # not built, or GNU Radio's libraries not loadable here
    pass

from . import acquire, cacode, nav, pvt, synth, track     # the engines: importable without GNU Radio

try:
    from .acquisition import acquisition
    from .channel import channel
    from .nav_decoder import nav_decoder
    from .pvt_solver import pvt_solver
    from .receiver import receiver
    from .status_sink import status_sink
except ImportError:                    # gnuradio not installed: the engines still work
    pass

try:
    from .sky_panel import sky_panel   # needs PyQt5; the rest of the module does not
except ImportError:
    pass
