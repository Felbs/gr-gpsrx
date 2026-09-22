# SPDX-License-Identifier: GPL-3.0-or-later
"""Run from the source tree with nothing installed: grafts python/ onto the gnuradio package."""
import os
import gnuradio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_py = os.path.join(ROOT, "python")
if _py not in gnuradio.__path__:
    gnuradio.__path__.append(_py)
