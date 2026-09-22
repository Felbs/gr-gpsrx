# SPDX-License-Identifier: GPL-3.0-or-later
"""Run from the source tree with nothing installed: grafts python/ onto the gnuradio package,
and if the C++ blocks have been built (build/cpp), makes their extension module importable
too - the .pyd/.so is copied beside the Python package and the shared library's directory is
registered, so `from gnuradio import gpsrx` gives both the Python and the C++ blocks."""
import glob
import os
import shutil

import gnuradio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_py = os.path.join(ROOT, "python")
if _py not in gnuradio.__path__:
    gnuradio.__path__.append(_py)

_build = os.path.join(ROOT, "build", "cpp")
_ext = glob.glob(os.path.join(_build, "python", "gpsrx", "bindings", "gpsrx_python*.pyd")) + \
    glob.glob(os.path.join(_build, "python", "gpsrx", "bindings", "gpsrx_python*.so"))
if _ext:
    dst = os.path.join(_py, "gpsrx", os.path.basename(_ext[0]))
    if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(_ext[0]):
        shutil.copy2(_ext[0], dst)
    _lib = os.path.join(_build, "lib")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_lib)                          # Windows: the DLL beside the build
    else:
        os.environ["LD_LIBRARY_PATH"] = _lib + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")
