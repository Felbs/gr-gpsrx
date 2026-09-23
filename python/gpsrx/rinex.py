# SPDX-License-Identifier: GPL-3.0-or-later
"""RINEX 3.04 output: the observables (C1C, L1C, D1C, S1C for GPS and Galileo) and the decoded
ephemerides (GPS LNAV, Galileo I/NAV), so anyone can process what this receiver measured with
RTKLIB, gnss-sdr's tools or their own - and check it against ours.

What goes in: per solve, `fix["observables"]` (pvt.fix_from_channels): the raw code pseudorange
with SV clock, ionosphere and troposphere still in it, the accumulated carrier phase in cycles,
the Doppler and C/N0, all referred to the solve's receive instant t_rx (GPS time; the receiver
clock term is therefore ~0, which RINEX allows - the clock offset field is optional). Epochs are
written in GPS time; the file's `TIME OF FIRST OBS` says so. Galileo satellites are E records
with the same four observables (E1-C pilot code and carrier, in RINEX's terms C1C/L1C on E1).

Nothing here is a position: the APPROX POSITION line is left at zero on purpose."""
import os
import time

import numpy as np

GPS_EPOCH_JD = 2444244.5             # 1980-01-06 00:00 UTC as a Julian day
OBS_TYPES = ("C1C", "L1C", "D1C", "S1C")


def full_week(wn_mod1024, year_hint=None):
    """The GPS week from the broadcast 10-bit week: the 2019-04-07 rollover is the current era
    (weeks 2048..3071); a receiver run before 2038 needs no other assumption."""
    return int(wn_mod1024) + 2048


def gps_to_calendar(week, tow):
    """GPS week + time of week (s) -> (year, month, day, hour, minute, second) in GPS time
    (no leap seconds applied: RINEX epochs for GPS/Galileo are in GPS time)."""
    days = week * 7 + int(tow // 86400)
    jd = GPS_EPOCH_JD + days + 0.5                       # noon-based day number for the calendar step
    a = int(jd + 0.5)
    b = a + 1537
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    sod = tow - int(tow // 86400) * 86400
    hour = int(sod // 3600)
    minute = int((sod - hour * 3600) // 60)
    second = sod - hour * 3600 - minute * 60
    return year, month, day, hour, minute, second


def _f(v, width=14, dec=3):
    return f"{v:{width}.{dec}f}"


def _e19(v):
    """RINEX nav field: D19.12 with a two-digit exponent, Fortran style."""
    s = f"{float(v):19.12E}"
    mant, exp = s.split("E")
    return f"{mant}D{int(exp):+03d}"


class ObsWriter:
    """RINEX 3.04 observation file, written epoch by epoch (header first, data appended)."""

    def __init__(self, path, marker="GPSRX", program="gr-gpsrx"):
        self.path = path
        self.marker = marker
        self.program = program
        self.n = 0
        self._f = None

    def _open(self, week, tow):
        os.makedirs(os.path.dirname(os.path.abspath(self.path)) or ".", exist_ok=True)
        self._f = open(self.path, "w", newline="\n")
        y, mo, d, h, mi, s = gps_to_calendar(week, tow)
        w = self._f.write
        w(f"{3.04:9.2f}           {'OBSERVATION DATA':<20}{'M':<20}RINEX VERSION / TYPE\n")
        w(f"{self.program:<20}{'':<20}{time.strftime('%Y%m%d %H%M%S UTC'):<20}PGM / RUN BY / DATE\n")
        w(f"{self.marker:<60}MARKER NAME\n")
        w(f"{'NON_PHYSICAL':<60}MARKER TYPE\n")
        w(f"{'':<20}{'':<40}OBSERVER / AGENCY\n")
        w(f"{'0':<20}{'gr-gpsrx':<20}{'':<20}REC # / TYPE / VERS\n")
        w(f"{'0':<20}{'UNKNOWN':<20}{'':<20}ANT # / TYPE\n")
        w(f"{0.0:14.4f}{0.0:14.4f}{0.0:14.4f}{'':<18}APPROX POSITION XYZ\n")
        w(f"{0.0:14.4f}{0.0:14.4f}{0.0:14.4f}{'':<18}ANTENNA: DELTA H/E/N\n")
        for sysc in ("G", "E"):
            w(f"{sysc}  {len(OBS_TYPES):3d}" + "".join(f" {t}" for t in OBS_TYPES) + " " * (60 - 6 - 4 * len(OBS_TYPES)) + "SYS / # / OBS TYPES\n")
        w(f"{y:6d}{mo:6d}{d:6d}{h:6d}{mi:6d}{s:13.7f}     {'GPS':<3}         TIME OF FIRST OBS\n")
        w(f"{'DBHZ':<60}SIGNAL STRENGTH UNIT\n")
        w(f"{'':<60}END OF HEADER\n")

    def add_epoch(self, week, tow, observables):
        """observables: [{sys, prn, pr_m, phase_cyc, doppler_hz, cn0_db}] at GPS time (week, tow)."""
        if not observables:
            return
        if self._f is None:
            self._open(week, tow)
        y, mo, d, h, mi, s = gps_to_calendar(week, tow)
        w = self._f.write
        w(f"> {y:4d} {mo:02d} {d:02d} {h:02d} {mi:02d}{s:11.7f}  0{len(observables):3d}\n")
        for o in sorted(observables, key=lambda o: (o["sys"] != "GPS", o["prn"])):
            sysc = "E" if o["sys"] == "GAL" else "G"
            w(f"{sysc}{o['prn']:02d}{_f(o['pr_m'])}  {_f(o['phase_cyc'])}  {_f(o['doppler_hz'])}  {_f(o['cn0_db'])}  \n")
        self.n += 1
        self._f.flush()

    def close(self):
        if self._f is not None:
            self._f.close()
            self._f = None


def write_nav(path, ephs, program="gr-gpsrx"):
    """RINEX 3.04 navigation file (mixed) from a list of ephemeris dicts as the decoders make
    them: GPS (nav.py keys) and Galileo (nav_gal.py keys; af0 has BGD(E1,E5b) folded in for the
    single-frequency user - it is unfolded here, RINEX carries the broadcast value plus the BGDs)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="\n") as f:
        w = f.write
        w(f"{3.04:9.2f}           {'N: GNSS NAV DATA':<20}{'M: MIXED':<20}RINEX VERSION / TYPE\n")
        w(f"{program:<20}{'':<20}{time.strftime('%Y%m%d %H%M%S UTC'):<20}PGM / RUN BY / DATE\n")
        w(f"{'':<60}END OF HEADER\n")
        for e in ephs:
            if e.get("sys", "GPS") == "GAL":
                _write_gal(w, e)
            else:
                _write_gps(w, e)


def _write_gps(w, e):
    need = ("WN", "af0", "af1", "af2", "toc", "Crs", "dn", "M0", "Cuc", "e", "Cus", "sqrtA", "toe",
            "Cic", "Omega0", "Cis", "i0", "Crc", "omega", "OmegaDot", "IDOT")
    if any(k not in e for k in need):
        return
    week = full_week(e["WN"])
    y, mo, d, h, mi, s = gps_to_calendar(week, e["toc"])
    iode = e.get("IODE2", e.get("IODE3", 0))
    w(f"G{int(e['prn']):02d} {y:4d} {mo:02d} {d:02d} {h:02d} {mi:02d} {int(s):02d}{_e19(e['af0'])}{_e19(e['af1'])}{_e19(e['af2'])}\n")
    w(f"    {_e19(iode)}{_e19(e['Crs'])}{_e19(e['dn'])}{_e19(e['M0'])}\n")
    w(f"    {_e19(e['Cuc'])}{_e19(e['e'])}{_e19(e['Cus'])}{_e19(e['sqrtA'])}\n")
    w(f"    {_e19(e['toe'])}{_e19(e['Cic'])}{_e19(e['Omega0'])}{_e19(e['Cis'])}\n")
    w(f"    {_e19(e['i0'])}{_e19(e['Crc'])}{_e19(e['omega'])}{_e19(e['OmegaDot'])}\n")
    w(f"    {_e19(e['IDOT'])}{_e19(0)}{_e19(week)}{_e19(0)}\n")                     # codes on L2, week, L2 P flag
    w(f"    {_e19(_ura_m(e.get('URA', 0)))}{_e19(e.get('health', 0))}{_e19(e.get('TGD', 0.0))}{_e19(e.get('IODC', iode))}\n")
    w(f"    {_e19(e['toe'] - 3600.0)}{_e19(4)}\n")                                  # transmission time (approx), fit interval


def _write_gal(w, e):
    need = ("toe", "M0", "e", "sqrtA", "Omega0", "i0", "omega", "IDOT", "OmegaDot", "dn", "Cuc", "Cus",
            "Crc", "Crs", "toc", "af0", "af1", "af2", "Cic", "Cis", "IODnav")
    if any(k not in e for k in need) or e.get("WN") is None:
        return
    week = int(e["WN"]) + 1024                          # Galileo week numbers count from the GPS 1980 epoch + 1024
    y, mo, d, h, mi, s = gps_to_calendar(week, e["toc"])
    af0 = e["af0"] + (e.get("BGD_E1E5b", 0.0) if e.get("bgd_folded") else 0.0)
    w(f"E{int(e['prn']):02d} {y:4d} {mo:02d} {d:02d} {h:02d} {mi:02d} {int(s):02d}{_e19(af0)}{_e19(e['af1'])}{_e19(e['af2'])}\n")
    w(f"    {_e19(e['IODnav'])}{_e19(e['Crs'])}{_e19(e['dn'])}{_e19(e['M0'])}\n")
    w(f"    {_e19(e['Cuc'])}{_e19(e['e'])}{_e19(e['Cus'])}{_e19(e['sqrtA'])}\n")
    w(f"    {_e19(e['toe'])}{_e19(e['Cic'])}{_e19(e['Omega0'])}{_e19(e['Cis'])}\n")
    w(f"    {_e19(e['i0'])}{_e19(e['Crc'])}{_e19(e['omega'])}{_e19(e['OmegaDot'])}\n")
    w(f"    {_e19(e['IDOT'])}{_e19(513)}{_e19(week)}{_e19(0)}\n")                   # data source 513 = I/NAV E1-B, af0 wrt E5b/E1
    # health, RINEX packing: E1-B DVS bit 0, E1-B HS bits 1-2, E5b DVS bit 3, E5b HS bits 4-5 (E5a: not in I/NAV)
    health = (int(e.get("E1BDVS", 0)) | (int(e.get("E1BHS", 0)) << 1) | (int(e.get("E5bDVS", 0)) << 3) | (int(e.get("E5bHS", 0)) << 4))
    w(f"    {_e19(_sisa_m(e.get('SISA', 255)))}{_e19(health)}{_e19(e.get('BGD_E1E5a', 0.0))}{_e19(e.get('BGD_E1E5b', 0.0))}\n")
    w(f"    {_e19(e['toe'] - 3600.0)}{_e19(0)}\n")


def _ura_m(idx):
    """The URA index as metres, RTKLIB's (and BKG's BNC's) convention 2^(1 + N/2) for N <= 6, 2^(N - 2)
    above - what the IGS broadcast files carry, so the field compares equal (checked against BRDC)."""
    i = int(idx)
    if i < 0 or i > 15:
        return -1.0
    return float(2 ** (1 + i / 2.0)) if i <= 6 else float(2 ** (i - 2))


def _sisa_m(idx):
    i = int(idx)
    if i < 50:
        return i * 0.01
    if i < 75:
        return 0.5 + (i - 50) * 0.02
    if i < 100:
        return 1.0 + (i - 75) * 0.04
    if i < 126:
        return 2.0 + (i - 100) * 0.16
    return -1.0


def read_obs(path):
    """A small reader for the tests: epochs of [(sys, prn, C1C, L1C, D1C, S1C)] with (week-less)
    GPS time of day from the epoch lines. Not a general RINEX parser."""
    epochs = []
    with open(path) as f:
        for line in f:
            if line.rstrip().endswith("END OF HEADER"):
                break
        cur = None
        for line in f:
            if line.startswith(">"):
                p = line.split()
                sod = int(p[4]) * 3600 + int(p[5]) * 60 + float(p[6])
                cur = (sod, [])
                epochs.append(cur)
            elif cur is not None and line[:1] in "GE":
                vals = [float(line[3 + 16 * i: 3 + 16 * i + 14]) for i in range(4)]
                cur[1].append((line[0], int(line[1:3]), *vals))
    return epochs
