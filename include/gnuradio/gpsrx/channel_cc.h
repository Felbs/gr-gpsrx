/* -*- c++ -*- */
/*
 * Copyright 2026 gr-gpsrx authors.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_GPSRX_CHANNEL_CC_H
#define INCLUDED_GPSRX_CHANNEL_CC_H

#include <gnuradio/block.h>
#include <gnuradio/gpsrx/api.h>

namespace gr {
namespace gpsrx {

/*!
 * \brief One satellite's tracking loop, in C++: the same Channel as python/gpsrx/channel.py,
 * for when eight of them must keep up with a radio.
 * \ingroup gpsrx
 *
 * Complex baseband in; one prompt correlation out per code period (1 kHz). Message ports:
 * 'assign' in (from Acquisition: slot, prn, doppler_hz, sample), 'obs' out once a second
 * (epochs, code_phase, carrier_hz, samples_in, epoch_sample, lock, sample_abs, slot),
 * 'status' out (tracking / lost / idle). A 'gpsrx_assign' tag marks the first prompt of each
 * assignment. Behaviour, ports, tags and messages are those of the Python block, so a Nav
 * Decoder or PVT cannot tell which one it is wired to; the arithmetic is python/gpsrx/track.py
 * line for line (NCO, E/P/L, Costas atan(Q/I), early-late DLL, carrier aiding, lock detector).
 */
class GPSRX_API channel_cc : virtual public gr::block
{
public:
    typedef std::shared_ptr<channel_cc> sptr;

    /*!
     * \param samp_rate input sample rate (2.048 MS/s tested)
     * \param slot this channel's number; Acquisition assigns by it
     * \param pll_bw Costas loop noise bandwidth, Hz
     * \param dll_bw code loop noise bandwidth, Hz
     * \param obs_every_ms code periods between 'obs' messages
     * \param pll_bw_narrow PLL bandwidth after bit sync (0 = stay in stage 1); 5 Hz static, 8-15 Hz moving
     * \param dll_bw_narrow DLL bandwidth after bit sync
     * \param coherent_ms coherent integration after bit sync, a divisor of 20
     */
    static sptr make(double samp_rate, int slot, double pll_bw = 18.0, double dll_bw = 2.0, int obs_every_ms = 1000,
                     double pll_bw_narrow = 15.0, double dll_bw_narrow = 0.5, int coherent_ms = 20);
};

} // namespace gpsrx
} // namespace gr

#endif /* INCLUDED_GPSRX_CHANNEL_CC_H */
