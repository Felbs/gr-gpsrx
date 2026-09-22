/* -*- c++ -*- */
/*
 * Copyright 2026 gr-gpsrx authors.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_GPSRX_CHANNEL_CC_IMPL_H
#define INCLUDED_GPSRX_CHANNEL_CC_IMPL_H

#include <gnuradio/gpsrx/channel_cc.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>
#include <mutex>

namespace gr {
namespace gpsrx {

// python/gpsrx/track.py, in C++. One code period per step().
class tracker
{
public:
    static constexpr double CODE_RATE = 1.023e6;
    static constexpr int CODE_LEN = 1023;
    static constexpr double L1_HZ = 1575.42e6;
    static constexpr double CARRIER_TO_CODE = CODE_RATE / L1_HZ;

    tracker(int prn, double fs, double doppler_hz, double code_phase_samples, double pll_bw, double dll_bw,
            double pll_bw_narrow = 15.0, double dll_bw_narrow = 0.5, int coherent_ms = 20);

    int samples_needed() const;
    // exactly samples_needed() samples in; the prompt out; state advanced
    std::complex<double> step(const std::complex<float>* x, int n);

    // state, named as in ChannelState
    int prn;
    double fs;
    double carrier_hz;
    double code_phase;
    int64_t epochs = 0;
    double carrier_phase = 0.0;
    double code_rate = CODE_RATE;
    int64_t samples_in = 0;
    double pll_e_prev = 0.0, dll_e_prev = 0.0;
    double lock = 0.0;
    double cn0_db = 0.0;            // dB-Hz, moment estimator over 20 prompts
    double carrier_cycles = 0.0;    // accumulated NCO phase in cycles: the carrier-phase observable
    int bit_offset = -1;            // period index (mod 20) at which a data bit begins; -1 = stage 1
    double epoch_sample() const { return samples_in - code_phase * fs / code_rate; }
    double epoch_cycles() const { return carrier_cycles - carrier_hz * (code_phase * fs / code_rate) / fs; }

private:
    std::array<int8_t, CODE_LEN> code_;
    double pll_t1_, pll_t2_, dll_t1_, dll_t2_;
    double doppler0_, carr_corr_ = 0.0, code_dop_ = 0.0;
    int n_nominal_;
    // stage 2 (track.py): bit sync -> narrow loops + coherent integration over a whole bit
    double pll_bw_narrow_, dll_bw_narrow_;
    int coh_;
    int64_t flips_[20] = { 0 };
    double last_ip_ = 0.0, f_avg_ = 0.0;
    bool have_f_avg_ = false, aligned_ = false;
    std::complex<double> acc_[3] = { 0, 0, 0 };
    double acc_dt_ = 0.0;
    double m2_ = 0.0, m4_ = 0.0;
    int mn_ = 0;
    void bit_sync(double ip);
    static constexpr double SPACING = 0.5;
};

class channel_cc_impl : public channel_cc
{
private:
    double fs_;
    int slot_;
    double pll_bw_, dll_bw_, pll_bw_narrow_, dll_bw_narrow_;
    int coherent_ms_;
    int obs_every_, batch_;
    std::mutex mtx_;
    int prn_ = 0;
    bool have_pending_ = false;
    int pend_prn_ = 0;
    double pend_dop_ = 0.0, pend_sample_ = 0.0;
    std::unique_ptr<tracker> eng_;
    int64_t t0_abs_ = 0;
    int64_t n_periods_ = 0;
    int lost_run_ = 0;
    static constexpr double LOST_LOCK = 0.2;
    static constexpr int LOST_PERIODS = 500;

    void on_assign(pmt::pmt_t msg);
    void publish_status(const char* what, pmt::pmt_t extra = pmt::PMT_NIL);
    int nominal() const { return (int)std::lround(fs_ * 1e-3); }

public:
    channel_cc_impl(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                    double pll_bw_narrow, double dll_bw_narrow, int coherent_ms);
    ~channel_cc_impl() override;

    void forecast(int noutput_items, gr_vector_int& ninput_items_required) override;
    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items) override;
};

} // namespace gpsrx
} // namespace gr

#endif /* INCLUDED_GPSRX_CHANNEL_CC_IMPL_H */
