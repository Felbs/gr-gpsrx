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
#include <string>
#include <vector>
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

    // signal: the code (chips, +-1), an optional data code (a pilot-aided signal: the loops run on
    // `code`, the data prompt comes from `data_code`), BOC(1,1) subcarrier or not, periods per
    // data bit (20 GPS; 25 = the pilot's secondary code length; 1 = no bit sync), the secondary
    // code (+-1 per period, known) and the correlator spacing. Empty code = GPS L1 C/A from prn.
    struct signal {
        std::vector<int8_t> code, data_code, secondary;
        bool boc = false;
        int bit_periods = 20;
        double spacing = 0.5;
        std::string name = "L1CA";
    };
    tracker(int prn, double fs, double doppler_hz, double code_phase_samples, double pll_bw, double dll_bw,
            double pll_bw_narrow = 15.0, double dll_bw_narrow = 0.5, int coherent_ms = 20, int pll_order = 3,
            const signal* sig = nullptr);

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
    int code_len = CODE_LEN;        // 1023 GPS, 4092 Galileo
    double period_s() const { return code_len / CODE_RATE; }
    std::string sys = "GPS";
    double epoch_sample() const { return samples_in - code_phase * fs / code_rate; }
    double epoch_cycles() const { return carrier_cycles - carrier_hz * (code_phase * fs / code_rate) / fs; }

private:
    std::vector<int8_t> code_, data_code_, secondary_;
    bool boc_ = false;
    int bit_periods_ = 20;
    double spacing_ = 0.5;
    std::vector<int8_t> sec_hist_;
    double sec_pol_ = 1.0;
    inline double chip(const std::vector<int8_t>& c, double phase_chips) const;
    void secondary_sync(double ip);
    void switch_to_narrow();
    double pll_t1_, pll_t2_, dll_t1_, dll_t2_;
    double doppler0_, carr_corr_ = 0.0, code_dop_ = 0.0;
    int n_nominal_;
    // stage 2 (track.py): bit sync -> narrow loops + coherent integration over a whole bit
    double pll_bw_narrow_, dll_bw_narrow_;
    int coh_;
    int64_t flips_[32] = { 0 };
    double last_ip_ = 0.0, f_avg_ = 0.0, cd_avg_ = 0.0;
    int pll_order_;
    double w3_ = 0.0, acc3_ = 0.0, vel3_ = 0.0; // third-order loop (narrow stage), Kaplan & Hegarty form
    bool have_f_avg_ = false, aligned_ = false;
    std::complex<double> acc_[3] = { 0, 0, 0 };
    double acc_dt_ = 0.0;
    double m2_ = 0.0, m4_ = 0.0;
    int mn_ = 0;
    void bit_sync(double ip);
};

class channel_cc_impl : public channel_cc
{
private:
    double fs_;
    int slot_;
    std::string signal_;
    tracker::signal sig_;             // built at assignment from signal_ (Galileo: the vendored tables)
    double pll_bw_, dll_bw_, pll_bw_narrow_, dll_bw_narrow_;
    int coherent_ms_, pll_order_;
    int obs_every_, batch_, obs_periods_ = 1000;
    std::mutex mtx_;
    int prn_ = 0;
    bool have_pending_ = false;
    int pend_prn_ = 0;
    double pend_dop_ = 0.0, pend_sample_ = 0.0;
    int64_t start_at_ = 0;
    std::unique_ptr<tracker> eng_;
    int64_t t0_abs_ = 0;
    int64_t n_periods_ = 0;
    int lost_run_ = 0;
    static constexpr double LOST_LOCK = 0.2;
    static constexpr int LOST_PERIODS = 500;

    void on_assign(pmt::pmt_t msg);
    void publish_status(const char* what, pmt::pmt_t extra = pmt::PMT_NIL);
    int nominal() const { return (int)std::lround(fs_ * (signal_ == "L1CA" ? 1e-3 : 4e-3)); }

public:
    channel_cc_impl(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                    double pll_bw_narrow, double dll_bw_narrow, int coherent_ms, int pll_order,
                    const std::string& signal);
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
