/* -*- c++ -*- */
/*
 * Copyright 2026 gr-gpsrx authors.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include "channel_cc_impl.h"
#include <gnuradio/io_signature.h>

#include <algorithm>
#include <cstring>

namespace gr {
namespace gpsrx {

// ---- the C/A code (python/gpsrx/cacode.py): IS-GPS-200 G2 phase-select taps --------------
static const int G2_TAPS[33][2] = {
    { 0, 0 },  { 2, 6 },  { 3, 7 },  { 4, 8 },  { 5, 9 },  { 1, 9 },  { 2, 10 }, { 1, 8 },  { 2, 9 },
    { 3, 10 }, { 2, 3 },  { 3, 4 },  { 5, 6 },  { 6, 7 },  { 7, 8 },  { 8, 9 },  { 9, 10 }, { 1, 4 },
    { 2, 5 },  { 3, 6 },  { 4, 7 },  { 5, 8 },  { 6, 9 },  { 1, 3 },  { 4, 6 },  { 5, 7 },  { 6, 8 },
    { 7, 9 },  { 8, 10 }, { 1, 6 },  { 2, 7 },  { 3, 8 },  { 4, 9 }
};

static void ca_code(int prn, std::array<int8_t, tracker::CODE_LEN>& out)
{
    int g1[10], g2[10];
    for (int i = 0; i < 10; i++)
        g1[i] = g2[i] = 1;
    const int s1 = G2_TAPS[prn][0] - 1, s2 = G2_TAPS[prn][1] - 1;
    for (int i = 0; i < tracker::CODE_LEN; i++) {
        const int bit = g1[9] ^ g2[s1] ^ g2[s2];
        out[i] = (int8_t)(1 - 2 * bit);
        const int fb1 = g1[2] ^ g1[9];
        const int fb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];
        for (int k = 9; k > 0; k--) {
            g1[k] = g1[k - 1];
            g2[k] = g2[k - 1];
        }
        g1[0] = fb1;
        g2[0] = fb2;
    }
}

// track.loop_gains: SoftGNSS / Kaplan & Hegarty second-order loop filter constants
static void loop_gains(double bw_hz, double k, double& t1, double& t2)
{
    const double zeta = 0.707;
    const double wn = bw_hz * 8 * zeta / (4 * zeta * zeta + 1);
    t1 = k / (wn * wn);
    t2 = 2 * zeta / wn;
}

static inline int chip_index(double phase_chips)
{
    const int i = (int)std::floor(phase_chips) % tracker::CODE_LEN;
    return i < 0 ? i + tracker::CODE_LEN : i;
}

tracker::tracker(int prn_, double fs_, double doppler_hz, double code_phase_samples, double pll_bw, double dll_bw,
                 double pll_bw_narrow, double dll_bw_narrow, int coherent_ms)
    : prn(prn_), fs(fs_), carrier_hz(doppler_hz), doppler0_(doppler_hz),
      pll_bw_narrow_(pll_bw_narrow), dll_bw_narrow_(dll_bw_narrow), coh_(pll_bw_narrow > 0 ? coherent_ms : 1)
{
    ca_code(prn, code_);
    // acquisition hands over the sample at which the code starts; the phase at sample 0 is
    // therefore -(that) chips/sample, modulo 1023
    const double cps = CODE_RATE / fs;
    code_phase = std::fmod(-code_phase_samples * cps, (double)CODE_LEN);
    if (code_phase < 0)
        code_phase += CODE_LEN;
    loop_gains(pll_bw, 0.25, pll_t1_, pll_t2_); // atan/2pi is +-1/4 cycle full scale
    loop_gains(dll_bw, 1.0, dll_t1_, dll_t2_);
    n_nominal_ = (int)std::lround(fs * CODE_LEN / CODE_RATE);
}

int tracker::samples_needed() const
{
    const double chips_left = CODE_LEN - code_phase;
    const int n = (int)std::ceil(chips_left * fs / code_rate);
    return n < 1 ? 1 : n;
}

std::complex<double> tracker::step(const std::complex<float>* x, int n)
{
    // 1. carrier wipe-off: a rotating phasor at the current Doppler from the current phase
    //    (the Python engine builds the same ramp as a vector; here it is one multiply per sample)
    const double w = -2.0 * M_PI * carrier_hz / fs;
    std::complex<double> ph(std::cos(-carrier_phase), std::sin(-carrier_phase));
    const std::complex<double> stepc(std::cos(w), std::sin(w));
    // 2. early / prompt / late replicas half a chip apart, nearest-chip lookup, accumulated
    const double cps = code_rate / fs;
    std::complex<double> E = 0, P = 0, L = 0;
    double phc = code_phase;
    for (int k = 0; k < n; k++) {
        const std::complex<double> xb = std::complex<double>(x[k].real(), x[k].imag()) * ph;
        ph *= stepc;
        E += xb * (double)code_[chip_index(phc + SPACING)];
        P += xb * (double)code_[chip_index(phc)];
        L += xb * (double)code_[chip_index(phc - SPACING)];
        phc += cps;
    }
    // advance phases by what this period consumed
    const double dt = n / fs;
    carrier_phase = std::fmod(carrier_phase + 2.0 * M_PI * carrier_hz * dt, 2.0 * M_PI);
    if (carrier_phase < 0)
        carrier_phase += 2.0 * M_PI;
    code_phase = std::fmod(code_phase + n * cps, (double)CODE_LEN);
    samples_in += n;
    epochs += 1;
    const double ip = P.real(), qp = P.imag();
    // lock indicator (I^2 - Q^2) / (I^2 + Q^2), smoothed over ~50 periods
    const double p = ip * ip + qp * qp;
    const double li = p > 0 ? (ip * ip - qp * qp) / p : 0.0;
    lock = 0.98 * lock + 0.02 * li;
    // stage 1: bit sync (where, mod 20, does the prompt's sign flip?) and the averaged frequency
    if (bit_offset < 0 && coh_ > 1) {
        f_avg_ = have_f_avg_ ? 0.99 * f_avg_ + 0.01 * carrier_hz : carrier_hz;
        have_f_avg_ = true;
        bit_sync(ip);
    }
    // stage 2: coherent window over a whole bit, one loop update per window (track.py)
    std::complex<double> Ew = E, Pw = P, Lw = L;
    double dt_loop = dt;
    if (bit_offset >= 0) {
        const bool at_edge = ((epochs - bit_offset) % coh_ + coh_) % coh_ == 0;
        if (!aligned_) { // no loop update until the first bit edge (a partial window kicks the loop)
            aligned_ = at_edge;
            return P;
        }
        acc_[0] += E; acc_[1] += P; acc_[2] += L;
        acc_dt_ += dt;
        if (!at_edge)
            return P;
        Ew = acc_[0]; Pw = acc_[1]; Lw = acc_[2];
        dt_loop = acc_dt_;
        acc_[0] = acc_[1] = acc_[2] = 0;
        acc_dt_ = 0.0;
    }
    const double ipw = Pw.real(), qpw = Pw.imag();
    // 3. PLL: Costas discriminator atan(Q/I) (NOT atan2: blind to the data flips), error in
    //    cycles, loop filter -> frequency correction; the NCO phase accumulation integrates it
    const double e_pll = ipw != 0.0 ? std::atan(qpw / ipw) / (2.0 * M_PI) : 0.0;
    carr_corr_ += (pll_t2_ * (e_pll - pll_e_prev) + dt_loop * e_pll) / pll_t1_;
    pll_e_prev = e_pll;
    carrier_hz = doppler0_ + carr_corr_;
    // 4. DLL: normalised early-minus-late envelope, carrier-aided code rate
    const double Em = std::abs(Ew), Lm = std::abs(Lw);
    const double e_dll = (Em + Lm) > 0 ? 0.5 * (Em - Lm) / (Em + Lm) : 0.0;
    code_dop_ += (dll_t2_ * (e_dll - dll_e_prev) + dt_loop * e_dll) / dll_t1_;
    dll_e_prev = e_dll;
    code_rate = CODE_RATE + carrier_hz * CARRIER_TO_CODE + code_dop_;
    return P;
}

void tracker::bit_sync(double ip)
{
    // this period's index is epochs - 1 (epochs already counts it)
    if (last_ip_ != 0.0 && (ip < 0) != (last_ip_ < 0) && lock > 0.5)
        flips_[(epochs - 1) % 20]++;
    last_ip_ = ip;
    if (epochs < 300 || lock <= 0.8)
        return;
    int64_t top = 0, second = 0, arg = 0;
    for (int i = 0; i < 20; i++) {
        if (flips_[i] > top) { second = top; top = flips_[i]; arg = i; }
        else if (flips_[i] > second) second = flips_[i];
    }
    if (top >= 8 && top >= 3 * std::max<int64_t>(second, 1)) {
        bit_offset = (int)arg;
        // k = 1 for the narrow stage (see track.py: with k = 0.25 the 20 ms loop is over unity gain)
        loop_gains(pll_bw_narrow_, 1.0, pll_t1_, pll_t2_);
        loop_gains(dll_bw_narrow_, 1.0, dll_t1_, dll_t2_);
        pll_e_prev = dll_e_prev = 0.0;
        acc_[0] = acc_[1] = acc_[2] = 0;
        acc_dt_ = 0.0;
        aligned_ = false;
        if (have_f_avg_) { // start from the averaged frequency, not the jittering instantaneous one
            carr_corr_ = f_avg_ - doppler0_;
            carrier_hz = f_avg_;
        }
    }
}

// ---- the block ------------------------------------------------------------------------------
channel_cc::sptr channel_cc::make(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                                  double pll_bw_narrow, double dll_bw_narrow, int coherent_ms)
{
    return gnuradio::make_block_sptr<channel_cc_impl>(samp_rate, slot, pll_bw, dll_bw, obs_every_ms,
                                                      pll_bw_narrow, dll_bw_narrow, coherent_ms);
}

channel_cc_impl::channel_cc_impl(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                                 double pll_bw_narrow, double dll_bw_narrow, int coherent_ms)
    : gr::block("gpsrx_channel_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      fs_(samp_rate),
      slot_(slot),
      pll_bw_(pll_bw),
      dll_bw_(dll_bw),
      pll_bw_narrow_(pll_bw_narrow),
      dll_bw_narrow_(dll_bw_narrow),
      coherent_ms_(coherent_ms),
      obs_every_(obs_every_ms),
      batch_(3)
{
    message_port_register_in(pmt::mp("assign"));
    set_msg_handler(pmt::mp("assign"), [this](pmt::pmt_t msg) { this->on_assign(msg); });
    message_port_register_out(pmt::mp("obs"));
    message_port_register_out(pmt::mp("status"));
    // the same batching as the Python block: whole batches or nothing (see channel.py)
    set_output_multiple(batch_);
    set_min_output_buffer(16 * batch_);
}

channel_cc_impl::~channel_cc_impl() {}

static double num(pmt::pmt_t d, const char* key, double dflt)
{
    pmt::pmt_t v = pmt::dict_ref(d, pmt::mp(key), pmt::PMT_NIL);
    if (pmt::is_integer(v))
        return (double)pmt::to_long(v);
    if (pmt::is_real(v))
        return pmt::to_double(v);
    return dflt;
}

void channel_cc_impl::on_assign(pmt::pmt_t msg)
{
    if (!pmt::is_dict(msg))
        return;
    if (pmt::dict_has_key(msg, pmt::mp("slot")) && (int)num(msg, "slot", -1) != slot_)
        return; // every channel hears every assignment
    std::lock_guard<std::mutex> g(mtx_);
    const int prn = (int)num(msg, "prn", 0);
    if (prn <= 0) {
        eng_.reset();
        prn_ = 0;
        have_pending_ = false;
        publish_status("idle");
        return;
    }
    pend_prn_ = prn;
    pend_dop_ = num(msg, "doppler_hz", 0.0);
    pend_sample_ = num(msg, "sample", 0.0);
    have_pending_ = true;
    eng_.reset();
    prn_ = prn;
}

void channel_cc_impl::publish_status(const char* what, pmt::pmt_t extra)
{
    pmt::pmt_t d = pmt::make_dict();
    d = pmt::dict_add(d, pmt::mp("slot"), pmt::from_long(slot_));
    d = pmt::dict_add(d, pmt::mp("prn"), pmt::from_long(prn_));
    d = pmt::dict_add(d, pmt::mp("what"), pmt::string_to_symbol(what));
    if (pmt::is_dict(extra)) {
        pmt::pmt_t items = pmt::dict_items(extra);
        while (pmt::is_pair(items)) {
            pmt::pmt_t kv = pmt::car(items);
            d = pmt::dict_add(d, pmt::car(kv), pmt::cdr(kv));
            items = pmt::cdr(items);
        }
    }
    message_port_pub(pmt::mp("status"), d);
}

void channel_cc_impl::forecast(int noutput_items, gr_vector_int& ninput_items_required)
{
    // a period is ~1 ms of samples, and can run a few samples long while the DLL pulls in
    ninput_items_required[0] = noutput_items * (nominal() + 8);
}

int channel_cc_impl::general_work(int noutput_items,
                                  gr_vector_int& ninput_items,
                                  gr_vector_const_void_star& input_items,
                                  gr_vector_void_star& output_items)
{
    const gr_complex* x = (const gr_complex*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];
    const int n_in = ninput_items[0];
    std::lock_guard<std::mutex> g(mtx_);
    if (prn_ == 0) {
        consume(0, n_in); // idle: eat the input, produce nothing
        return 0;
    }
    if (!eng_) {
        if (!have_pending_) {
            consume(0, n_in);
            return 0;
        }
        const int64_t start = nitems_read(0);
        // acquisition's code start is `sample` absolute, seconds old by now: wrap it forward
        // with the DOPPLER-SHIFTED code period (nominal would be 3 chips/s wrong at 5 kHz)
        const double period = fs_ * tracker::CODE_LEN / (tracker::CODE_RATE * (1.0 + pend_dop_ / tracker::L1_HZ));
        double rel = std::fmod(pend_sample_ - (double)start, period);
        if (rel < 0)
            rel += period;
        eng_.reset(new tracker(pend_prn_, fs_, pend_dop_, rel, pll_bw_, dll_bw_, pll_bw_narrow_, dll_bw_narrow_, coherent_ms_));
        t0_abs_ = start;
        lost_run_ = 0;
        have_pending_ = false;
        pmt::pmt_t tag = pmt::make_dict();
        tag = pmt::dict_add(tag, pmt::mp("prn"), pmt::from_long(pend_prn_));
        tag = pmt::dict_add(tag, pmt::mp("slot"), pmt::from_long(slot_));
        add_item_tag(0, nitems_written(0), pmt::mp("gpsrx_assign"), tag);
        pmt::pmt_t ex = pmt::make_dict();
        ex = pmt::dict_add(ex, pmt::mp("doppler_hz"), pmt::from_double(pend_dop_));
        publish_status("tracking", ex);
    }
    int consumed = 0, produced = 0;
    const int want = std::min(noutput_items, batch_);
    if (n_in < want * (nominal() + 8)) {
        consume(0, 0);
        return 0;
    }
    while (produced < want) {
        const int need = eng_->samples_needed();
        if (consumed + need > n_in)
            break;
        const std::complex<double> P = eng_->step(x + consumed, need);
        out[produced] = gr_complex((float)P.real(), (float)P.imag());
        produced++;
        consumed += need;
        n_periods_++;
        if (eng_->lock < LOST_LOCK)
            lost_run_++;
        else
            lost_run_ = 0;
        if (n_periods_ % obs_every_ == 0) {
            pmt::pmt_t d = pmt::make_dict();
            d = pmt::dict_add(d, pmt::mp("prn"), pmt::from_long(eng_->prn));
            d = pmt::dict_add(d, pmt::mp("epochs"), pmt::from_long((long)eng_->epochs));
            d = pmt::dict_add(d, pmt::mp("code_phase"), pmt::from_double(eng_->code_phase));
            d = pmt::dict_add(d, pmt::mp("carrier_hz"), pmt::from_double(eng_->carrier_hz));
            d = pmt::dict_add(d, pmt::mp("samples_in"), pmt::from_long((long)eng_->samples_in));
            d = pmt::dict_add(d, pmt::mp("epoch_sample"), pmt::from_double((double)t0_abs_ + eng_->epoch_sample()));
            d = pmt::dict_add(d, pmt::mp("lock"), pmt::from_double(eng_->lock));
            d = pmt::dict_add(d, pmt::mp("sample_abs"), pmt::from_long((long)(t0_abs_ + eng_->samples_in)));
            d = pmt::dict_add(d, pmt::mp("slot"), pmt::from_long(slot_));
            message_port_pub(pmt::mp("obs"), d);
        }
        if (lost_run_ >= LOST_PERIODS) {
            pmt::pmt_t ex = pmt::make_dict();
            ex = pmt::dict_add(ex, pmt::mp("periods"), pmt::from_long((long)n_periods_));
            publish_status("lost", ex);
            eng_.reset();
            prn_ = 0;
            break;
        }
    }
    consume(0, consumed);
    if (produced % batch_) { // a short batch (end of input, or lost): pad it
        const int pad = batch_ - produced % batch_;
        std::memset(out + produced, 0, pad * sizeof(gr_complex));
        produced += pad;
    }
    return produced;
}

} // namespace gpsrx
} // namespace gr
