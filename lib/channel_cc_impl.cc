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
#include "galileo_e1_codes.h"
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

inline double tracker::chip(const std::vector<int8_t>& c, double phase_chips) const
{
    // nearest-chip lookup; with BOC(1,1) the chip is +1 in its first half and -1 in its second
    const double fl = std::floor(phase_chips);
    int i = (int)fl % code_len;
    if (i < 0)
        i += code_len;
    double v = (double)c[i];
    if (boc_ && (phase_chips - fl) >= 0.5)
        v = -v;
    return v;
}

tracker::tracker(int prn_, double fs_, double doppler_hz, double code_phase_samples, double pll_bw, double dll_bw,
                 double pll_bw_narrow, double dll_bw_narrow, int coherent_ms, int pll_order, const signal* sig)
    : prn(prn_), fs(fs_), carrier_hz(doppler_hz), doppler0_(doppler_hz),
      pll_bw_narrow_(pll_bw_narrow), dll_bw_narrow_(dll_bw_narrow), pll_order_(pll_order)
{
    if (sig && !sig->code.empty()) {
        code_ = sig->code;
        data_code_ = sig->data_code;
        secondary_ = sig->secondary;
        boc_ = sig->boc;
        bit_periods_ = sig->bit_periods;
        spacing_ = sig->spacing;
        code_len = (int)code_.size();
        sys = "GAL";
    } else {
        std::array<int8_t, CODE_LEN> ca;
        ca_code(prn, ca);
        code_.assign(ca.begin(), ca.end());
        code_len = CODE_LEN;
    }
    // coherent_ms is MILLISECONDS; the window is counted in periods and must divide the bit /
    // secondary-code length (track.py)
    const double period_ms = 1e3 * code_len / CODE_RATE;
    int cohp = std::max(1, std::min((int)std::lround(coherent_ms / period_ms), bit_periods_));
    while (bit_periods_ % cohp)
        cohp--;
    coh_ = (pll_bw_narrow > 0 && bit_periods_ > 1) ? cohp : 1;
    // acquisition hands over the sample at which the code starts; the phase at sample 0 is
    // therefore -(that) chips/sample, modulo the code length
    const double cps = CODE_RATE / fs;
    code_phase = std::fmod(-code_phase_samples * cps, (double)code_len);
    if (code_phase < 0)
        code_phase += code_len;
    loop_gains(pll_bw, 0.25, pll_t1_, pll_t2_); // atan/2pi is +-1/4 cycle full scale
    loop_gains(dll_bw, 1.0, dll_t1_, dll_t2_);
    n_nominal_ = (int)std::lround(fs * code_len / CODE_RATE);
    sc_target_ = (int64_t)std::lround(60.0 / period_s());
}

int tracker::samples_needed() const
{
    const double chips_left = code_len - code_phase;
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
    // 2. early / prompt / late replicas `spacing` chips apart, nearest-chip lookup, accumulated;
    //    a pilot-aided signal adds a fourth correlator on the data code
    const double cps = code_rate / fs;
    const bool pilot = !data_code_.empty();
    std::complex<double> E = 0, P = 0, L = 0, D = 0;
    double phc = code_phase;
    for (int k = 0; k < n; k++) {
        const std::complex<double> xb = std::complex<double>(x[k].real(), x[k].imag()) * ph;
        ph *= stepc;
        E += xb * chip(code_, phc + spacing_);
        P += xb * chip(code_, phc);
        L += xb * chip(code_, phc - spacing_);
        if (pilot)
            D += xb * chip(data_code_, phc);
        phc += cps;
    }
    const std::complex<double> Pdata = pilot ? D : P;
    // advance phases by what this period consumed
    const double dt = n / fs;
    carrier_phase = std::fmod(carrier_phase + 2.0 * M_PI * carrier_hz * dt, 2.0 * M_PI);
    if (carrier_phase < 0)
        carrier_phase += 2.0 * M_PI;
    carrier_cycles += carrier_hz * dt;
    code_phase = std::fmod(code_phase + n * cps, (double)code_len);
    samples_in += n;
    epochs += 1;
    double ip = P.real(), qp = P.imag();
    // lock indicator (I^2 - Q^2) / (I^2 + Q^2), smoothed over ~50 periods
    const double p = ip * ip + qp * qp;
    const double li = p > 0 ? (ip * ip - qp * qp) / p : 0.0;
    lock = 0.98 * lock + 0.02 * li;
    // C/N0, moment method over 20 prompts (track.py)
    m2_ += p;
    m4_ += p * p;
    if (++mn_ == 20) {
        const double m2 = m2_ / 20, m4 = m4_ / 20;
        const double pd = std::sqrt(std::max(2 * m2 * m2 - m4, 0.0)), pn = m2 - pd;
        if (pd > 0 && pn > 0) {
            const double cn0 = 10 * std::log10(pd / pn / period_s());
            cn0_db = cn0_db == 0.0 ? cn0 : 0.9 * cn0_db + 0.1 * cn0;
        }
        m2_ = m4_ = 0.0;
        mn_ = 0;
    }
    // stage 1: bit sync (where, mod 20, does the prompt's sign flip? - or, for a pilot with a known
    // secondary code, where in the sequence are we?) and the averaged loop state for the handover
    if (bit_offset < 0 && coh_ > 1) {
        f_avg_ = have_f_avg_ ? 0.99 * f_avg_ + 0.01 * carrier_hz : carrier_hz;
        have_f_avg_ = true;
        cd_avg_ = 0.99 * cd_avg_ + 0.01 * code_dop_;
        if (!secondary_.empty())
            secondary_sync(ip);
        else
            bit_sync(ip);
    } else if (bit_offset >= 0) {
        monitor_grid(ip);
    }
    // stage 2: coherent window over a whole bit, one loop update per window (track.py)
    std::complex<double> Ew = E, Pw = P, Lw = L;
    double dt_loop = dt;
    if (bit_offset >= 0) {
        if (!secondary_.empty()) { // wipe the known secondary chip off this period's pilot correlations
            const int m = (int)secondary_.size();
            const double c = sec_pol_ * (double)secondary_[(int)(((epochs - 1 - bit_offset) % m + m) % m)];
            E *= c; P *= c; L *= c;
        }
        const bool at_edge = ((epochs - bit_offset) % coh_ + coh_) % coh_ == 0;
        if (!aligned_) { // no loop update until the first bit edge (a partial window kicks the loop)
            aligned_ = at_edge;
            return Pdata;
        }
        acc_[0] += E; acc_[1] += P; acc_[2] += L;
        acc_dt_ += dt;
        if (!at_edge)
            return Pdata;
        Ew = acc_[0]; Pw = acc_[1]; Lw = acc_[2];
        dt_loop = acc_dt_;
        acc_[0] = acc_[1] = acc_[2] = 0;
        acc_dt_ = 0.0;
    }
    const double ipw = Pw.real(), qpw = Pw.imag();
    // 3. PLL: Costas discriminator atan(Q/I) (NOT atan2: blind to the data flips), error in
    //    cycles, loop filter -> frequency correction; the NCO phase accumulation integrates it
    // a pilot with its secondary code wiped carries no data: a pure four-quadrant PLL, +-1/2 cycle
    // of pull-in and no half-cycle ambiguity (track.py)
    const bool pure_pll = !secondary_.empty() && bit_offset >= 0;
    const double e_pll = pure_pll ? std::atan2(qpw, ipw) / (2.0 * M_PI)
                                  : (ipw != 0.0 ? std::atan(qpw / ipw) / (2.0 * M_PI) : 0.0);
    if (w3_ > 0.0) { // third order: acceleration integrator, rate integrator, direct term (track.py)
        acc3_ += w3_ * w3_ * w3_ * e_pll * dt_loop;
        vel3_ += (acc3_ + 1.1 * w3_ * w3_ * e_pll) * dt_loop;
        carr_corr_ = vel3_ + 2.4 * w3_ * e_pll;
    } else {
        carr_corr_ += (pll_t2_ * (e_pll - pll_e_prev) + dt_loop * e_pll) / pll_t1_;
    }
    pll_e_prev = e_pll;
    if (bit_offset >= 0) {
        // scintillation per 60 s block from the coherent-window power and the residual phase error,
        // the thermal-noise part of S4 taken out from the measured C/N0 (track.py)
        const double pw = ipw * ipw + qpw * qpw;
        sc_p1_ += pw;
        sc_p2_ += pw * pw;
        sc_n_ += coh_;
        sc_e1_ += e_pll * 2.0 * M_PI;
        sc_e2_ += (e_pll * 2.0 * M_PI) * (e_pll * 2.0 * M_PI);
        sc_ne_++;
        if (sc_n_ >= sc_target_) {
            const double nw = (double)sc_ne_;
            const double m = sc_p1_ / nw;
            const double v = std::max(sc_p2_ / nw - m * m, 0.0);
            const double cn0_lin = cn0_db > 0 ? std::pow(10.0, cn0_db / 10.0) : 1.0;
            const double tw = coh_ * period_s();
            const double s4n2 = (1.0 / (cn0_lin * tw)) * (1.0 + 1.0 / (2.0 * cn0_lin * tw));
            s4 = m > 0 ? std::sqrt(std::max(v / (m * m) - s4n2, 0.0)) : -1.0;
            const double me = sc_e1_ / nw;
            sigma_phi = std::sqrt(std::max(sc_e2_ / nw - me * me, 0.0));
            sc_p1_ = sc_p2_ = sc_e1_ = sc_e2_ = 0.0;
            sc_n_ = sc_ne_ = 0;
        }
    }
    carrier_hz = doppler0_ + carr_corr_;
    // 4. DLL: normalised early-minus-late envelope, carrier-aided code rate
    const double Em = std::abs(Ew), Lm = std::abs(Lw);
    const double e_dll = (Em + Lm) > 0 ? 0.5 * (Em - Lm) / (Em + Lm) : 0.0;
    code_dop_ += (dll_t2_ * (e_dll - dll_e_prev) + dt_loop * e_dll) / dll_t1_;
    dll_e_prev = e_dll;
    code_rate = CODE_RATE + carrier_hz * CARRIER_TO_CODE + code_dop_;
    return Pdata;
}

void tracker::switch_to_narrow()
{
    // bandwidth x window capped at 0.3 (15 Hz x 20 ms; a 100 ms pilot window gets 3 Hz): at 1.5
    // the loop lost every satellite within a second (track.py)
    const double bw = std::min(pll_bw_narrow_, 0.3 / (coh_ * period_s()));
    // k = 1 for the narrow stage (see track.py: with k = 0.25 the 20 ms loop is over unity gain)
    loop_gains(bw, 1.0, pll_t1_, pll_t2_);
    loop_gains(dll_bw_narrow_, 1.0, dll_t1_, dll_t2_);
    pll_e_prev = dll_e_prev = 0.0;
    acc_[0] = acc_[1] = acc_[2] = 0;
    acc_dt_ = 0.0;
    aligned_ = false;
    if (have_f_avg_) { // start from the averaged frequency, not the jittering instantaneous one
        carr_corr_ = f_avg_ - doppler0_;
        carrier_hz = f_avg_;
    }
    if (pll_order_ == 3) {
        w3_ = bw / 0.7845;
        acc3_ = 0.0;
        vel3_ = have_f_avg_ ? f_avg_ - doppler0_ : carr_corr_;
    }
    code_dop_ = cd_avg_; // and the averaged code-rate correction (track.py)
    code_rate = CODE_RATE + carrier_hz * CARRIER_TO_CODE + code_dop_;
}

void tracker::secondary_sync(double ip)
{
    // the prompt's signs over the last 100 periods against the known secondary code at every
    // cyclic offset and both polarities; 90% agreement names the offset (track.py)
    sec_hist_.push_back(ip > 0 ? 1 : (ip < 0 ? -1 : 0));
    const int m = (int)secondary_.size();
    const int N = 4 * m;
    if ((int)sec_hist_.size() > N)
        sec_hist_.erase(sec_hist_.begin(), sec_hist_.begin() + (sec_hist_.size() - N));
    if ((int)sec_hist_.size() < N || epochs < 300 || lock <= 0.5)
        return;
    double best = 0.0;
    int best_off = 0;
    for (int off = 0; off < m; off++) {
        double sc = 0.0;
        for (int j = 0; j < N; j++) {
            const int64_t k = epochs - N + 1 + j; // the epoch count this sign belongs to (last = this period = epochs)
            sc += sec_hist_[j] * secondary_[(int)(((k - 1 - off) % m + m) % m)];
        }
        sc /= N;
        if (std::fabs(sc) > std::fabs(best)) {
            best = sc;
            best_off = off;
        }
    }
    if (std::fabs(best) >= 0.9) {
        bit_offset = best_off;
        sec_pol_ = best > 0 ? 1.0 : -1.0;
        switch_to_narrow();
    }
}

void tracker::bit_sync(double ip)
{
    // this period's index is epochs - 1 (epochs already counts it)
    if (last_ip_ != 0.0 && (ip < 0) != (last_ip_ < 0) && lock > 0.5)
        flips_[(epochs - 1) % bit_periods_]++;
    last_ip_ = ip;
    if (epochs < 300 || lock <= 0.5)
        return;
    int64_t top = 0, second = 0, arg = 0;
    for (int i = 0; i < bit_periods_; i++) {
        if (flips_[i] > top) { second = top; top = flips_[i]; arg = i; }
        else if (flips_[i] > second) second = flips_[i];
    }
    if (top >= 8 && top >= 3 * std::max<int64_t>(second, 1)) {
        bit_offset = (int)arg;
        switch_to_narrow();
    }
}

void tracker::monitor_grid(double ip)
{
    // stage 2, every period: is the bit grid (or the secondary-code alignment) still where the sync
    // put it? A whole number of code periods missing from the stream moves nothing the loops or the
    // solver can see (the code repeats; every count and the sample counter skip the same
    // millisecond; the decoders re-anchor on the shifted grid, consistently and 1 ms wrong for
    // ever - measured). The one thing that moves is where the data bits flip. (track.py)
    const int n = bit_periods_;
    if (!secondary_.empty()) {
        sec_hist_.push_back(ip > 0 ? 1 : (ip < 0 ? -1 : 0));
        const int N = 4 * n;
        if ((int)sec_hist_.size() > N)
            sec_hist_.erase(sec_hist_.begin(), sec_hist_.begin() + (sec_hist_.size() - N));
        if (++mon_count_ < N || (int)sec_hist_.size() < N)
            return;
        mon_count_ = 0;
        double best = 0.0, cur = 0.0;
        int best_off = 0;
        for (int off = 0; off < n; off++) {
            double sc = 0.0;
            for (int j = 0; j < N; j++) {
                const int64_t k = epochs - N + 1 + j;
                sc += sec_hist_[j] * secondary_[(int)(((k - 1 - off) % n + n) % n)];
            }
            sc /= N;
            if (off == bit_offset)
                cur = sc;
            if (std::fabs(sc) > std::fabs(best)) {
                best = sc;
                best_off = off;
            }
        }
        if (best_off != bit_offset && std::fabs(best) >= 0.9 && std::fabs(cur) < 0.5) {
            bit_offset = best_off;
            sec_pol_ = best > 0 ? 1.0 : -1.0;
            acc_[0] = acc_[1] = acc_[2] = 0;
            acc_dt_ = 0.0;
            aligned_ = false;
            slips++;
        }
        return;
    }
    if (last_ip_ != 0.0 && (ip < 0) != (last_ip_ < 0))
        flips2_[(epochs - 1) % n]++;
    last_ip_ = ip;
    if (++mon_count_ < 300)
        return;
    mon_count_ = 0;
    int top = 0;
    for (int i = 1; i < n; i++)
        if (flips2_[i] > flips2_[top])
            top = i;
    if (top != bit_offset && flips2_[top] >= 8 && flips2_[top] >= 3 * std::max<int64_t>(flips2_[bit_offset], 1)) {
        bit_offset = top;
        acc_[0] = acc_[1] = acc_[2] = 0;
        acc_dt_ = 0.0;
        aligned_ = false;
        slips++;
    }
    for (int i = 0; i < n; i++)
        flips2_[i] = 0;
}

// ---- the block ------------------------------------------------------------------------------
channel_cc::sptr channel_cc::make(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                                  double pll_bw_narrow, double dll_bw_narrow, int coherent_ms, int pll_order,
                                  const std::string& signal)
{
    return gnuradio::make_block_sptr<channel_cc_impl>(samp_rate, slot, pll_bw, dll_bw, obs_every_ms,
                                                      pll_bw_narrow, dll_bw_narrow, coherent_ms, pll_order, signal);
}

channel_cc_impl::channel_cc_impl(double samp_rate, int slot, double pll_bw, double dll_bw, int obs_every_ms,
                                 double pll_bw_narrow, double dll_bw_narrow, int coherent_ms, int pll_order,
                                 const std::string& signal)
    : gr::block("gpsrx_channel_cc",
                gr::io_signature::make(1, 1, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      fs_(samp_rate),
      slot_(slot),
      signal_(signal),
      pll_bw_(pll_bw),
      dll_bw_(dll_bw),
      pll_bw_narrow_(pll_bw_narrow),
      dll_bw_narrow_(dll_bw_narrow),
      coherent_ms_(coherent_ms),
      pll_order_(pll_order),
      obs_every_(obs_every_ms),
      batch_(3)
{
    // obs_every is milliseconds; a Galileo period is 4 ms
    obs_periods_ = std::max(1, signal_ == "L1CA" ? obs_every_ : obs_every_ / 4);
    message_port_register_in(pmt::mp("assign"));
    set_msg_handler(pmt::mp("assign"), [this](pmt::pmt_t msg) { this->on_assign(msg); });
    message_port_register_out(pmt::mp("obs"));
    message_port_register_out(pmt::mp("status"));
    // the same batching as the Python block: whole batches or nothing (see channel.py)
    set_output_multiple(batch_);
    set_min_output_buffer(16 * batch_);
}

channel_cc_impl::~channel_cc_impl() {}

// the E1 primary code from its hex string (galileo_e1_codes.h): 4 chips per character, bit 1 = -1
static std::vector<int8_t> e1_code(const char* hex)
{
    std::vector<int8_t> out;
    out.reserve(4092);
    for (const char* c = hex; *c; c++) {
        const int v = (*c >= 'A') ? *c - 'A' + 10 : *c - '0';
        for (int b = 3; b >= 0; b--)
            out.push_back((v >> b) & 1 ? -1 : 1);
    }
    return out;
}

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
    start_at_ = (int64_t)num(msg, "start_sample", 0.0); // replay: start exactly here (deterministic)
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
    d = pmt::dict_add(d, pmt::mp("system"), pmt::string_to_symbol(signal_ == "L1CA" ? "GPS" : "GAL"));
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
        if (start_at_ > start) { // not yet: eat input up to the start sample, no more
            consume(0, (int)std::min<int64_t>(n_in, start_at_ - start));
            return 0;
        }
        // acquisition's code start is `sample` absolute, seconds old by now: wrap it forward
        // with the DOPPLER-SHIFTED code period (nominal would be 3 chips/s wrong at 5 kHz)
        const double period = fs_ * (signal_ == "L1CA" ? tracker::CODE_LEN : 4092) / (tracker::CODE_RATE * (1.0 + pend_dop_ / tracker::L1_HZ));
        double rel = std::fmod(pend_sample_ - (double)start, period);
        if (rel < 0)
            rel += period;
        const tracker::signal* sigp = nullptr;
        if (signal_ != "L1CA") {
            const int i = pend_prn_ - 1;
            if (i < 0 || i >= 36) {
                publish_status("idle");   // no such Galileo PRN
                prn_ = 0;
                have_pending_ = false;
                consume(0, n_in);
                return 0;
            }
            sig_ = tracker::signal();
            sig_.name = signal_;
            sig_.boc = true;
            sig_.spacing = 0.25;
            if (signal_ == "E1") {         // pilot-aided: loops on E1-C, data from E1-B
                sig_.code = e1_code(GALILEO_E1C_HEX[i]);
                sig_.data_code = e1_code(GALILEO_E1B_HEX[i]);
                static const char* CS25 = "0011100000001010110110010";
                for (const char* c = CS25; *c; c++)
                    sig_.secondary.push_back(*c == '1' ? -1 : 1);
                sig_.bit_periods = 25;
            } else {                         // E1B: data channel alone, one symbol per period
                sig_.code = e1_code(GALILEO_E1B_HEX[i]);
                sig_.bit_periods = 1;
            }
            sigp = &sig_;
        }
        eng_.reset(new tracker(pend_prn_, fs_, pend_dop_, rel, pll_bw_, dll_bw_, pll_bw_narrow_, dll_bw_narrow_, coherent_ms_, pll_order_, sigp));
        t0_abs_ = start;
        lost_run_ = 0;
        slips_seen_ = 0;
        have_pending_ = false;
        pmt::pmt_t tag = pmt::make_dict();
        tag = pmt::dict_add(tag, pmt::mp("prn"), pmt::from_long(pend_prn_));
        tag = pmt::dict_add(tag, pmt::mp("slot"), pmt::from_long(slot_));
        tag = pmt::dict_add(tag, pmt::mp("system"), pmt::string_to_symbol(signal_ == "L1CA" ? "GPS" : "GAL"));
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
        if (eng_->slips != slips_seen_) { // the bit grid moved: whole code periods missing from the stream
            slips_seen_ = eng_->slips;
            pmt::pmt_t ex = pmt::make_dict();
            ex = pmt::dict_add(ex, pmt::mp("periods"), pmt::from_long((long)n_periods_));
            ex = pmt::dict_add(ex, pmt::mp("bit_offset"), pmt::from_long(eng_->bit_offset));
            publish_status("slip", ex);
            pmt::pmt_t tag = pmt::make_dict(); // the Nav Decoder downstream finds its bit grid again
            tag = pmt::dict_add(tag, pmt::mp("prn"), pmt::from_long(eng_->prn));
            tag = pmt::dict_add(tag, pmt::mp("slot"), pmt::from_long(slot_));
            add_item_tag(0, nitems_written(0) + produced - 1, pmt::mp("gpsrx_slip"), tag);
        }
        if (eng_->lock < LOST_LOCK)
            lost_run_++;
        else
            lost_run_ = 0;
        if (n_periods_ % obs_periods_ == 0) {
            pmt::pmt_t d = pmt::make_dict();
            d = pmt::dict_add(d, pmt::mp("prn"), pmt::from_long(eng_->prn));
            d = pmt::dict_add(d, pmt::mp("epochs"), pmt::from_long((long)eng_->epochs));
            d = pmt::dict_add(d, pmt::mp("code_phase"), pmt::from_double(eng_->code_phase));
            d = pmt::dict_add(d, pmt::mp("carrier_hz"), pmt::from_double(eng_->carrier_hz));
            d = pmt::dict_add(d, pmt::mp("samples_in"), pmt::from_long((long)eng_->samples_in));
            d = pmt::dict_add(d, pmt::mp("epoch_sample"), pmt::from_double((double)t0_abs_ + eng_->epoch_sample()));
            d = pmt::dict_add(d, pmt::mp("lock"), pmt::from_double(eng_->lock));
            d = pmt::dict_add(d, pmt::mp("cn0_db"), pmt::from_double(eng_->cn0_db));
            d = pmt::dict_add(d, pmt::mp("carrier_cycles"), pmt::from_double(eng_->epoch_cycles()));
            d = pmt::dict_add(d, pmt::mp("sample_abs"), pmt::from_long((long)(t0_abs_ + eng_->samples_in)));
            d = pmt::dict_add(d, pmt::mp("slot"), pmt::from_long(slot_));
            d = pmt::dict_add(d, pmt::mp("period_s"), pmt::from_double(eng_->period_s()));
            d = pmt::dict_add(d, pmt::mp("s4"), eng_->s4 >= 0 ? pmt::from_double(eng_->s4) : pmt::PMT_NIL);
            d = pmt::dict_add(d, pmt::mp("sigma_phi"), eng_->sigma_phi >= 0 ? pmt::from_double(eng_->sigma_phi) : pmt::PMT_NIL);
            d = pmt::dict_add(d, pmt::mp("sys"), pmt::string_to_symbol(eng_->sys));
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
