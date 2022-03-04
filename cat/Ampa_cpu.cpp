#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <arbor/mechanism_abi.h>
#include <arbor/math.hpp>

namespace kernel_Ampa {

using ::arb::math::exprelr;
using ::arb::math::safeinv;
using ::std::abs;
using ::std::cos;
using ::std::exp;
using ::std::log;
using ::std::max;
using ::std::min;
using ::std::pow;
using ::std::sin;

static constexpr unsigned simd_width_ = 1;
static constexpr unsigned min_align_ = std::max(alignof(arb_value_type), alignof(arb_index_type));

#define PPACK_IFACE_BLOCK \
[[maybe_unused]] auto  _pp_var_width             = pp->width;\
[[maybe_unused]] auto  _pp_var_n_detectors       = pp->n_detectors;\
[[maybe_unused]] auto* _pp_var_vec_ci            = pp->vec_ci;\
[[maybe_unused]] auto* _pp_var_vec_di            = pp->vec_di;\
[[maybe_unused]] auto* _pp_var_vec_t             = pp->vec_t;\
[[maybe_unused]] auto* _pp_var_vec_dt            = pp->vec_dt;\
[[maybe_unused]] auto* _pp_var_vec_v             = pp->vec_v;\
[[maybe_unused]] auto* _pp_var_vec_i             = pp->vec_i;\
[[maybe_unused]] auto* _pp_var_vec_g             = pp->vec_g;\
[[maybe_unused]] auto* _pp_var_temperature_degC  = pp->temperature_degC;\
[[maybe_unused]] auto* _pp_var_diam_um           = pp->diam_um;\
[[maybe_unused]] auto* _pp_var_time_since_spike  = pp->time_since_spike;\
[[maybe_unused]] auto* _pp_var_node_index        = pp->node_index;\
[[maybe_unused]] auto* _pp_var_peer_index        = pp->peer_index;\
[[maybe_unused]] auto* _pp_var_multiplicity      = pp->multiplicity;\
[[maybe_unused]] auto* _pp_var_weight            = pp->weight;\
[[maybe_unused]] auto& _pp_var_events            = pp->events;\
[[maybe_unused]] auto& _pp_var_mechanism_id      = pp->mechanism_id;\
[[maybe_unused]] auto& _pp_var_index_constraints = pp->index_constraints;\
\
[[maybe_unused]] auto _pp_var_u0 = pp->globals[0];\
\
[[maybe_unused]] auto* _pp_var_C = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_O = pp->state_vars[1];\
[[maybe_unused]] auto* _pp_var_D = pp->state_vars[2];\
[[maybe_unused]] auto* _pp_var_delay = pp->state_vars[3];\
[[maybe_unused]] auto* _pp_var_v = pp->state_vars[4];\
[[maybe_unused]] auto* _pp_var_g = pp->state_vars[5];\
[[maybe_unused]] auto* _pp_var_T = pp->state_vars[6];\
[[maybe_unused]] auto* _pp_var_r1 = pp->state_vars[7];\
[[maybe_unused]] auto* _pp_var_r6 = pp->state_vars[8];\
[[maybe_unused]] auto* _pp_var_Trelease = pp->state_vars[9];\
[[maybe_unused]] auto* _pp_var_x = pp->state_vars[10];\
[[maybe_unused]] auto* _pp_var_tsyn = pp->state_vars[11];\
[[maybe_unused]] auto* _pp_var_Mres = pp->state_vars[12];\
[[maybe_unused]] auto* _pp_var_NTdiffusion = pp->state_vars[13];\
[[maybe_unused]] auto* _pp_var_numpulses = pp->state_vars[14];\
[[maybe_unused]] auto* _pp_var_xview = pp->state_vars[15];\
[[maybe_unused]] auto* _pp_var_yview = pp->state_vars[16];\
[[maybe_unused]] auto* _pp_var_zview = pp->state_vars[17];\
[[maybe_unused]] auto* _pp_var_Pview = pp->state_vars[18];\
[[maybe_unused]] auto* _pp_var_on = pp->state_vars[19];\
[[maybe_unused]] auto* _pp_var_nspike = pp->state_vars[20];\
[[maybe_unused]] auto* _pp_var_t0 = pp->state_vars[21];\
[[maybe_unused]] auto* _pp_var_y = pp->state_vars[22];\
[[maybe_unused]] auto* _pp_var_z = pp->state_vars[23];\
[[maybe_unused]] auto* _pp_var_u = pp->state_vars[24];\
[[maybe_unused]] auto* _pp_var_tspike = pp->state_vars[25];\
[[maybe_unused]] auto* _pp_var_PRE    = pp->state_vars[75];\
\
[[maybe_unused]] auto* _pp_var_gmax = pp->parameters[0];\
[[maybe_unused]] auto* _pp_var_Cdur = pp->parameters[1];\
[[maybe_unused]] auto* _pp_var_Erev = pp->parameters[2];\
[[maybe_unused]] auto* _pp_var_kB = pp->parameters[3];\
[[maybe_unused]] auto* _pp_var_r1FIX = pp->parameters[4];\
[[maybe_unused]] auto* _pp_var_r6FIX = pp->parameters[5];\
[[maybe_unused]] auto* _pp_var_r2 = pp->parameters[6];\
[[maybe_unused]] auto* _pp_var_r5 = pp->parameters[7];\
[[maybe_unused]] auto* _pp_var_tau_1 = pp->parameters[8];\
[[maybe_unused]] auto* _pp_var_tau_rec = pp->parameters[9];\
[[maybe_unused]] auto* _pp_var_tau_facil = pp->parameters[10];\
[[maybe_unused]] auto* _pp_var_U = pp->parameters[11];\
[[maybe_unused]] auto* _pp_var_Tmax = pp->parameters[12];\
[[maybe_unused]] auto* _pp_var_M = pp->parameters[13];\
[[maybe_unused]] auto* _pp_var_R = pp->parameters[14];\
[[maybe_unused]] auto* _pp_var_Diff = pp->parameters[15];\
[[maybe_unused]] auto* _pp_var_lamd = pp->parameters[16];\
//End of IFACEBLOCK

// procedure prototypes

// interface methods
static void init(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        auto vec_dii_     = _pp_var_vec_di[node_indexi_];
        auto t            = _pp_var_vec_t[vec_dii_];
        _pp_var_C[i_] =  1.0;
        _pp_var_O[i_] =  0.;
        _pp_var_D[i_] =  0.;
        _pp_var_delay[i_] =  0.;
        _pp_var_T[i_] =  0.;
        _pp_var_Trelease[i_] =  0.;
        _pp_var_Mres[i_] =  1.6605778811026237e-06*_pp_var_M[i_];
        _pp_var_numpulses[i_] =  0.;
        _pp_var_zview[i_] =  0.;
        _pp_var_on[i_] =  0.;
        _pp_var_y[i_] =  0.;
        _pp_var_z[i_] =  0.;
        _pp_var_u[i_] = _pp_var_u0;
        _pp_var_tsyn[i_] = t;
        _pp_var_tspike[i_] =  1e12; // EDITTED
    }
    if (!_pp_var_multiplicity) return;
    for (arb_size_type ix = 0; ix < 4; ++ix) {
        for (arb_size_type iy = 0; iy < _pp_var_width; ++iy) {
            pp->state_vars[ix][iy] *= _pp_var_multiplicity[iy];
        }
    }
}

static void advance_state(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        
        // Read 
        auto O  = _pp_var_O[i_];
        auto D  = _pp_var_D[i_];
        auto r2 = _pp_var_r2[i_];
        auto r5 = _pp_var_r5[i_];

        auto tr = _pp_var_Trelease[i_];
        auto k  = _pp_var_kB[i_];
        auto ratio = std::pow(tr,2)/std::pow((tr+k),2);

        auto r1 = _pp_var_r1FIX[i_] * ratio;
        auto r6 = _pp_var_r6FIX[i_] * ratio;

        // Solve ODEs 
        auto t0  =  -r6 * dt;
        auto t1  =  -r1 * dt;
        auto t2  =  1.0 + r5*dt;
        auto t3  =  1.0 + r2*dt;
        auto t4  = -t3 * t1;
        auto t5  = t3 - O;
        auto t6  = (t2 * t4) - (t3 * t0);
        auto t7  = (t2 * t5) - (t3 * D);
        auto t8  = t6 * t2;
        auto t9  = (t6 * D) - (t0 * t7);
        auto t10 = t6 * t3;
        auto t11 = (t6 * O) - (t1 * t7);

        // Update 
        _pp_var_C[i_] = t7 / t6;
        _pp_var_D[i_] = t9 / t8;
        _pp_var_O[i_] = t11 / t10;
        _pp_var_delay[i_] -= dt;
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        auto vec_dii_ = _pp_var_vec_di[node_indexi_];

        arb_value_type t = _pp_var_vec_t[vec_dii_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];

        // START EDIT
        // Read
        const auto mres = _pp_var_Mres[i_];
        const auto r = _pp_var_R[i_];
        const auto diff = _pp_var_Diff[i_];
        const auto lamd = _pp_var_lamd[i_];
        const auto numpulses = (int)_pp_var_numpulses[i_];

        // Calculate
        const auto rsq = r*r;  
        const auto diff_4 = diff*4;
        const auto lamd_scaled = (1e-3)*lamd;

        auto NTdiffWave = _pp_var_T[i_];
        const auto max_pulses = std::min(numpulses, 50); 
	for (unsigned pulse = 0; pulse < max_pulses; ++pulse) {
            auto offset = pulse*_pp_var_width + i_;
            auto ts     = _pp_var_tspike[offset]; 

            auto delta_t = t - ts; 
            if (delta_t > 0.) {
                auto pre = _pp_var_PRE[offset]; 
                auto invariant = delta_t*diff_4; 
                NTdiffWave += pre*mres*std::exp(-rsq/invariant)/(3.14159*invariant*lamd_scaled);
            }
        }

        // Update
        _pp_var_Trelease[i_] = NTdiffWave;
        // END EDIT

        // Reset 
        if (_pp_var_delay[i_]== 0.) {
            _pp_var_T[i_] =  0.;
            _pp_var_on[i_] =  0.;
        }

        // Update
        const auto gmax = _pp_var_gmax[i_];
        const auto O    = _pp_var_O[i_];
        const auto Erev = _pp_var_Erev[i_];
        const auto weight = _pp_var_weight[i_];

        auto conductivity_ =  1e-06*gmax*O;
        auto current       =  conductivity_*(v-Erev);
        _pp_var_vec_g[node_indexi_] = fma(weight, conductivity_, _pp_var_vec_g[node_indexi_]);
        _pp_var_vec_i[node_indexi_] = fma(weight, current,       _pp_var_vec_i[node_indexi_]);
    }
}

static void write_ions(arb_mechanism_ppack* pp) {}

static void apply_events(arb_mechanism_ppack* pp, arb_deliverable_event_stream* stream_ptr) {
    PPACK_IFACE_BLOCK;
    auto ncell = stream_ptr->n_streams;
    for (arb_size_type c = 0; c<ncell; ++c) {

        auto begin  = stream_ptr->events + stream_ptr->begin[c];
        auto end    = stream_ptr->events + stream_ptr->end[c];

        for (auto p = begin; p<end; ++p) {
            auto i_     = p->mech_index;
            auto weight = p->weight;

            if (p->mech_id == _pp_var_mechanism_id) {
                auto node_indexi_ = _pp_var_node_index[i_]; // EDITTED - This line was missing.
                auto vec_dii_ = _pp_var_vec_di[node_indexi_];

                arb_value_type t = _pp_var_vec_t[vec_dii_];
                if (!_pp_var_on[i_]) {
                    _pp_var_on[i_] =  1.0;

                    // Read
                    auto z = _pp_var_z[i_];
                    auto y = _pp_var_y[i_];
                    auto x = _pp_var_x[i_];
                    auto u = _pp_var_u[i_];

                    const auto tsyn      = _pp_var_tsyn[i_];
                    const auto tau_rec   = _pp_var_tau_rec[i_];
                    const auto tau_1     = _pp_var_tau_1[i_];
                    const auto tau_facil = _pp_var_tau_facil[i_];
                    const auto U         = _pp_var_U[i_];
                    const auto Tmax      = _pp_var_Tmax[i_];
                    const auto numpulses = (int)_pp_var_numpulses[i_];

                    // Modify
                    z = z*exp(-(t-tsyn)/tau_rec);
                    z = z+y*(exp(-(t-tsyn)/tau_1) - exp(-(t-tsyn)/tau_rec)) / (tau_1/tau_rec - 1.0);
                    y = y*exp(-(t-tsyn)/tau_1);
                    x =  1.0-y-z;

                    if (tau_facil> 0.) {
                        u = u*exp(-(t-tsyn)/tau_facil);
                        u = u + U*(1.0-u);
                    }
                    else {
                        u = U;
                    }
                    y = y + x * u;

                    // Update
                    _pp_var_T[i_] = Tmax*y;
                    _pp_var_z[i_] = z;
                    _pp_var_y[i_] = y;
                    _pp_var_x[i_] = x;
                    _pp_var_u[i_] = u;

                    // START EDIT
                    auto offset = (numpulses%50)*_pp_var_width; // rolling window update
                    _pp_var_tspike[offset+i_] = t; 
                    _pp_var_PRE[offset+i_]    = y; 
                    // END EDIT

                    _pp_var_numpulses[i_] = numpulses + 1.;
                    _pp_var_tsyn[i_] = t;
                }
                _pp_var_delay[i_] = _pp_var_Cdur[i_];
            }
        }
    }
}

static void post_event(arb_mechanism_ppack*) {}

// Procedure definitions
#undef PPACK_IFACE_BLOCK
} // namespace kernel_Ampa

extern "C" {
  arb_mechanism_interface* make_arb_syn_catalogue_Ampa_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = kernel_Ampa::simd_width_;
    result.backend = arb_backend_kind_cpu;
    result.alignment = kernel_Ampa::min_align_;
    result.init_mechanism = kernel_Ampa::init;
    result.compute_currents = kernel_Ampa::compute_currents;
    result.apply_events = kernel_Ampa::apply_events;
    result.advance_state = kernel_Ampa::advance_state;
    result.write_ions = kernel_Ampa::write_ions;
    result.post_event = kernel_Ampa::post_event;
    return &result;
  }}

