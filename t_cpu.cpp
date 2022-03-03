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
[[maybe_unused]] auto _pp_var_u0 = pp->globals[0];\
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
[[maybe_unused]] auto* _pp_var_tsyn = pp->state_vars[25];\
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
        auto vec_dii_ = _pp_var_vec_di[node_indexi_];
        arb_value_type t = _pp_var_vec_t[vec_dii_];
        _pp_var_C[i_] =  1.0;
        _pp_var_O[i_] =  0.;
        _pp_var_D[i_] =  0.;
        _pp_var_delay[i_] =  0.;
        _pp_var_T[i_] =  0.;
        _pp_var_Trelease[i_] =  0.;
        _pp_var_Mres[i_] =  1.6605778811026237e-06*_pp_var_M[i_];
        _pp_var_numpulses[i_] =  0.;
        _pp_var_xview[i_] =  1.0;
        _pp_var_yview[i_] =  0.;
        _pp_var_zview[i_] =  0.;
        _pp_var_Pview[i_] =  0.;
        _pp_var_on[i_] =  0.;
        _pp_var_y[i_] =  0.;
        _pp_var_z[i_] =  0.;
        _pp_var_u[i_] = _pp_var_u0;
        _pp_var_tsyn[i_] = t;
        _pp_var_nspike[i_] =  1.0;
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
        arb_value_type t_8_, t_7_, t_6_, t_5_, a_0_, a_3_, a_7_, t_1_, a_1_, a_4_, b_0_, t_2_, a_2_, a_5_, t_4_, a_6_, t_0_, t_3_;
        _pp_var_r1[i_] = _pp_var_r1FIX[i_]*pow(_pp_var_Trelease[i_],  2.0)/pow(_pp_var_Trelease[i_]+_pp_var_kB[i_],  2.0);
        _pp_var_r6[i_] = _pp_var_r6FIX[i_]*pow(_pp_var_Trelease[i_],  2.0)/pow(_pp_var_Trelease[i_]+_pp_var_kB[i_],  2.0);
        a_0_ =  1.0;
        a_1_ =  1.0;
        a_2_ =  1.0;
        a_3_ =  1.0;
        a_4_ =  -( -1.0* -_pp_var_r6[i_]*dt);
        a_5_ =  1.0- -1.0*_pp_var_r5[i_]*dt;
        a_6_ =  -(_pp_var_r1[i_]*dt);
        a_7_ =  1.0- -_pp_var_r2[i_]*dt;
        t_0_ = a_7_*a_0_-a_2_*a_6_;
        t_1_ = a_7_*a_1_;
        t_2_ = a_7_*a_3_-a_2_*_pp_var_O[i_];
        t_3_ = a_5_*t_0_-t_1_*a_4_;
        t_4_ = a_5_*t_2_-t_1_*_pp_var_D[i_];
        t_5_ = t_3_*a_5_;
        t_6_ = t_3_*_pp_var_D[i_]-a_4_*t_4_;
        t_7_ = t_3_*a_7_;
        t_8_ = t_3_*_pp_var_O[i_]-a_6_*t_4_;
        _pp_var_C[i_] = t_4_/t_3_;
        _pp_var_D[i_] = t_6_/t_5_;
        _pp_var_O[i_] = t_8_/t_7_;
        b_0_ =  -1.0;
        _pp_var_delay[i_] = _pp_var_delay[i_]+b_0_*dt;
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        auto vec_dii_ = _pp_var_vec_di[node_indexi_];
        arb_value_type conductivity_ = 0;
        arb_value_type current_ = 0;
        arb_value_type t = _pp_var_vec_t[vec_dii_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type i = 0;
        _pp_var_g[i_] = _pp_var_gmax[i_]*_pp_var_O[i_];
        i =  9.9999999999999995e-07*_pp_var_g[i_]*(v-_pp_var_Erev[i_]);
        if (_pp_var_delay[i_]== 0.) {
            _pp_var_t0[i_] = t;
            _pp_var_T[i_] =  0.;
            _pp_var_on[i_] =  0.;
        }
        current_ = i;
        conductivity_ =  9.9999999999999995e-07*_pp_var_g[i_];
        _pp_var_vec_g[node_indexi_] = fma(_pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
        _pp_var_vec_i[node_indexi_] = fma(_pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
    }
}

static void write_ions(arb_mechanism_ppack* pp) {
}

static void apply_events(arb_mechanism_ppack* pp, arb_deliverable_event_stream* stream_ptr) {
    PPACK_IFACE_BLOCK;
    auto ncell = stream_ptr->n_streams;
    for (arb_size_type c = 0; c<ncell; ++c) {
        auto begin  = stream_ptr->events + stream_ptr->begin[c];
        auto end    = stream_ptr->events + stream_ptr->end[c];
        for (auto p = begin; p<end; ++p) {
            auto i_     = p->mech_index;
            auto weight = p->weight;
            if (p->mech_id==_pp_var_mechanism_id) {
                auto vec_dii_ = _pp_var_vec_di[node_indexi_];
                arb_value_type t = _pp_var_vec_t[vec_dii_];
                _pp_var_nspike[i_] = _pp_var_nspike[i_]+ 1.0;
                if (_pp_var_on[i_]!= 1.0) {
                    _pp_var_t0[i_] = t;
                    _pp_var_on[i_] =  1.0;
                    _pp_var_z[i_] = _pp_var_z[i_]*exp( -(t-_pp_var_tsyn[i_])/_pp_var_tau_rec[i_]);
                    _pp_var_z[i_] = _pp_var_z[i_]+_pp_var_y[i_]*(exp( -(t-_pp_var_tsyn[i_])/_pp_var_tau_1[i_])-exp( -(t-_pp_var_tsyn[i_])/_pp_var_tau_rec[i_]))/(_pp_var_tau_1[i_]/_pp_var_tau_rec[i_]- 1.0);
                    _pp_var_y[i_] = _pp_var_y[i_]*exp( -(t-_pp_var_tsyn[i_])/_pp_var_tau_1[i_]);
                    _pp_var_x[i_] =  1.0-_pp_var_y[i_]-_pp_var_z[i_];
                    if (_pp_var_tau_facil[i_]> 0.) {
                        _pp_var_u[i_] = _pp_var_u[i_]*exp( -(t-_pp_var_tsyn[i_])/_pp_var_tau_facil[i_]);
                        _pp_var_u[i_] = _pp_var_u[i_]+_pp_var_U[i_]*( 1.0-_pp_var_u[i_]);
                    }
                    else {
                        _pp_var_u[i_] = _pp_var_U[i_];
                    }
                    _pp_var_y[i_] = _pp_var_y[i_]+_pp_var_x[i_]*_pp_var_u[i_];
                    _pp_var_xview[i_] = _pp_var_x[i_];
                    _pp_var_yview[i_] = _pp_var_y[i_];
                    _pp_var_Pview[i_] = _pp_var_u[i_];
                    _pp_var_T[i_] = _pp_var_Tmax[i_]*_pp_var_y[i_];
                    _pp_var_numpulses[i_] = _pp_var_numpulses[i_]+ 1.0;
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
  arb_mechanism_interface* make__Ampa_interface_multicore() {
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

