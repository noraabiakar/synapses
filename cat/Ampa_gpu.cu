#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>
#include <arbor/gpu/reduce_by_key.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {
namespace Ampa_catalogue {

#define PPACK_IFACE_BLOCK \
auto  _pp_var_width             __attribute__((unused)) = params_.width;\
auto  _pp_var_n_detectors       __attribute__((unused)) = params_.n_detectors;\
auto* _pp_var_vec_ci            __attribute__((unused)) = params_.vec_ci;\
auto* _pp_var_vec_di            __attribute__((unused)) = params_.vec_di;\
auto* _pp_var_vec_t             __attribute__((unused)) = params_.vec_t;\
auto* _pp_var_vec_dt            __attribute__((unused)) = params_.vec_dt;\
auto* _pp_var_vec_v             __attribute__((unused)) = params_.vec_v;\
auto* _pp_var_vec_i             __attribute__((unused)) = params_.vec_i;\
auto* _pp_var_vec_g             __attribute__((unused)) = params_.vec_g;\
auto* _pp_var_temperature_degC  __attribute__((unused)) = params_.temperature_degC;\
auto* _pp_var_diam_um           __attribute__((unused)) = params_.diam_um;\
auto* _pp_var_time_since_spike  __attribute__((unused)) = params_.time_since_spike;\
auto* _pp_var_node_index        __attribute__((unused)) = params_.node_index;\
auto* _pp_var_peer_index        __attribute__((unused)) = params_.peer_index;\
auto* _pp_var_multiplicity      __attribute__((unused)) = params_.multiplicity;\
auto* _pp_var_state_vars        __attribute__((unused)) = params_.state_vars;\
auto* _pp_var_weight            __attribute__((unused)) = params_.weight;\
auto& _pp_var_events            __attribute__((unused)) = params_.events;\
auto& _pp_var_mechanism_id      __attribute__((unused)) = params_.mechanism_id;\
auto& _pp_var_index_constraints __attribute__((unused)) = params_.index_constraints;\
\
auto _pp_var_u0 __attribute__((unused)) = params_.globals[0];\
\
auto* _pp_var_C __attribute__((unused)) = params_.state_vars[0];\
auto* _pp_var_O __attribute__((unused)) = params_.state_vars[1];\
auto* _pp_var_D __attribute__((unused)) = params_.state_vars[2];\
auto* _pp_var_delay __attribute__((unused)) = params_.state_vars[3];\
auto* _pp_var_v __attribute__((unused)) = params_.state_vars[4];\
auto* _pp_var_g __attribute__((unused)) = params_.state_vars[5];\
auto* _pp_var_T __attribute__((unused)) = params_.state_vars[6];\
auto* _pp_var_Trelease __attribute__((unused)) = params_.state_vars[7];\
auto* _pp_var_x __attribute__((unused)) = params_.state_vars[8];\
auto* _pp_var_tsyn __attribute__((unused)) = params_.state_vars[9];\
auto* _pp_var_Mres __attribute__((unused)) = params_.state_vars[10];\
auto* _pp_var_NTdiffusion __attribute__((unused)) = params_.state_vars[11];\
auto* _pp_var_numpulses __attribute__((unused)) = params_.state_vars[12];\
auto* _pp_var_on __attribute__((unused)) = params_.state_vars[13];\
auto* _pp_var_t0 __attribute__((unused)) = params_.state_vars[14];\
auto* _pp_var_y __attribute__((unused)) = params_.state_vars[15];\
auto* _pp_var_z __attribute__((unused)) = params_.state_vars[16];\
auto* _pp_var_u __attribute__((unused)) = params_.state_vars[17];\
\
auto* _pp_var_gmax __attribute__((unused)) = params_.parameters[0];\
auto* _pp_var_Cdur __attribute__((unused)) = params_.parameters[1];\
auto* _pp_var_Erev __attribute__((unused)) = params_.parameters[2];\
auto* _pp_var_kB __attribute__((unused)) = params_.parameters[3];\
auto* _pp_var_r1FIX __attribute__((unused)) = params_.parameters[4];\
auto* _pp_var_r6FIX __attribute__((unused)) = params_.parameters[5];\
auto* _pp_var_r2 __attribute__((unused)) = params_.parameters[6];\
auto* _pp_var_r5 __attribute__((unused)) = params_.parameters[7];\
auto* _pp_var_tau_1 __attribute__((unused)) = params_.parameters[8];\
auto* _pp_var_tau_rec __attribute__((unused)) = params_.parameters[9];\
auto* _pp_var_tau_facil __attribute__((unused)) = params_.parameters[10];\
auto* _pp_var_U __attribute__((unused)) = params_.parameters[11];\
auto* _pp_var_Tmax __attribute__((unused)) = params_.parameters[12];\
auto* _pp_var_M __attribute__((unused)) = params_.parameters[13];\
auto* _pp_var_R __attribute__((unused)) = params_.parameters[14];\
auto* _pp_var_Diff __attribute__((unused)) = params_.parameters[15];\
auto* _pp_var_lamd __attribute__((unused)) = params_.parameters[16];\
\
arb_value_type* __attribute__((unused)) _pp_var_tspike[50] = {\
    pp->state_vars[18], pp->state_vars[19], pp->state_vars[20], pp->state_vars[21], pp->state_vars[22],\
    pp->state_vars[23], pp->state_vars[24], pp->state_vars[25], pp->state_vars[26], pp->state_vars[27],\
    pp->state_vars[28], pp->state_vars[29], pp->state_vars[30], pp->state_vars[31], pp->state_vars[32],\
    pp->state_vars[33], pp->state_vars[34], pp->state_vars[35], pp->state_vars[36], pp->state_vars[37],\
    pp->state_vars[38], pp->state_vars[39], pp->state_vars[40], pp->state_vars[41], pp->state_vars[42],\
    pp->state_vars[43], pp->state_vars[44], pp->state_vars[45], pp->state_vars[46], pp->state_vars[47],\
    pp->state_vars[48], pp->state_vars[49], pp->state_vars[50], pp->state_vars[51], pp->state_vars[52],\
    pp->state_vars[53], pp->state_vars[54], pp->state_vars[55], pp->state_vars[56], pp->state_vars[57],\
    pp->state_vars[58], pp->state_vars[59], pp->state_vars[60], pp->state_vars[61], pp->state_vars[62],\
    pp->state_vars[63], pp->state_vars[64], pp->state_vars[65], pp->state_vars[66], pp->state_vars[67],\
};\
\
arb_value_type* ___attribute__((unused)) pp_var_PRE[50] = {\
    pp->state_vars[68], pp->state_vars[69], pp->state_vars[70], pp->state_vars[71], pp->state_vars[72],\
    pp->state_vars[73], pp->state_vars[74], pp->state_vars[75], pp->state_vars[76], pp->state_vars[77],\
    pp->state_vars[78], pp->state_vars[79], pp->state_vars[80], pp->state_vars[81], pp->state_vars[82],\
    pp->state_vars[83], pp->state_vars[84], pp->state_vars[85], pp->state_vars[86], pp->state_vars[87],\
    pp->state_vars[88], pp->state_vars[89], pp->state_vars[90], pp->state_vars[91], pp->state_vars[92],\
    pp->state_vars[93], pp->state_vars[94], pp->state_vars[95], pp->state_vars[96], pp->state_vars[97],\
    pp->state_vars[98], pp->state_vars[99], pp->state_vars[100], pp->state_vars[101], pp->state_vars[102],\
    pp->state_vars[103], pp->state_vars[104], pp->state_vars[105], pp->state_vars[106], pp->state_vars[107],\
    pp->state_vars[108], pp->state_vars[109], pp->state_vars[110], pp->state_vars[111], pp->state_vars[112],\
    pp->state_vars[113], pp->state_vars[114], pp->state_vars[115], pp->state_vars[116], pp->state_vars[117],\
};\
//End of IFACEBLOCK

namespace {

using ::arb::gpu::exprelr;
using ::arb::gpu::safeinv;
using ::arb::gpu::min;
using ::arb::gpu::max;

__global__
void init(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto node_indexi_ = _pp_var_node_index[tid_];
        auto vec_dii_     = _pp_var_vec_di[node_indexi_];
        arb_value_type t = _pp_var_vec_t[vec_dii_];
        _pp_var_C[tid_] =  1.0;
        _pp_var_O[tid_] =  0.;
        _pp_var_D[tid_] =  0.;
        _pp_var_delay[tid_] =  0.;
        _pp_var_T[tid_] =  0.;
        _pp_var_Trelease[tid_] =  0.;
        _pp_var_Mres[tid_] =  1.6605778811026237e-06*_pp_var_M[tid_];
        _pp_var_numpulses[tid_] =  0.;
        _pp_var_on[tid_] =  0.;
        _pp_var_y[tid_] =  0.;
        _pp_var_z[tid_] =  0.;
        _pp_var_u[tid_] = _pp_var_u0;
        _pp_var_tsyn[tid_] = t;
        _pp_var_tspike[0][tid_] =  1e12; // EDITTED
    }
}

__global__
void multiply(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    auto tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    auto idx_ = blockIdx.y;    if(tid_<_pp_var_width) {
        _pp_var_state_vars[idx_][tid_] *= _pp_var_multiplicity[tid_];
    }
}

__global__
void advance_state(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto node_indexi_ = _pp_var_node_index[tid_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];

        // Read 
        double O  = _pp_var_O[tid_];
        double D  = _pp_var_D[tid_];
        double r2 = _pp_var_r2[tid_];
        double r5 = _pp_var_r5[tid_];

        double tr = _pp_var_Trelease[i_];
        double k  = _pp_var_kB[i_];
        double ratio = std::pow(tr,2)/std::pow((tr+k),2);

        double r1 = _pp_var_r1FIX[i_] * ratio;
        double r6 = _pp_var_r6FIX[i_] * ratio;

        // Solve ODEs 
        double t0  =  -r6 * dt;
        double t1  =  -r1 * dt;
        double t2  =  1.0 + r5*dt;
        double t3  =  1.0 + r2*dt;
        double t4  = t3 - t1;
        double t5  = t3 - O;
        double t6  = (t2 * t4) - (t3 * t0);
        double t7  = (t2 * t5) - (t3 * D);
        double t8  = t6 * t2;
        double t9  = (t6 * D) - (t0 * t7);
        double t10 = t6 * t3;


        // Update 
        if (tr > 0) {
            _pp_var_C[tid_] = t7 / t6;
            _pp_var_D[tid_] = t9 / t8;
            _pp_var_O[tid_] = t11 / t10;
        }
        _pp_var_delay[tid_] -= dt;
    }
}

__global__
void compute_currents(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned lane_mask_ = arb::gpu::ballot(0xffffffff, tid_<n_);
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto node_indexi_ = _pp_var_node_index[tid_];
        auto vec_dii_ = _pp_var_vec_di[node_indexi_];

        arb_value_type t = _pp_var_vec_t[vec_dii_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];

        // START EDIT
        // Read
        const auto mres = _pp_var_Mres[tid_];
        const auto r = _pp_var_R[tid_];
        const auto diff = _pp_var_Diff[tid_];
        const auto lamd = _pp_var_lamd[tid_];
        const auto numpulses = (int)_pp_var_numpulses[tid_];

        // Calculate
        const auto rsq = r*r;  
        const auto diff_4 = diff*4;
        const auto lamd_scaled = (1e-3)*lamd;

        auto NTdiffWave = _pp_var_T[tid_];
        const auto max_pulses = std::min(numpulses, 50); 
	for (unsigned pulse = 0; pulse < max_pulses; ++pulse) {
            auto ts     = _pp_var_tspike[pulse][tid_]; 

            auto delta_t = t - ts; 
            if (delta_t > 0.) {
                auto pre = _pp_var_PRE[pulse][tid_]; 
                auto invariant = delta_t*diff_4; 
                NTdiffWave += pre*mres*std::exp(-rsq/invariant)/(3.14159*invariant*lamd_scaled);
            }
        }

        // Update
        _pp_var_Trelease[tid_] = NTdiffWave;
        // END EDIT

        // Reset 
        if (_pp_var_delay[tid_]< 0.) {
            _pp_var_T[tid_] =  0.;
            _pp_var_on[tid_] =  0.;
        }

        // Update
        const auto gmax = _pp_var_gmax[tid_];
        const auto O    = _pp_var_O[tid_];
        const auto Erev = _pp_var_Erev[tid_];
        const auto weight = _pp_var_weight[tid_];

        auto conductivity_ =  1e-06*gmax*O;
        auto current       =  conductivity_*(v-Erev);
        ::arb::gpu::reduce_by_key(_pp_var_weight[tid_]*conductivity_,_pp_var_vec_g, node_indexi_, lane_mask_);
        ::arb::gpu::reduce_by_key(_pp_var_weight[tid_]*current_,_pp_var_vec_i, node_indexi_, lane_mask_);
    }
}

__global__
void apply_events(arb_mechanism_ppack params_, arb_deliverable_event_stream stream) {
    PPACK_IFACE_BLOCK;
    auto tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid_<stream.n_streams) {
        auto begin = stream.events + stream.begin[tid_];
        auto end   = stream.events + stream.end[tid_];
        for (auto p = begin; p<end; ++p) {
            if (p->mech_id==_pp_var_mechanism_id) {
                auto tid_ = p->mech_index;
                auto weight = p->weight;
		auto node_indexi_ = _pp_var_node_index[tid_];
                auto vec_dii_     = _pp_var_vec_di[node_indexi_];

                arb_value_type t = _pp_var_vec_t[vec_dii_];
                if (!_pp_var_on[tid_]) {

                    // Read
                    auto z = _pp_var_z[tid_];
                    auto y = _pp_var_y[tii_];
                    auto x = _pp_var_x[tii_];
                    auto u = _pp_var_u[tii_];

                    const auto tsyn      = _pp_var_tsyn[tid_];
                    const auto tau_rec   = _pp_var_tau_rec[tid_];
                    const auto tau_1     = _pp_var_tau_1[tid_];
                    const auto tau_facil = _pp_var_tau_facil[tid_];
                    const auto U         = _pp_var_U[tid_];
                    const auto Tmax      = _pp_var_Tmax[tid_];
                    const auto numpulses = (int)_pp_var_numpulses[tid_];

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
                    _pp_var_T[tid_] = Tmax*y;
                    _pp_var_z[tid_] = z;
                    _pp_var_y[tii_] = y;
                    _pp_var_x[tii_] = x;
                    _pp_var_u[tii_] = u;

                    // START EDIT
                    auto pulse = (numpulses%50); // rolling window update
                    _pp_var_tspike[pulse][tid_] = t; 
                    _pp_var_PRE[pulse][tid_]    = y; 

                    // END EDIT

                    _pp_var_numpulses[tid_] = numpulses + 1.;
                    _pp_var_tsyn[tid_] = t;
                }
                _pp_var_delay[tid_] = _pp_var_Cdur[tid_];
            }
        }
    }
}
} // namespace

void mechanism_Ampa_gpu_init_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    init<<<grid_dim, block_dim>>>(*p);
    if (!p->multiplicity) return;
    multiply<<<dim3{grid_dim, 4}, block_dim>>>(*p);
}

void mechanism_Ampa_gpu_compute_currents_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}

void mechanism_Ampa_gpu_advance_state_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    advance_state<<<grid_dim, block_dim>>>(*p);
}

void mechanism_Ampa_gpu_write_ions_(arb_mechanism_ppack* p) {}

void mechanism_Ampa_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_Ampa_gpu_apply_events_(arb_mechanism_ppack* p, arb_deliverable_event_stream* stream_ptr) {
    auto n = stream_ptr->n_streams;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    apply_events<<<grid_dim, block_dim>>>(*p, *stream_ptr);
}

} // namespace Ampa_catalogue
} // namespace arb
