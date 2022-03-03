#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
  arb_mechanism_type make__Ampa() {
    // Tables
    static arb_field_info globals[] = {{ "u0", "1", 0, 0, 1 } };
    static arb_size_type n_globals = 1;
    static arb_field_info state_vars[] = {
        { "C", "", NAN, -1e10, 1e10 },
        { "O", "", NAN, -1e10, 1e10 },
        { "D", "", NAN, -1e10, 1e10 },
        { "delay", "ms", NAN, -1e10, 1e10 },
        { "v", "mV", NAN, -1e10, 1e10 },
        { "g", "pS", NAN, -1e10, 1e10 },
        { "T", "mM", NAN, -1e10, 1e10 },
        { "r1", "/ ms", NAN, -1e10, 1e10 },
        { "r6", "/ ms", NAN, -1e10, 1e10 },
        { "Trelease", "mM", NAN, -1e10, 1e10 },
        { "x", "", NAN, -1e10, 1e10 },
        { "tsyn", "ms", NAN, -1e10, 1e10 },
        { "Mres", "mM", NAN, -1e10, 1e10 },
        { "NTdiffusion", "mM", NAN, -1e10, 1e10 },
        { "numpulses", "", NAN, -1e10, 1e10 },
        { "xview", "", NAN, -1e10, 1e10 },
        { "yview", "", NAN, -1e10, 1e10 },
        { "zview", "", NAN, -1e10, 1e10 },
        { "Pview", "", NAN, -1e10, 1e10 },
        { "on", "", NAN, -1e10, 1e10 },
        { "nspike", "", NAN, -1e10, 1e10 },
        { "t0", "ms", NAN, -1e10, 1e10 },
        { "y", "", NAN, -1e10, 1e10 },
        { "z", "", NAN, -1e10, 1e10 },
        { "u", "", NAN, -1e10, 1e10 },
        { "tsyn", "ms", NAN, -1e10, 1e10}, 
        { "tspike", "", NAN, -1e10, 1e10, 50} // 50 refers to size of array 
        { "PRE"   , "", NAN, -1e10, 1e10, 50} // 50 refers to size of array 
    };
    static arb_size_type n_state_vars = 126;  // add 100
    static arb_field_info parameters[] = {
        { "gmax", "pS", 1200, -1e10, 1e10 },
        { "Cdur", "ms", 0.3, -1e10, 1e10 },
        { "Erev", "mV", 0, -1e10, 1e10 },
        { "kB", "mM", 0.44, -1e10, 1e10 },
        { "r1FIX", "/ ms / mM", 5.4, -1e10, 1e10 },
        { "r6FIX", "/ ms / mM", 1.12, -1e10, 1e10 },
        { "r2", "/ ms", 0.82, -1e10, 1e10 },
        { "r5", "/ ms", 0.013, -1e10, 1e10 },
        { "tau_1", "ms", 3, 1e-9, 1e9 },
        { "tau_rec", "ms", 35.1, 1e-9, 1e9 },
        { "tau_facil", "ms", 10.8, 0, 1e9 },
        { "U", "1", 0.416, 0, 1 },
        { "Tmax", "mM", 1, -1e10, 1e10 },
        { "M", "", 21500, -1e10, 1e10 },
        { "R", "um", 1.033, -1e10, 1e10 },
        { "Diff", "um2 / ms", 0.223, -1e10, 1e10 },
        { "lamd", "nm", 20, -1e10, 1e10} 
    };
    static arb_size_type n_parameters = 17;
    static arb_ion_info* ions = NULL;
    static arb_size_type n_ions = 0;

    arb_mechanism_type result;
    result.abi_version=ARB_MECH_ABI_VERSION;
    result.fingerprint="<placeholder>";
    result.name="Ampa";
    result.kind=arb_mechanism_kind_point;
    result.is_linear=false;
    result.has_post_events=false;
    result.globals=globals;
    result.n_globals=n_globals;
    result.ions=ions;
    result.n_ions=n_ions;
    result.state_vars=state_vars;
    result.n_state_vars=n_state_vars;
    result.parameters=parameters;
    result.n_parameters=n_parameters;
    return result;
  }

  arb_mechanism_interface* make__Ampa_interface_multicore();
  arb_mechanism_interface* make__Ampa_interface_gpu() { return nullptr; }
}
