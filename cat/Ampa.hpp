#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
  arb_mechanism_type make_arb_Ampa_catalogue_Ampa() {
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
        { "tspike0", "", NAN, -1e10, 1e10},
        { "tspike1", "", NAN, -1e10, 1e10},
        { "tspike2", "", NAN, -1e10, 1e10},
        { "tspike3", "", NAN, -1e10, 1e10},
        { "tspike4", "", NAN, -1e10, 1e10},
        { "tspike5", "", NAN, -1e10, 1e10},
        { "tspike6", "", NAN, -1e10, 1e10},
        { "tspike7", "", NAN, -1e10, 1e10},
        { "tspike8", "", NAN, -1e10, 1e10},
        { "tspike9", "", NAN, -1e10, 1e10},
        { "tspike10", "", NAN, -1e10, 1e10},
        { "tspike11", "", NAN, -1e10, 1e10},
        { "tspike12", "", NAN, -1e10, 1e10},
        { "tspike13", "", NAN, -1e10, 1e10},
        { "tspike14", "", NAN, -1e10, 1e10},
        { "tspike15", "", NAN, -1e10, 1e10},
        { "tspike16", "", NAN, -1e10, 1e10},
        { "tspike17", "", NAN, -1e10, 1e10},
        { "tspike18", "", NAN, -1e10, 1e10},
        { "tspike19", "", NAN, -1e10, 1e10},
        { "tspike20", "", NAN, -1e10, 1e10},
        { "tspike21", "", NAN, -1e10, 1e10},
        { "tspike22", "", NAN, -1e10, 1e10},
        { "tspike23", "", NAN, -1e10, 1e10},
        { "tspike24", "", NAN, -1e10, 1e10},
        { "tspike25", "", NAN, -1e10, 1e10},
        { "tspike26", "", NAN, -1e10, 1e10},
        { "tspike27", "", NAN, -1e10, 1e10},
        { "tspike28", "", NAN, -1e10, 1e10},
        { "tspike29", "", NAN, -1e10, 1e10},
        { "tspike30", "", NAN, -1e10, 1e10},
        { "tspike31", "", NAN, -1e10, 1e10},
        { "tspike32", "", NAN, -1e10, 1e10},
        { "tspike33", "", NAN, -1e10, 1e10},
        { "tspike34", "", NAN, -1e10, 1e10},
        { "tspike35", "", NAN, -1e10, 1e10},
        { "tspike36", "", NAN, -1e10, 1e10},
        { "tspike37", "", NAN, -1e10, 1e10},
        { "tspike38", "", NAN, -1e10, 1e10},
        { "tspike39", "", NAN, -1e10, 1e10},
        { "tspike40", "", NAN, -1e10, 1e10},
        { "tspike41", "", NAN, -1e10, 1e10},
        { "tspike42", "", NAN, -1e10, 1e10},
        { "tspike43", "", NAN, -1e10, 1e10},
        { "tspike44", "", NAN, -1e10, 1e10},
        { "tspike45", "", NAN, -1e10, 1e10},
        { "tspike46", "", NAN, -1e10, 1e10},
        { "tspike47", "", NAN, -1e10, 1e10},
        { "tspike48", "", NAN, -1e10, 1e10},
        { "tspike49", "", NAN, -1e10, 1e10},
        { "PRE0"   , "", NAN, -1e10, 1e10},
        { "PRE1"   , "", NAN, -1e10, 1e10},
        { "PRE2"   , "", NAN, -1e10, 1e10},
        { "PRE3"   , "", NAN, -1e10, 1e10},
        { "PRE4"   , "", NAN, -1e10, 1e10},
        { "PRE5"   , "", NAN, -1e10, 1e10},
        { "PRE6"   , "", NAN, -1e10, 1e10},
        { "PRE7"   , "", NAN, -1e10, 1e10},
        { "PRE8"   , "", NAN, -1e10, 1e10},
        { "PRE9"   , "", NAN, -1e10, 1e10},
        { "PRE10"   , "", NAN, -1e10, 1e10},
        { "PRE11"   , "", NAN, -1e10, 1e10},
        { "PRE12"   , "", NAN, -1e10, 1e10},
        { "PRE13"   , "", NAN, -1e10, 1e10},
        { "PRE14"   , "", NAN, -1e10, 1e10},
        { "PRE15"   , "", NAN, -1e10, 1e10},
        { "PRE16"   , "", NAN, -1e10, 1e10},
        { "PRE17"   , "", NAN, -1e10, 1e10},
        { "PRE18"   , "", NAN, -1e10, 1e10},
        { "PRE19"   , "", NAN, -1e10, 1e10},
        { "PRE20"   , "", NAN, -1e10, 1e10},
        { "PRE21"   , "", NAN, -1e10, 1e10},
        { "PRE22"   , "", NAN, -1e10, 1e10},
        { "PRE23"   , "", NAN, -1e10, 1e10},
        { "PRE24"   , "", NAN, -1e10, 1e10},
        { "PRE25"   , "", NAN, -1e10, 1e10},
        { "PRE26"   , "", NAN, -1e10, 1e10},
        { "PRE27"   , "", NAN, -1e10, 1e10},
        { "PRE28"   , "", NAN, -1e10, 1e10},
        { "PRE29"   , "", NAN, -1e10, 1e10},
        { "PRE30"   , "", NAN, -1e10, 1e10},
        { "PRE31"   , "", NAN, -1e10, 1e10},
        { "PRE32"   , "", NAN, -1e10, 1e10},
        { "PRE33"   , "", NAN, -1e10, 1e10},
        { "PRE34"   , "", NAN, -1e10, 1e10},
        { "PRE35"   , "", NAN, -1e10, 1e10},
        { "PRE36"   , "", NAN, -1e10, 1e10},
        { "PRE37"   , "", NAN, -1e10, 1e10},
        { "PRE38"   , "", NAN, -1e10, 1e10},
        { "PRE39"   , "", NAN, -1e10, 1e10},
        { "PRE40"   , "", NAN, -1e10, 1e10},
        { "PRE41"   , "", NAN, -1e10, 1e10},
        { "PRE42"   , "", NAN, -1e10, 1e10},
        { "PRE43"   , "", NAN, -1e10, 1e10},
        { "PRE44"   , "", NAN, -1e10, 1e10},
        { "PRE45"   , "", NAN, -1e10, 1e10},
        { "PRE46"   , "", NAN, -1e10, 1e10},
        { "PRE47"   , "", NAN, -1e10, 1e10},
        { "PRE48"   , "", NAN, -1e10, 1e10},
        { "PRE49"   , "", NAN, -1e10, 1e10},
    };
    static arb_size_type n_state_vars = 125;  // add 100
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

  arb_mechanism_interface* make_arb_Ampa_catalogue_Ampa_interface_multicore();
  arb_mechanism_interface* make_arb_Ampa_catalogue_Ampa_interface_gpu() { return nullptr; }
}
