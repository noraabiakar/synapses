This repository contains 2 versions of the same original mechanism [Ampa](https://senselab.med.yale.edu/ModelDB/ShowModel?model=128446&file=%2fShortTermPlasticityMFGRC_Nieus2006%2fAmpa.mod#tabs-2) obtained from modeldb.
1. cat_gt

   This directory contains the "ground truth" version of the mechanism, written against arbor's mechanism ABI in C++. It can be used to generate a mechanism catalogue (which can then be loaded in arbor) using `arbor-build-catalogue Ampa cat_gt/ --raw Ampa -g cuda -C`.
   
   `cat_gt/used_nmodl/AMPA.mod` was used as the starting point for writing the ground truth mechanism. It differs from the original mechanism in the following (major) ways:
   - Arrays `tspike` and `PRE` are commented out, along with the definition and usage of the `NTdiffusion` function. This is because arbor's modcc doesn't support arrays and for loops.
   - `on`, `t0`, `x`, `y`, `u` are ASSIGNED instead of arguments to the `NET_RECEIVE` function. This means that these variables are not shared among connections to the same synapse. This is needed because arbor doesn't support connection-local variables.
   - `delay` is added as a new state variable. It is used to replace the `net_send` functionality. It decays linearly and acts as a timer delaying the effect of a spike, similar to `net_send` in the original code.

   `cat_gt/used_nmodl/AMPA.mod` was first compiled using `modcc` by running `modcc -t gpu -t cpu -N arb::Ampa_catalogue -o Ampa AMPA.mod`, then the generated files were editted in the following (major) ways:
   - The arrays `tspike` and `PRE` as well as the definition and usage of the `NTdiffusion` function were added again.
   - The code was simplified for readability.


   Writing mechanisms this way is not recommended as it requires C++ knowledge as well as sufficient knowledge of the mechanism ABI. However, it can be used to verify the performance of the third version of the mechanism `AMPA_approx.mod`

3. cat_approx

   This directory contains the "approximate" version of the mechanism, written in nmodl. It can be used to generate a mechanism catalogue (which can then be loaded in arbor) using `arbor-build-catalogue Ampa cat_approx/ -g cuda -C`.
   
   `cat_approx/AMPA.mod` differs from the original mechanism in the following (major) ways:
   - A kinetic scheme `gstate` is added to approximate the effect of the `NTdiffusion` function. It works on 3 state variables `g_active`, `g_emission`
   ```
    g_emission <-> g_active   (alpha, beta)
    g_active   <-> (gamma, 0)
   ```
   The first equation represents the diffusion of glutamate between the emission area and the active area.
   The second equation represents the diffusion of glutamate between the active area and the surrounds.

   Values for `alpha`, `beta`, and `gamma`need to be determined using parameter optimization. `ampa_cell.py` and `curve_fit.py` are used to determine these parameters. 
