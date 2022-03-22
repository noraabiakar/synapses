TITLE AMPA

COMMENT
ENDCOMMENT

NEURON {
    POINT_PROCESS Ampa

    NONSPECIFIC_CURRENT i
    RANGE Cdur, Erev, g, gmax, kB
    RANGE r1FIX, r6FIX, r2, r5
    RANGE tau_1, tau_rec, tau_facil, U
    RANGE T, Tmax
    RANGE M, Diff, R, lamd
}

UNITS {
    (nA)     = (nanoamp)
    (mV)     = (millivolt)
    (umho)   = (micromho)
    (mM)     = (milli/liter)
    (pS)     = (picosiemens)
    (nS)     = (nanosiemens)
    (um)     = (micrometer)
}

CONSTANT {
    PI = 3.14159
}

PARAMETER {
    : postsynaptic parameters
    gmax        = 1200  (pS)
    Cdur        = 0.3   (ms)
    Erev        = 0     (mV)
    kB          = 0.44  (mM)

    r1FIX       = 5.4   (/ms/mM)
    r6FIX       = 1.12  (/ms/mM)
    r2          = 0.82  (/ms)
    r5          = 0.013 (/ms)

    : presynaptic parameters
    tau_1         = 3    (ms)     < 1e-9, 1e9 >
    tau_rec       = 35.1 (ms)     < 1e-9, 1e9 >
    tau_facil     = 10.8 (ms)     < 0, 1e9 >

    U             = 0.416 (1)     < 0, 1 >
    u0            = 0     (1)     < 0, 1 >
    Tmax          = 1     (mM)

    : Diffusion
    M            = 21500
    R             = 1.033 (um)
    Diff          = 0.223 (um2/ms)
    lamd          = 20    (nm)
    t
}


ASSIGNED {
    v         (mV)
    g         (pS)
    T         (mM)

    x
    tsyn      (ms)

    Mres      (mM)
    on
    t0        (ms)
    y
    z
    u
    tsyn      (ms)
}

STATE {
    C
    O
    D
    delay      (ms)

    g_emission : glutamate concentration at the emission area
    g_active   : glutamate concentration at the active area
    g_surround : glutamate concentration at the surroundings (elsewhere)
}

INITIAL {
    C = 1
    O = 0
    D = 0
    delay = 0

    T = 0
    Mres = (1e3 * 1e15 / 6.022e23 * M)   : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15

    on = 0
    y = 0
    z = 0
    u = u0
    tsyn = t
}

BREAKPOINT {
    SOLVE gstates METHOD sparse
    SOLVE kstates METHOD sparse
    SOLVE sdelay  METHOD cnexp

    g = gmax * O
    i = (1e-6) * g * (v-Erev)

    if (delay == 0) {
        t0 = t
        T = 0
        on = 0
    }
}

KINETIC kstates {
    LOCAL Trelease, r1, r6

    Trelease = T + g_active
    r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
    r6 = r6FIX * Trelease^2 / (Trelease + kB)^2

    ~ C  <-> O (r1,r2)
    ~ D  <-> C (r5,r6)
    CONSERVE C+O+D = 1
}

KINETIC gstates {
    LOCAL alpha, beta, gamma, delta, epsilon

    alpha = 0
    beta = 0
    gamma = 0
    delta = 0
    epsilon = 0

    ~ g_emission <-> g_active   (alpha, beta)
    ~ g_active   <-> g_surround (gamma, delta)
    ~ g_surround <-> (epsilon, 0)
}

DERIVATIVE sdelay {
   delay' = -1
}

NET_RECEIVE(weight) {
    : presynaptic modulation
    if (on != 1) {
        t0 = t
        on = 1

        z = z*exp( - (t - tsyn) / tau_rec )
        z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec))/((tau_1/tau_rec)-1) )
        y = y*exp(-(t - tsyn)/tau_1)
        x = 1-y-z

        if (tau_facil > 0) {
            u = u*exp(-(t - tsyn)/tau_facil)
            u = u + U * ( 1 - u )
        }
        else {
            u = U
        }
        y = y + x * u

        T = Tmax * y
        tsyn = t
    }
    delay = Cdur
}
