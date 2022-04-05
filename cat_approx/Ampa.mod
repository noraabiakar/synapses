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
    RANGE M, r, h, a_ratio, b_ratio, c_ratio
    RANGE g_emission, g_active, delay 
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
    t

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
    M             = 21500
    h             = 20              (um)    : Synaptic cleft height
    r             = 1.033           (um)    <0, 3>       : needs to be obtained by BluPyOpt
    a_ratio       = 3.92417130e-01  (um)    <-1e9, 1e9>  : needs to be obtained by BluPyOpt
    b_ratio       = -1.66705011e-01 (um)    <-1e9, 1e9>  : needs to be obtained by BluPyOpt
    c_ratio       = 1.07668617e-05  (um)    <-1e9, 1e9>  : needs to be obtained by BluPyOpt
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
}

STATE {
    C
    O
    D
    delay      (ms)

    g_emission : glutamate concentration at the emission area
    g_active   : glutamate concentration at the active area

    Trelease
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

    g_emission = 0
    g_active = 0
}

BREAKPOINT {
    SOLVE gstates METHOD sparse
    SOLVE kstates METHOD sparse
    SOLVE sdelay  METHOD cnexp

    g = gmax * O
    i = (1e-6) * g * (v-Erev)

    if (delay <= 0) {
        t0 = t
        T = 0
        on = 0
    }
}

KINETIC kstates {
    LOCAL r1, r6

    Trelease = T + g_active
    r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
    r6 = r6FIX * Trelease^2 / (Trelease + kB)^2

    ~ C  <-> O (r1,r2)
    ~ D  <-> C (r5,r6)
    CONSERVE C+O+D = 1
}

KINETIC gstates {
    LOCAL emission_area, active_area, alpha, beta, gamma, delta, epsilon

    emission_area = 4*PI*r^2
    active_area   = 4*PI*(r+h)^2

    alpha = a_ratio*emission_area
    beta  = b_ratio*emission_area
    gamma = c_ratio*active_area

    ~ g_emission <-> g_active (alpha, beta)
    ~ g_active   <->          (gamma, 0)
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
        g_emission = g_emission + Mres
        tsyn = t
    }
    delay = Cdur
}
