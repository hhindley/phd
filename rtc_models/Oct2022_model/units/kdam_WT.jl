include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

M = 1e-15

r0 = 0.087
p = 0.76

Φ0 = p*r0

tlr_el = g_max*atp/(θtlr+atp)

kdam = (tlr_el*kin/M)*(1/Φ0)

kdam = (kin/M)/r0 

kdam = kin/Φ0

kdam = kin/r0