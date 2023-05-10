include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

v = 1e-15

r0 = 0.087
p = 0.76
r_mw = 2.7e6
tot_prot = 15e-12

Φ0 = p*r0

g_ribos = Φ0*tot_prot # unit grams

mol_ribos = g_ribos/r_mw # unit mol

conc_ribos = mol_ribos/v # unit M

rh0 = conc_ribos*1e6 # unit μM

kdam = kin/rh0



