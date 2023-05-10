using DifferentialEquations, StaticArrays, DataFrames, Plots, OrderedCollections, Statistics#, PlotlyJS
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

# translation rate calculation for kdeg 
params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]
solu = sol(rtc_model, init, tspan, params)
rh = get_ssval(solu, :rh)
rm_a = get_ssval(solu, :rm_a)
rt = get_ssval(solu, :rt)
rtcr = get_ssval(solu, :rtcr)
alpha = rt/kr 
fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
ra = fa*rtcr
Vinit = ra*Vmax_init*atp/(Km_init+atp)
tscr_el_ab = ω_ab*atp/(θtscr+atp)
tscr_ab = Vinit*tscr_el_ab
tlr_el = g_max*atp/(θtlr+atp)
tlr = rh*rm_a*tlr_el # using ss values here but don't know parameters? 

tscr_ab+tlr

d





