using DifferentialEquations, StaticArrays, DataFrames, Plots, OrderedCollections, Statistics#, PlotlyJS
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

θ = 4.38; max = 1260; thr = 7;

# to go from molecules/cell to μM/min
SF = 1e6/(6.022e23*1e-15)

a_SF = 22000

θtscr = θ*SF*a_SF
θtlr = thr*SF*a_SF

# max rate of translation 
g_max = max*SF


# ribosomes
ribo_cell =  [6800,13500,26300,45100,72000]
ribo_conc = []
for i in ribo_cell
    push!(ribo_conc, i*SF)
end
print(ribo_conc)

#checking gr_c values from real data 
grs = [0.6,1,1.5,2,2.5]/60
print(grs)
cs = []
for (i,j) in zip(grs, ribo_conc)
    push!(cs, i/j)
end
print(cs)
growth_rate_constant = mean(cs)

# growth rate constant
lam = 2.77/60 # taken from colD data in google colab notebook - need more WT data to base this off I think
rh = ribo_conc[end] # taken from ribosome table at the highest growth rate 
gr_c = lam/rh

#initial value of rh 
growth_rate = 2.5/60 # /60 to get from per hour to per minute
rh_0 = growth_rate/growth_rate_constant

# influx of healthy ribosomes
kin = lam*rh/g_max

# translation rate calculation for kdeg 



params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]
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


# steady state value of RtcR to use as initial value 
rtcr_init = ω_r*g_max/((rh_0*nr*gr_c^2)+d*nr*gr_c)
ω_r = rtcr_0/(g_max/((rh_0*nr*gr_c^2)+d*nr*gr_c))
