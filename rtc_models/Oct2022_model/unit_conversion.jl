using DifferentialEquations, StaticArrays, DataFrames, Plots, OrderedCollections, Statistics#, PlotlyJS
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

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
lambda = grs[end]
rh =  ribo_conc[end] # taken from ribosome table at the highest growth rate 
gr_c = lam/rh

#initial value of rh 
growth_rate = 2.5/60 # /60 to get from per hour to per minute
rh_0 = growth_rate/growth_rate_constant

# influx of healthy ribosomes
kin = lambda*rh/g_max

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




# steady state value of RtcR to use as initial value 
rtcr_init = ω_r*g_max/((rh_0*nr*gr_c^2)+d*nr*gr_c)
ω_r = rtcr_0/(g_max/((rh_0*nr*gr_c^2)+d*nr*gr_c))



# extrapolate to get value of OD_0
using Plots, Dierckx
time1 = dfe[:,1]
spl = Spline1D(time1, WT2; k=1, bc="extrapolate")
OD_0 = evaluate(spl, 0)

time2 = dff[:,1]
spl2 = Spline1D(time2, WT3; k=1, bc="extrapolate")
OD_02 = evaluate(spl2, 0)

time3 = df2[:,1]
spl3 = Spline1D(time3, WT4; k=1, bc="extrapolate")
OD_03 = evaluate(spl3, 0)

spl4 = Spline1D(time2, hpx_rtcon; k=1, bc="extrapolate")
OD_0_hpxon = evaluate(spl4, 0)

spl5 = Spline1D(time3, WT_colD; k=1, bc="extrapolate")
OD_0_wtcolD = evaluate(spl5, 0)

spl6 = Spline1D(time3, nA_colD; k=1, bc="extrapolate")
OD_04_nAcolD = evaluate(spl6, 0)

spl7 = Spline1D(time3, nB_colD; k=1, bc="extrapolate")
OD_04_nB_colD = evaluate(spl7, 0)

spl8 = Spline1D(time3, nB_B_colD; k=1, bc="extrapolate")
OD_04_nBBcolD = evaluate(spl8, 0)

spl9 = Spline1D(time3, nB_Bmut_colD; k=1, bc="extrapolate")
OD_04_bBBmutcolD = evaluate(spl9, 0)

spl10 = Spline1D(time3, nR_colD; k=1, bc="extrapolate")
OD_04_nRcolD = evaluate(spl0, 0)


function check_OD_0(time, data, OD)
    p1 = plot(time*60, data, markershape=:circle)
    scatter!((0,OD));
    time2 = [[0];time]
    data2 = [[OD];data]
    p2 = plot(time2*60, data2, markershape=:circle)
    return p1, p2
end

p1, p2 = check_OD_0(time1, WT2, OD_0)
plot(p1,p2, layout=(2,1))

p3, p4 = check_OD_0(time2, WT3, OD_02)
plot(p3,p4, layout=(2,1))

p5, p6 = check_OD_0(time3, WT4, OD_03)
plot(p5,p6, layout=(2,1))
