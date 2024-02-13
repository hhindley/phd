using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl")
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_inhibition_model.jl");

@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

# STILL HAVE BISTABILITY SO MOST LIKELY MEANS THAT THE RH POSITIVE FEEDBACK LOOP IS THE ONE CONTRIBUTING TO BISTABILITY 

tspan = (0,1e9)

params2 = @LArray [0, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu3 = sol(rtc_model, initial, tspan, params2)
p3 = plotly_plot_sol(solu3, "", "", "L = 0, no damage");

atp_range = range(0,5000, length=1000)
res = change_param(atp_range, :atp, rtc_model, initial, get_ssval, all_species, params2)

plot([scatter(x=atp_range, y=res[:rh]),scatter(x=atp_range, y=res[:rm_a]),scatter(x=atp_range, y=res[:rtca]),scatter(x=atp_range, y=res1[:rh]),scatter(x=atp_range, y=res1[:rm_a]),scatter(x=atp_range, y=res1[:rtca])])


params1 = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "orig, higher damage")

res1 = change_param(atp_range, :atp, rtc_model, initial, get_ssval, all_species, params1)

plot([scatter(x=atp_range, y=res1[:rh]),scatter(x=atp_range, y=res1[:rm_a]),scatter(x=atp_range, y=res1[:rtca])])





params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

params3 = (L = 0., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

# params4 = (L = 10., c = 0.001, kr = 0.9, Vmax_init = 39.51, Km_init = 250.,
# θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
# krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
# kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
# kdam =  0.01, lam = 0.014)

br = get_br(rtc_mod, params2, initial, 3.)
bs = plot_all_curves_bistable(br, colors2, colorsr, "Original, L = 10")

br2 = get_br(rtc_mod, params3, initial, 3.)
bs2 = plot_all_curves_bistable(br2, colors2, colorsr, "L = $(params3.L)")

# br3 = get_br(rtc_mod, params4, initial, 3.)
# bs3 = plot_all_curves_bistable(br3, colors2, colorsr, "Large L, L = $(params3.L)")
nocarf = [bs bs2]
