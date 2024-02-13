using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
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

tspan = (0,1e9)
k_inhib1 = 0.4#1 #10 #1
k_inhib2 = 0.025
inhib = 0.1 #1.5 #0.2
params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014, k_inhib1, k_inhib2, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :k_inhib1, :k_inhib2, :inhib)

initial_i = [0,0,0,0,0,0,11.29,0,0,0]


solu_rtcb = sol(rtc_inhib_model_rtcb, initial_i, tspan, params_inhib)
p1_rtcb = plotly_plot_sol(solu_rtcb, "", "", "RtcB inhib");

solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu_rtca, "", "", "RtcA inhib");

solu_rtcr = sol(rtc_inhib_model_rtcr, initial_i, tspan, params_inhib)
p1_rtcr = plotly_plot_sol(solu_rtcr, "", "", "RtcR inhib");


params1 = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "orig, higher damage")

params2 = @LArray [0, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu3 = sol(rtc_model, initial, tspan, params2)
p3 = plotly_plot_sol(solu3, "", "", "L = 0, higher damage");

nocarf = [p2 p3]
p_species_inhib = [p2 p1_rtca p1_rtcb p1_rtcr]

params_for_ssval_setup_inhib = (L = 50000., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib)

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
largeL = [bs bs2]

br_rtca = get_br(rtc_inhib_mod_rtca, params_for_ssval_setup_inhib, initial_i, 3.)
rtca_bs = plot_all_curves_bistable(br_rtca, colors2, colorsr, "RtcA inhibited");

br_rtcb = get_br(rtc_inhib_mod_rtcb, params_for_ssval_setup_inhib, initial_i, 3.)
rtcb_bs = plot_all_curves_bistable(br_rtcb, colors2, colorsr, "RtcB inhibited");

br_rtcr = get_br(rtc_inhib_mod_rtcr, params_for_ssval_setup_inhib, initial_i, 3.)
rtcr_bs = plot_all_curves_bistable(br_rtcr, colors2, colorsr, "RtcR inhibited")

br_rt = get_br(rtc_inhib_mod_rt, params_for_ssval_setup_inhib, initial_i, 3.)
rt_bs = plot_all_curves_bistable(br_rt, colors2, colorsr, "Rt inhibited");

init_i = [0,0,0,0,0,0,11.29,0,0,0,0]

br_rtcab = get_br(rtc_inhib_mod_rtcab, params_for_ssval_setup_inhib, init_i, 3.)
rtcab_bs = plot_all_curves_bistable(br_rtcab, colors2, colorsr, "RtcAB inhibited")

br_rtcbr = get_br(rtc_inhib_mod_rtcbr, params_for_ssval_setup_inhib, init_i, 3.)
rtcbr_bs = plot_all_curves_bistable(br_rtcbr, colors2, colorsr, "RtcBR inhibited")

br_rtcar = get_br(rtc_inhib_mod_rtcar, params_for_ssval_setup_inhib, init_i, 3.)
rtcar_bs = plot_all_curves_bistable(br_rtcar, colors2, colorsr, "RtcAR inhibited")

br_rtcbt = get_br(rtc_inhib_mod_rtcbt, params_for_ssval_setup_inhib, init_i, 3.)
rtcbt_bs = plot_all_curves_bistable(br_rtcbt, colors2, colorsr, "RtcBT inhibited")

br_rtcat = get_br(rtc_inhib_mod_rtcat, params_for_ssval_setup_inhib, init_i, 3.)
rtcat_bs = plot_all_curves_bistable(br_rtcat, colors2, colorsr, "RtcAT inhibited")

br_rtcrt = get_br(rtc_inhib_mod_rtcrt, params_for_ssval_setup_inhib, init_i, 3.)
rtcrt_bs = plot_all_curves_bistable(br_rtcrt, colors2, colorsr, "RtcRT inhibited")



all_bs = [bs rtca_bs rtcb_bs rtcr_bs rt_bs]
combos = [bs rtcab_bs rtcbr_bs rtcar_bs; rtcbt_bs rtcat_bs rtcrt_bs]
open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_species_inhib_mid.html", "w") do io
    PlotlyBase.to_html(io, p_species_inhib.plot)
end
open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_all_bs_with_rt.html", "w") do io
    PlotlyBase.to_html(io, all_bs.plot)
end

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/no_carf_higherdam.html", "w") do io
    PlotlyBase.to_html(io, nocarf.plot)
end

savefig(ribos, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/ribos.png")

# bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt, first_r, middle_r, last_r, first, middle, last = get_all_curves_for_bistab_plotting(br, colors2, colors_r)

rf, rm, rl = plot_species_separately_ribosomes(br, colors_r, "Original")
rf1, rm1, rl1 = plot_species_separately_ribosomes(br2, colors_r1, "large L")
rf2, rm2, rl2 = plot_species_separately_ribosomes(br_rtca, colors_r2, "RtcA inhib")
rf3, rm3, rl3 = plot_species_separately_ribosomes(br_rtcb, colors_r3, "RtcB inhib")
rf4, rm4, rl4 = plot_species_separately_ribosomes(br_rtcr, colors_r4, "RtcR inhib")
rf5, rm5, rl5 = plot_species_separately_ribosomes(br_rtcab, colors_r5, "RtcAB inhib")
rf6, rm6, rl6 = plot_species_separately_ribosomes(br_rt, colors_r6, "Rt inhib")

rf7, rm7, rl7 = plot_species_separately_ribosomes(br_rtcar, colors_r7, "RtcAR inhib")
rf8, rm8, rl8 = plot_species_separately_ribosomes(br_rtcat, colors_r8, "RtcAT inhib")
rf9, rm9, rl9 = plot_species_separately_ribosomes(br_rtcbr, colors_r9, "RtcBR inhib")
rf10, rm10, rl10 = plot_species_separately_ribosomes(br_rtcbt, colors_r10, "RtcBT inhib")
rf11, rm11, rl11 = plot_species_separately_ribosomes(br_rtcrt, colors_r11, "RtcRT inhib")

plot([rf[1],rm[1],rl[1],rf2[1],rm2[1],rl2[1]])

rt = plot([rf[3],rm[3],rl[3],rf1[3],rm1[3],rl1[3],rf2[3],rm2[3],rl2[3],rf3[3],rm3[3],rl3[3],rf4[3],rm4[3],rl4[3],rf5[3],rm5[3],rl5[3],rf6[3],rm6[3],rl6[3],rf7[3],rm7[3],rl7[3],rf8[3],rm8[3],rl8[3],rf9[3],rm9[3],rl9[3],rf10[3],rm10[3],rl10[3], rf11[3], rm11[3],rl11[3]], Layout(title="Rt", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))
rh = plot([rf[1],rm[1],rl[1],rf1[1],rm1[1],rl1[1],rf2[1],rm2[1],rl2[1],rf3[1],rm3[1],rl3[1],rf4[1],rm4[1],rl4[1],rf6[1],rm6[1],rl6[1],rf5[1],rm5[1],rl5[1],rf7[1],rm7[1],rl7[1],rf8[1],rm8[1],rl8[1],rf9[1],rm9[1],rl9[1],rf10[1],rm10[1],rl10[1], rf11[1], rm11[1],rl11[1]], Layout(title="Healthy ribosomes", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))#, xaxis_type="log"))
rd = plot([rf[2],rm[2],rl[2],rf1[2],rm1[2],rl1[2],rf2[2],rm2[2],rl2[2],rf3[2],rm3[2],rl3[2],rf4[2],rm4[2],rl4[2],rf5[2],rm5[2],rl5[2],rf6[2],rm6[2],rl6[2],rf7[2],rm7[2],rl7[2],rf8[2],rm8[2],rl8[2],rf9[2],rm9[2],rl9[2],rf10[2],rm10[2],rl10[2], rf11[2], rm11[2],rl11[2]], Layout(title="Rd", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))
[rt rh rd]
sf, sm, sl = plot_species_separately(br, colors_s, "orig")
sf1, sm1, sl1 = plot_species_separately(br2, colors_s1, "large L")
sf2, sm2, sl2 = plot_species_separately(br_rtca, colors_s2, "rtca inhib")
sf3, sm3, sl3 = plot_species_separately(br_rtcb, colors_s3, "rtcb inhib")
sf4, sm4, sl4 = plot_species_separately(br_rtcr, colors_s4, "rtcr inhib")
sf5, sm5, sl5 = plot_species_separately(br_rtcab, colors_s5, "rtcAB inhib")
sf6, sm6, sl6 = plot_species_separately(br_rt, colors_s6, "Rt inhib")

sf7, sm7, sl7 = plot_species_separately(br_rtcar, colors_s7, "rtcAR inhib")
sf8, sm8, sl8 = plot_species_separately(br_rtcat, colors_s8, "rtcAT inhib")
sf9, sm9, sl9 = plot_species_separately(br_rtcbr, colors_s9, "rtcBR inhib")
sf10, sm10, sl10 = plot_species_separately(br_rtcbt, colors_s10, "rtcBT inhib")
sf11, sm11, sl11 = plot_species_separately(br_rtcrt, colors_s11, "rtcRT inhib")

rma = plot([sf[3],sm[3],sl[3],sf1[3],sm1[3],sl1[3],sf2[3],sm2[3],sl2[3],sf3[3],sm3[3],sl3[3],sf4[3],sm4[3],sl4[3],sf5[3],sm5[3],sl5[3],sf6[3],sm6[3],sl6[3],sf7[3],sm7[3],sl7[3],sf8[3],sm8[3],sl8[3],sf9[3],sm9[3],sl9[3],sf10[3],sm10[3],sl10[3],sf11[3],sm11[3],sl11[3]], Layout(title="mRNA RtcAB", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))
rtca = plot([sf[2],sm[2],sl[2],sf1[2],sm1[2],sl1[2],sf2[2],sm2[2],sl2[2],sf3[2],sm3[2],sl3[2],sf4[2],sm4[2],sl4[2],sf5[2],sm5[2],sl5[2],sf6[2],sm6[2],sl6[2],sf7[2],sm7[2],sl7[2],sf8[2],sm8[2],sl8[2],sf9[2],sm9[2],sl9[2],sf10[2],sm10[2],sl10[2],sf11[2],sm11[2],sl11[2]], Layout(title="RtcA", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))
rtcb = plot([sf[4],sm[4],sl[4],sf1[4],sm1[4],sl1[4],sf2[4],sm2[4],sl2[4],sf3[4],sm3[4],sl3[4],sf4[4],sm4[4],sl4[4],sf5[4],sm5[4],sl5[4],sf6[4],sm6[4],sl6[4],sf7[4],sm7[4],sl7[4],sf8[4],sm8[4],sl8[4],sf9[4],sm9[4],sl9[4],sf10[4],sm10[4],sl10[4],sf11[4],sm11[4],sl11[4]], Layout(title="RtcB", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))
rtcr = plot([sf[6],sm[6],sl[6],sf1[6],sm1[6],sl1[6],sf2[6],sm2[6],sl2[6],sf3[6],sm3[6],sl3[6],sf4[6],sm4[6],sl4[6],sf5[6],sm5[6],sl5[6],sf6[6],sm6[6],sl6[6],sf7[6],sm7[6],sl7[6],sf8[6],sm8[6],sl8[6],sf9[6],sm9[6],sl9[6],sf10[6],sm10[6],sl10[6],sf11[6],sm11[6],sl11[6]], Layout(title="RtcR", yaxis_title="Steady-state concentration (μM)", xaxis_title="Damage rate (min<sup>-1</sup>)"))


savefig(rt, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/rt_inhib.png")
savefig(rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/rh_inhib_all.svg")

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/inhib_rtcr.html", "w") do io
    PlotlyBase.to_html(io, rtcr.plot)
end

k_inhib1_1 = 1#1 #10 #1
k_inhib2_1 = 0.0025
inhib_1 = 0.1 #
params_for_ssval_setup_inhib1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1_1, k_inhib2=k_inhib2_1, inhib=inhib_1)

k_inhib1_2 = 0.4#1 #10 #1
k_inhib2_2 = 0.025
inhib_2 = 8 #
params_for_ssval_setup_inhib2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1_2, k_inhib2=k_inhib2_2, inhib=inhib_2)

br_rtca1 = get_br(rtc_inhib_mod_rtca, params_for_ssval_setup_inhib1, initial_i, 3.)
br_rtca2 = get_br(rtc_inhib_mod_rtca, params_for_ssval_setup_inhib2, initial_i, 3.)

rf2_1, rm2_1, rl2_1 = plot_species_separately_ribosomes(br_rtca1, colors_r2, "RtcA inhib")
rf2_2, rm2_2, rl2_2 = plot_species_separately_ribosomes(br_rtca2, colors_r2, "RtcA inhib")

plot([rf[1],rm[1],rl[1],rf2_2[1],rm2_2[1],rl2_2[1],rf3[1],rm3[1],rl3[1]])

plot([rf[1],rm[1],rl[1],rf2[1],rm2[1],rl2[1],rf2_1[1],rm2_1[1],rl2_1[1],rf2_2[1],rm2_2[1],rl2_2[1]])

sf2, sm2, sl2 = plot_species_separately(br_rtca, colors_s2, "rtca inhib")
sf2_1, sm2_1, sl2_1 = plot_species_separately(br_rtca1, colors_s2, "rtca inhib")
sf2_2, sm2_2, sl2_2 = plot_species_separately(br_rtca2, colors_s2, "rtca inhib")

plot([sf[2],sm[2],sl[2],sf2[2], sm2[2], sl2[2],sf2_1[2], sm2_1[2], sl2_1[2],sf2_2[2], sm2_2[2], sl2_2[2]])


plot([first2[3],middle2[3],last2[3]])
plot([first3[3],middle3[3],last3[3]])
plot([first4[3],middle4[3],last4[3]])

plot([first[3],middle[3],last[3],first1[3],middle1[3],last1[3],first2[3],middle2[3],last2[3]])
ribos = plot([first_r[1],middle_r[1],last_r[1], bfp_rh,
first_r[2],middle_r[2],last_r[2], bfp_rd,
first_r[3],middle_r[3],last_r[3], bfp_rt],
Layout(yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)", title="RtcA inhibited"))


df = create_br_df(br)
df_rtcr = create_br_df(br_rtcr)
df_l = create_br_df(br2)
df_kr = create_br_df(br3)
df_rt = create_br_df(br_rt)

alpha = @. df.rt/kr
fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
ra = @. fa*df.rtcr

alpha1 = @. df_rtcr.rt/kr
fa1 = @. (1+alpha1)^6/(L*((1+c*alpha1)^6)+(1+alpha1)^6)
ra1 = @. fa1*df_rtcr.rtcr

alpha2 = @. df_l.rt/kr
fa2 = @. (1+alpha2)^6/(params3.L*((1+c*alpha2)^6)+(1+alpha2)^6)
ra2 = @. fa2*df_l.rtcr

alpha3 = @. df_kr.rt/params4.kr
fa3 = @. (1+alpha3)^6/(params4.L*((1+c*alpha3)^6)+(1+alpha3)^6)
ra3 = @. fa3*df_kr.rtcr

alpha4 = @. df_rt.rt/kr
fa4 = @. (1+alpha4)^6/(L*((1+c*alpha4)^6)+(1+alpha4)^6)
ra4 = @. fa4*df_rt.rtcr


p1 = plot([scatter(x=df.kdam, y=fa, name="orig", legendgroup="1", line=attr(color=:blue)), scatter(x=df_rtcr.kdam, y=fa1, name="rtcr inhib", legendgroup="2", line=attr(color=:green)), scatter(x=df_l.kdam, y=fa2, name="large L", legendgroup="3", line=attr(color=:purple)), scatter(x=df_kr.kdam, y=fa3, name="large kr", legendgroup="4", line=attr(color=:gold)), scatter(x=df_rt.kdam, y=fa4, name="Rt inhib", legendgroup="5", line=attr(color=:pink))], Layout(title="fraction of active RtcR"))
p2 = plot([scatter(x=df.kdam, y=ra, name="orig", legendgroup="1", line=attr(color=:blue)), scatter(x=df_rtcr.kdam, y=ra1, name="rtcr inhib", legendgroup="2", line=attr(color=:green)), scatter(x=df_l.kdam, y=ra2, name="large L", legendgroup="3", line=attr(color=:purple)), scatter(x=df_kr.kdam, y=ra3, name="large kr", legendgroup="4", line=attr(color=:gold)), scatter(x=df_rt.kdam, y=ra4, name="Rt inhib", legendgroup="5", line=attr(color=:pink))], Layout(title="amount of active RtcR"))
p3 = plot([scatter(x=df.kdam, y=df.rtcr, name="orig", legendgroup="1", line=attr(color=:blue)), scatter(x=df_rtcr.kdam, y=df_rtcr.rtcr, name="rtcr inhib", legendgroup="2", line=attr(color=:green)), scatter(x=df_l.kdam, y=df_l.rtcr, name="large L", legendgroup="3", line=attr(color=:purple)), scatter(x=df_kr.kdam, y=df_kr.rtcr, name="large kr", legendgroup="4", line=attr(color=:gold)), scatter(x=df_rt.kdam, y=df_rt.rtcr, name="Rt inhib", legendgroup="5", line=attr(color=:pink))], Layout(title="Total RtcR"))
[p1 p2 p3]


p_alpha, p_fa, p_ra, p_vinit, p_tscr_a, p_tscr_b, p_tlr_a, p_tlr_b, p_tlr_r, p_rtca1, p_rtcb1, p_vrep, p_vdam, p_vtag = plot_all_vars(solu_orig)
p_alpha1, p_fa1, p_ra1, p_vinit1, p_tscr_a1, p_tscr_b1, p_tlr_a1, p_tlr_b1, p_tlr_r1, p_rtca11, p_rtcb11, p_vrep1, p_vdam1, p_vtag1 = plot_all_vars(solu_rtca)
p_alpha2, p_fa2, p_ra2, p_vinit2, p_tscr_a2, p_tscr_b2, p_tlr_a2, p_tlr_b2, p_tlr_r2, p_rtca12, p_rtcb12, p_vrep2, p_vdam2, p_vtag2 = plot_all_vars(solu_rtcb)

p_orig_all = plot([p_alpha, p_fa, p_ra, p_vinit, p_tscr_a, p_tscr_b, p_tlr_a, p_tlr_b, p_tlr_r, p_rtca1, p_rtcb1, p_vrep, p_vdam, p_vtag], Layout(title="orig",xaxis=attr(range=(0, 800))))
p_rtca_all = plot([p_alpha1, p_fa1, p_ra1, p_vinit1, p_tscr_a1, p_tscr_b1, p_tlr_a1, p_tlr_b1, p_tlr_r1, p_rtca11, p_rtcb11, p_vrep1, p_vdam1, p_vtag1], Layout(title="RtcA",xaxis=attr(range=(0, 800))))
p_rtcb_all = plot([p_alpha2, p_fa2, p_ra2, p_vinit2, p_tscr_a2, p_tscr_b2, p_tlr_a2, p_tlr_b2, p_tlr_r2, p_rtca12, p_rtcb12, p_vrep2, p_vdam2, p_vtag2], Layout(title="RtcB",xaxis=attr(range=(0, 800))))

p_all_inhib = [p_orig_all p_rtca_all p_rtcb_all]

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_species_inhib.html", "w") do io
    PlotlyBase.to_html(io, p_species_inhib.plot)
end
open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_all_inhib.html", "w") do io
    PlotlyBase.to_html(io, p_all_inhib.plot)
end

solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu, "", "", "RtcA inhib")


solu_rtcb = sol(rtc_inhib_model_rtcb, initial_i, tspan, params_inhib)
p1_rtcb = plotly_plot_sol(solu, "", "", "RtcB inhib");

atp = 3578.9473684210525
rtcb = get_curve(solu, :rtcb)
rt = get_curve(solu, :rt)
rtcb_i = get_curve(solu, :rt)
rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt))
Vrep = @. krep*rtcb1*rt

rd = get_curve(solu, :rd)
rtca = get_curve(solu, :rtca)
rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
Vtag = @. ktag*rtca1*rd


params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "orig");

rtcb_orig = get_curve(solu_orig, :rtcb)
rt_orig = get_curve(solu_orig, :rt)
rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
Vrep_orig = @. krep*rtcb1_orig*rt_orig

rd_orig = get_curve(solu_orig, :rd)
rtca_orig = get_curve(solu_orig, :rtca)
rtca1_orig = @. (atp*rtca_orig)/(atp+(km_a*rd_orig)) 
Vtag_orig = @. ktag*rtca1_orig*rd_orig

[p1 p2]

plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB"), scatter(x=solu_orig.t, y=rtcb_orig, name="RtcB orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=rtcb1, name="RtcB1 inhib"), scatter(x=solu_orig.t, y=rtcb1_orig, name="RtcB1 orig")], Layout(xaxis=attr(range=(0, 800))))
plot([scatter(x=solu.t, y=Vrep, name="Vrep inhib"), scatter(x=solu_orig.t, y=Vrep_orig, name="Vrep orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t,y=rt, name="Rt inhib"), scatter(x=solu_orig.t, y=rt_orig, name="Rt orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=Vtag, name="Vtag inhib"), scatter(x=solu_orig.t, y=Vtag_orig, name="Vtag orig")], Layout(xaxis=attr(range=(0, 800))))



rtcb_ss=[]
rt_ss=[]
rtcb1_ss=[]
Vrep_ss=[]
rtcb_ss_orig=[]
rt_ss_orig=[]
rtcb1_ss_orig=[]
Vrep_ss_orig=[]
for i in range(0,3,length=50)
    params_inhib_new = deepcopy(params_inhib)
    params1_new = deepcopy(params1)
    params_inhib_new.kdam = i
    params1_new.kdam = i
    solu = sol(rtc_inhib_model, initial_i, tspan, params_inhib_new)
    solu_orig = sol(rtc_model, initial, tspan, params1_new)
    rtcb = get_ssval(solu, :rtcb)
    rt = get_ssval(solu, :rt)
    rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt)) 
    Vrep = @. krep*rtcb1*rt
    rtcb_orig = get_ssval(solu_orig, :rtcb)
    rt_orig = get_ssval(solu_orig, :rt)
    rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
    Vrep_orig = @. krep*rtcb1_orig*rt_orig
    push!(rtcb_ss, rtcb)
    push!(rt_ss, rt)
    push!(rtcb1_ss, rtcb1)
    push!(Vrep_ss, Vrep)
    push!(rtcb_ss_orig, rtcb_orig)
    push!(rt_ss_orig, rt_orig)
    push!(rtcb1_ss_orig, rtcb1_orig)
    push!(Vrep_ss_orig, Vrep_orig)
end

rtcb_ss

plot([scatter(x=range(0,3,length=50),y=rtcb_ss,name="RtcB inhib"), scatter(x=range(0,3,length=50),y=rtcb_ss_orig,name="RtcB orig")])
plot([scatter(x=range(0,3,length=50),y=rt_ss,name="Rt inhib"), scatter(x=range(0,3,length=50),y=rt_ss_orig,name="Rt orig")])
plot([scatter(x=range(0,3,length=50),y=rtcb1_ss,name="RtcB1 inhib"), scatter(x=range(0,3,length=50),y=rtcb1_ss_orig,name="RtcB1 orig")])
plot([scatter(x=range(0,3,length=50),y=Vrep_ss,name="Vrep inhib"), scatter(x=range(0,3,length=50),y=Vrep_ss_orig,name="Vrep orig")])



kinhib_range = 10 .^range(log10(0.0001),log10(1),length=5)
inhib_range = 10 .^range(log10(0.001),log10(10),length=5)

sweep_paramx2_new(rtc_inhib_model, :rtcb, get_ssval, :k_inhib, :inhib, kinhib_range, inhib_range, initial, params_inhib, "", "")

kdam_r = range(0,3,length=100)
res=[]
for i in kinhib_range
    p = deepcopy(params_inhib)
    p.k_inhib = i
    rtcb_ss=[]
    for kdam in kdam_r
        p.kdam = kdam
        solu = sol(rtc_inhib_model, initial, tspan, p)
        push!(rtcb_ss, get_ssval(solu, :rtcb))
    end
    push!(res, rtcb_ss)
end
plot([scatter(x=kdam_r, y=res[1], name="$(kinhib_range[1])"),scatter(x=kdam_r, y=res[2], name="$(kinhib_range[2])"),scatter(x=kdam_r, y=res[3], name="$(kinhib_range[3])"),scatter(x=kdam_r, y=res[4], name="$(kinhib_range[4])"),scatter(x=kdam_r, y=res[5], name="$(kinhib_range[5])")])

res2=[]
for i in inhib_range
    p = deepcopy(params_inhib)
    p.inhib = i
    rtcb_ss=[]
    for kdam in kdam_r
        p.kdam = kdam
        solu = sol(rtc_inhib_model, initial, tspan, p)
        push!(rtcb_ss, get_ssval(solu, :rtcb))
    end
    push!(res2, rtcb_ss)
end
res
plot([scatter(x=kdam_r, y=res2[1]),scatter(x=kdam_r, y=res2[2]),scatter(x=kdam_r, y=res2[3]),scatter(x=kdam_r, y=res2[4]),scatter(x=kdam_r, y=res2[5])])







rtcb_ss2=[]
for i in inhib_range
    p = deepcopy(params_inhib)
    p.inhib = i
    solu = sol(rtc_inhib_model, initial, tspan, p)
    push!(rtcb_ss2, get_ssval(solu, :rtcb))
end

plot(scatter(x=inhib_range, y=rtcb_ss2))









[p1_rtca p1_rtcb p2]



solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu, "", "", "RtcA inhib")

atp = 3578.9473684210525
rtcb = get_curve(solu, :rtcb)
rt = get_curve(solu, :rt)
rtcb_i = get_curve(solu, :rt)
rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt))
Vrep = @. krep*rtcb1*rt

rd = get_curve(solu, :rd)
rtca = get_curve(solu, :rtca)
rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
Vtag = @. ktag*rtca1*rd


params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "");

rtcb_orig = get_curve(solu_orig, :rtcb)
rt_orig = get_curve(solu_orig, :rt)
rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
Vrep_orig = @. krep*rtcb1_orig*rt_orig

rd_orig = get_curve(solu_orig, :rd)
rtca_orig = get_curve(solu_orig, :rtca)
rtca1_orig = @. (atp*rtca_orig)/(atp+(km_a*rd_orig)) 
Vtag_orig = @. ktag*rtca1_orig*rd_orig

[p1 p2]

plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB"), scatter(x=solu_orig.t, y=rtcb_orig, name="RtcB orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=rtcb1, name="RtcB1 inhib"), scatter(x=solu_orig.t, y=rtcb1_orig, name="RtcB1 orig")], Layout(xaxis=attr(range=(0, 800))))
plot([scatter(x=solu.t, y=Vrep, name="Vrep inhib"), scatter(x=solu_orig.t, y=Vrep_orig, name="Vrep orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t,y=rt, name="Rt inhib"), scatter(x=solu_orig.t, y=rt_orig, name="Rt orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=Vtag, name="Vtag inhib"), scatter(x=solu_orig.t, y=Vtag_orig, name="Vtag orig")], Layout(xaxis=attr(range=(0, 800))))
