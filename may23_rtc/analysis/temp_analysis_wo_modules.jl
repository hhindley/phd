using CSV, DataFrames, DifferentialEquations, LabelledArrays, StaticArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl"); include("/home/holliehindley/phd/may23_rtc/parameters/init.jl"); include("/home/holliehindley/phd/may23_rtc/functions/solving.jl");
include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); include("/home/holliehindley/phd/may23_rtc/analysis/t_param_setup.jl");

initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

# kin = 0.19
kdam = 0.4
p_none_0 = solvePlot_time(rtc_model, 0.1, 4000, 0.033, initial, "kin = 0.1", "log", "")
kdam = 0.01
p_none_01 = solvePlot_time(rtc_model, 1, 4000, 0.033, initial, "kin = 1", "log", "")
kdam = 0.1
p_none_1 = solvePlot_time(rtc_model, 1.5, 4000, 0.033, initial, "kin = 1.5", "log", "")
kdam = 1
p_none_10 = solvePlot_time(rtc_model, 2, 4000, 0.033, initial, "kin = 2", "log", "")

[p_none_0 p_none_01; p_none_1 p_none_10]
solu_none = sol(rtc_model, init, tspan, params)
init_ss = ss_init_vals(solu_none)
p_all = solvePlot_time(rtc_all_t!, kin_t, atp_t, lam_t, init_ss, "all varied over time", "log", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
L_range = 10 .^(range(1,stop=3,length=101))
# L_range = collect(0:0.1:10)
L_results = change_param(L_range, :L, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(L_range, L_results, "L", "log")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
c_range = collect(0:0.01:1)
c_results = change_param(c_range, :c, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(c_range, c_results, "c", "")
# plot_all_change_param(c_range, c_results)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
wab_range = collect(0:0.005:0.5)
wab_results = change_param(wab_range, :ω_ab, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(wab_range, wab_results, "ω_ab", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
wr_range = collect(0:0.0000001:1e-5)
wr_results = change_param(wr_range, :ω_r, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(wr_range, wr_results, "ω_r", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
atp_range = collect(0:50:5000)
atp_results = change_param(atp_range, :atp, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(atp_range, atp_results, "ATP", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
kin_range = collect(0:0.1:10)
kin_results = change_param(kin_range, :kin, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(kin_range, kin_results, "kin", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
lam_range = collect(0.1:0.01:1.1)
lam_results = change_param(lam_range, :lam, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(lam_range, lam_results, "λ", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
kdam_range = collect(0:0.01:1)
kdam_results = change_param(kdam_range, :kdam, rtc_model, initial, all_species, lam, atp, kin)
plot_change_param_sols(kdam_range, kdam_results, "kdam", "")

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :L, :c, L_range, c_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :L, :c, L_range, c_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :ω_r, :ω_ab, wr_range, wab_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :ω_r, :ω_ab, wr_range, wab_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :c, kdam_range, c_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :L, kdam_range, L_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :ω_ab, kdam_range, wab_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :ω_r, kdam_range, wr_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :atp, kdam_range, atp_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :kin, kdam_range, kin_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_a, get_ssval, :kdam, :lam, kdam_range, lam_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :kdam, :ω_r, kdam_range, wr_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :kdam, :atp, kdam_range, atp_range)
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :kdam, :lam, kdam_range, lam_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rh, get_ssval, :kdam, :c, kdam_range, c_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rh, get_ssval, :kdam, :L, kdam_range, L_range)

include("/home/holliehindley/phd/may23_rtc/parameters/params.jl");
sweep_paramx2_new(rtc_model, lam, atp, kin, :rh, get_ssval, :kdam, :kin, kdam_range, kin_range)


wr_range = collect(0:0.00000001:1e-6)

sweep_paramx2_new(rtc_model, lam, atp, kin, :rm_r, get_ssval, :atp, :ω_r, atp_range, wr_range)
