using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


ω_ab_range = collect(range(0, 0.1, length=10))
ω_r_range = collect(range(0, 0.1, length=10))




kdam_range = collect(0:0.001:0.2)
results_kdam = change_param(kdam_range, :kdam, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(kdam_range, results_kdam, "kdam")

ωab_range = collect(0:0.01:1)
results_ab = change_param(ωab_range, :ω_ab, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(ωab_range, results_ab, "ω_ab")

ωr_range = collect(0:0.01:1)
results_r = change_param(ωr_range, :ω_r, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(ωr_range, results_r, "ω_r")

L_range = collect(0:0.01:1)
results_L = change_param(L_range, :L, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(L_range, results_L, "L")

c_range = collect(0:0.01:1)
results_c = change_param(c_range, :c, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(c_range, results_c, "c")

kr_range = collect(0:0.01:1)
results_kr = change_param(kr_range, :kr, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(kr_range, results_kr, "kr")

vmax_range = collect(0:0.01:1)
results_vmax = change_param(vmax_range, :Vmax_init, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(vmax_range, results_vmax, "Vmax_init")

km_range = collect(0:0.01:1)
results_km = change_param(km_range, :Km_init, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(km_range, results_km, "km_init")

tscr_range = collect(0:0.01:1)
results_tscr = change_param(tscr_range, :θtscr, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(tscr_range, results_tscr, "θtscr")

tlr_range = collect(0:0.01:1)
results_tlr = change_param(tlr_range, :θtlr, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(tlr_range, results_tlr, "θtlr")

kb_range = collect(0:0.01:1)
results_kb = change_param(kb_range, :k_b, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(kb_range, results_kb, "kb")

d_range = collect(0:0.01:1)
results_d = change_param(d_range, :d, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(d_range, results_d, "d")

krep_range = collect(0:0.01:1)
results_krep = change_param(krep_range, :krep, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(krep_range, results_krep, "krep")

ktag_range = collect(0:0.01:1)
results_ktag = change_param(ktag_range, :ktag, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(ktag_range, results_ktag, "ktag")

atp_range = collect(0:0.01:1)
results_atp = change_param(atp_range, :atp, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(atp_range, results_atp, "atp")

km_a_range = collect(0:0.01:1)
results_km_a = change_param(km_a_range, :km_a, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(km_a_range, results_km_a, "km_a")

km_b_range = collect(0:0.01:1)
results_km_b = change_param(km_b_range, :km_b, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(km_b_range, results_km_b, "km_b")

gmax_range = collect(0:0.01:1)
results_gmax = change_param(gmax_range, :g_max, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(gmax_range, results_gmax, "gmax")

grc_range = collect(0:0.01:1)
results_grc = change_param(grc_range, :gr_c, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(grc_range, results_grc, "grc")

kdeg_range = collect(0:0.01:1)
results_kdeg = change_param(kdeg_range, :kdeg, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(kdeg_range, results_kdeg, "kdeg")

kin_range = collect(0:0.01:1)
results_kin = change_param(kin_range, :kin, rtc_model1!, initial, all_species, lam_colD)
plot_change_param_sols(kin_range, results_kin, "kin")







results_kdam_OD = change_param_OD(kdam_range, :kdam, rtc_model_OD_t!, all_species_OD, lam_colD, WT4)
plot_change_param_sols_OD(kdam_range, results_kdam_OD, "kdam")

results_ab_OD = change_param_OD(ωab_range, :ω_ab, rtc_model_OD_t!, all_species_OD, lam_colD, WT4)
plot_change_param_sols_OD(ωab_range, results_ab_OD, "ω_ab")

results_r_OD = change_param_OD(ωr_range, :ω_r, rtc_model_OD_t!, all_species_OD, lam_colD, WT4)#
plot_change_param_sols_OD(ωr_range, results_r_OD, "ω_r")





pc = sweep_paramx2(rtc_model1!, lam_colD, :rtca, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p1c = sweep_paramx2(rtc_model1!, lam_colD, :rtcb, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p2c = sweep_paramx2(rtc_model1!, lam_colD, :rtcr, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p3c = sweep_paramx2(rtc_model1!, lam_colD, :rm_a, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p4c = sweep_paramx2(rtc_model1!, lam_colD, :rm_b, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p5c = sweep_paramx2(rtc_model1!, lam_colD, :rm_r, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p6c = sweep_paramx2(rtc_model1!, lam_colD, :rh, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p7c = sweep_paramx2(rtc_model1!, lam_colD, :rd, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p8c = sweep_paramx2(rtc_model1!, lam_colD, :rt, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)

ph = sweep_paramx2(rtc_model1!, lam_wt, :rtca, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p1h = sweep_paramx2(rtc_model1!, lam_wt, :rtcb, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p2h = sweep_paramx2(rtc_model1!, lam_wt, :rtcr, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p3h = sweep_paramx2(rtc_model1!, lam_wt, :rm_a, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p4h = sweep_paramx2(rtc_model1!, lam_wt, :rm_b, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p5h = sweep_paramx2(rtc_model1!, lam_wt, :rm_r, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p6h = sweep_paramx2(rtc_model1!, lam_wt, :rh, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p7h = sweep_paramx2(rtc_model1!, lam_wt, :rd, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p8h = sweep_paramx2(rtc_model1!, lam_wt, :rt, get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)


# std to check steady state
pc_std = sweep_paramx2(rtc_model1!, lam_colD, :rtca, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p1c_std = sweep_paramx2(rtc_model1!, lam_colD, :rtcb, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p2c_std = sweep_paramx2(rtc_model1!, lam_colD, :rtcr, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p3c_std = sweep_paramx2(rtc_model1!, lam_colD, :rm_a, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p4c_std = sweep_paramx2(rtc_model1!, lam_colD, :rm_b, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p5c_std = sweep_paramx2(rtc_model1!, lam_colD, :rm_r, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p6c_std = sweep_paramx2(rtc_model1!, lam_colD, :rh, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p7c_std = sweep_paramx2(rtc_model1!, lam_colD, :rd, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p8c_std = sweep_paramx2(rtc_model1!, lam_colD, :rt, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)

ph_std = sweep_paramx2(rtc_model1!, lam_wt, :rtca, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p1h_std = sweep_paramx2(rtc_model1!, lam_wt, :rtcb, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p2h_std = sweep_paramx2(rtc_model1!, lam_wt, :rtcr, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p3h_std = sweep_paramx2(rtc_model1!, lam_wt, :rm_a, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p4h_std = sweep_paramx2(rtc_model1!, lam_wt, :rm_b, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p5h_std = sweep_paramx2(rtc_model1!, lam_wt, :rm_r, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p6h_std = sweep_paramx2(rtc_model1!, lam_wt, :rh, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p7h_std = sweep_paramx2(rtc_model1!, lam_wt, :rd, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)
p8h_std = sweep_paramx2(rtc_model1!, lam_wt, :rt, check_get_ssval, :ω_r, :ω_ab, ω_r_range, ω_ab_range)

# 3d plots
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

p3d1 = sweep_paramx3(rtc_model1!, lam_wt, :rtca, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d2 = sweep_paramx3(rtc_model1!, lam_wt, :rtcb, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d3 = sweep_paramx3(rtc_model1!, lam_wt, :rtcr, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d4 = sweep_paramx3(rtc_model1!, lam_wt, :rm_a, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d5 = sweep_paramx3(rtc_model1!, lam_wt, :rm_b, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d6 = sweep_paramx3(rtc_model1!, lam_wt, :rm_r, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d7 = sweep_paramx3(rtc_model1!, lam_wt, :rh, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d8 = sweep_paramx3(rtc_model1!, lam_wt, :rd, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)
p3d9 = sweep_paramx3(rtc_model1!, lam_wt, :rt, get_ssval, :ω_r, :ω_ab, :kdam, ω_r_range)



