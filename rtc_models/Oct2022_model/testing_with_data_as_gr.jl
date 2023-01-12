using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
plot(scatter(x=collect(0: length(csv."t")), y=csv."gr"), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT colD data"))

csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
plot(scatter(x=collect(0: length(csv_wt."t")), y=csv_wt."gr"), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT hpx data"))

  
lam_colD, new_df = extend_gr_curve(csv)
plot(scatter(x=new_df."t", y=new_df."gr"[1:68]), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT colD data"))

lam_wt, new_df_wt = extend_gr_curve(csv_wt)
plot(scatter(x=new_df_wt."t", y=new_df_wt."gr"[1:15]), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT hpx data"))


tspan_wt = (0, lam_wt[end])
params_wt = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_wt] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_wt = sol(rtc_model1!, initial, tspan_wt, params_wt)
p = plotly_plot_sol(solu_wt, "log")

tspan_colD = (0, lam_colD[end])
params_colD = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_colD = sol(rtc_model1!, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_colD, "log")


# # p_tp = plotly_plot_sol_timepoints(solu)

# lambda_colD = get_lambda(solu, lam_colD)
# plot(scatter(x=solu_colD.t, y=lambda_colD), Layout(xaxis_type="log"))

# lambda_wt = get_lambda(solu, lam_wt)
# plot(scatter(x=solu_wt.t, y=lambda_wt), Layout(xaxis_type="log"))




# plot_dilution(solu, lam_colD)
# plot_degradation(solu)
# plot_all_variables(solu, lam_colD)

# scale_lam(csv, :rtca)

ω_ab_range = collect(range(0, 1, length=10))
ω_r_range = collect(range(0, 1, length=10))

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



