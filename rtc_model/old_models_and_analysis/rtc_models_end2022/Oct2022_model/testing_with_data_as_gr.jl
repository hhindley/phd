using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")
include("$PATHrtc_models/params_init_tspan.jl")
include("$PATHParam_inf/inf_setup.jl")

csv = DataFrame(CSV.File("$PATHdata/results_colD_grfit.csv")) # read csv to a datafram
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
plot(scatter(x=collect(0: length(csv."t")), y=csv."gr"), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT colD data"))

csv_wt = DataFrame(CSV.File("$PATHdata/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
plot(scatter(x=collect(0: length(csv_wt."t")), y=csv_wt."gr"), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT hpx data"))
  
lam_colD, new_df = extend_gr_curve(csv)
plot(scatter(x=new_df."t", y=new_df."gr"[1:68]), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT colD data"))

lam_wt, new_df_wt = extend_gr_curve(csv_wt)
plot(scatter(x=new_df_wt."t", y=new_df_wt."gr"[1:15]), Layout(xaxis_title="t", yaxis_title="growth rate", title="WT hpx data"))


tspan_wt = (0, lam_wt[end])
params_wt = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_wt] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_wt = sol(rtc_model1!, initial, tspan_wt, params_wt)
p = plotly_plot_sol(solu_wt, "log", "")

tspan_colD = (0, lam_colD[end])
params_colD = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_colD = sol(rtc_model1!, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_colD, "log", "log")


k = set_k(WT4)
OD_0 = set_OD0(WT4)
init = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0];
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, k, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :k, :lam)
solu = sol(rtc_model_OD_t!, init, tspan_colD, params)

plotly_plot_sol_OD(solu, "log")
# # p_tp = plotly_plot_sol_timepoints(solu)

# lambda_colD = get_lambda(solu, lam_colD)
# plot(scatter(x=solu_colD.t, y=lambda_colD), Layout(xaxis_type="log"))

# lambda_wt = get_lambda(solu, lam_wt)
# plot(scatter(x=solu_wt.t, y=lambda_wt), Layout(xaxis_type="log"))




# plot_dilution(solu, lam_colD)
# plot_degradation(solu)
# plot_all_variables(solu, lam_colD)

# scale_lam(csv, :rtca)
