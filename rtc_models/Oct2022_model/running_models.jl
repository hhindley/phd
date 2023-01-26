using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
  
lam_colD, new_df = extend_gr_curve(csv)
lam_wt, new_df_wt = extend_gr_curve(csv_wt)

tspan_wt = (0, lam_wt[end])
params_wt = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_wt] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_wt = sol(rtc_model1!, initial, tspan_wt, params_wt)
p = plotly_plot_sol(solu_wt, "log", "")

tspan_colD = (0, lam_colD[end])
params_colD = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_colD = sol(rtc_model1!, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_colD, "log", "log")

params_kin_kdam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, ktag, kdeg, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :ktag, :kdeg, :atp, :na, :nb, :nr, :lam)


solu_kin_kdam = sol(rtc_model_kin_kdam_t!, initial, tspan_colD, params_kin_kdam)
p = plotly_plot_sol(solu_kin_kdam, "log", "log")
res = get_all_curves(solu_kin_kdam, all_species)
lam = get_lambda(solu_kin_kdam, lam_colD)
kdam = @. (kin/res[:rh][1]) - lam
kin1 = @. lam*res[:rh][1]/g_max
plot(scatter(x=solu_kin_kdam.t, y=kdam), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam.t, y=kin1), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam.t, y=lam), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam.t, y=res[:rh][1]), Layout(xaxis_type="log"))


solu_kin_kdam_tlr = sol(rtc_model_kin_kdam_tlr!, initial, tspan_colD, params_kin_kdam)
p_tlr = plotly_plot_sol(solu_kin_kdam_tlr, "log", "log")
res_tlr = get_all_curves(solu_kin_kdam_tlr, all_species)
lam_tlr = get_lambda(solu_kin_kdam_tlr, lam_colD)
kdam_tlr = @. (kin/res_tlr[:rh][1]) - lam_tlr
kin_tlr = @. lam_tlr*res_tlr[:rh][1]/((g_max*atp/(θtlr+atp)))
plot(scatter(x=solu_kin_kdam_tlr.t, y=kdam_tlr), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam_tlr.t, y=kin_tlr), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam_tlr.t, y=lam_tlr), Layout(xaxis_type="log"))
plot(scatter(x=solu_kin_kdam_tlr.t, y=res_tlr[:rh][1]), Layout(xaxis_type="log"))


g_max*atp/(θtlr+atp)
g_max