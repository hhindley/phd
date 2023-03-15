using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

csv_gr = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
csv_gr = select!(csv_gr, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_atp = DataFrame(CSV.File("/home/holliehindley/phd/data/atp_curve_from_growth_model.csv"))

lam, new_df = extend_gr_curve(csv_gr)
lam[lam.< 0] .= 0 #zero(eltype(lam_colD))

# atp = QuadraticInterpolation(csv_atp."atp", csv_atp."t")
atp, new_atp = extend_atp_curve(csv_atp)

plot(scatter(x=new_atp."t", y=atp), Layout(xaxis_type="log"))

# atp and lam varied over time 
tspan = (0, lam[end])
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_atp_lam = sol(atp_gr_t!, initial, tspan, params)
p = plotly_plot_sol(solu_atp_lam, "log", "")

lam_model = get_lambda(solu_atp_lam, lam)
atp_model = get_atp(solu_atp_lam, atp)
plot(scatter(x=solu_atp_lam.t, y=lam_model), Layout(xaxis_type="log"))
plot(scatter(x=solu_atp_lam.t, y=atp_model), Layout(xaxis_type="log"))

# just lam varied over time
atp = 4000
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_lam = sol(rtc_model1!, initial, tspan, params)
p = plotly_plot_sol(solu_lam, "log", "log")

# neither atp or lam varied over time 
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr)
solu = sol(rtc_model, initial, tspan, params)
p = plotly_plot_sol(solu, "log", "log")
