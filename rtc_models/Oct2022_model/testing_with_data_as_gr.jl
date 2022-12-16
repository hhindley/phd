using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a dataframe
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

plot(scatter(x=collect(0: length(csv."t")), y=csv."gr"))
mean_gr = mean(csv."gr"[45:end])

df = DataFrame(t=Float64[], gr=Float64[])
for t in collect(1300:5000:1e9)
    push!(df, [t, mean_gr])
end
df

new_df = vcat(csv, df)
plot(scatter(x=new_df."t", y=new_df."gr"[1:68]))

lam_colD = QuadraticInterpolation(csv."gr",csv."t")
plot(scatter(x=csv."t", y=lam_colD))

long_lam_colD = QuadraticInterpolation(new_df."gr",new_df."t")

# csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
# gr_wt = csv_wt."gr"
# t1_wt = csv_wt."t"

# lam_wt = LinearInterpolation(gr_wt, t1_wt)
# @show lam_wt

params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

solu = sol(rtc_model1!, initial, tspan, params)

p = plotly_plot_sol(solu, "log")

p_tp = plotly_plot_sol_timepoints(solu)

lambda = get_lambda(solu, long_lam_colD)
plot(scatter(x=solu.t, y=lambda), Layout(xaxis_type="log"))

plot_sol_and_lam(solu, lam_colD)

rtca = get_curve(solu, :rtca)
rtca[83:end]
# lambda_order = sort(lambda)
# perm = sortperm(lambda)
# rtca_sort = (res[:rtca][1])[perm]
# plot(scatter(x=lambda_order, y=rtca_sort))





plot_dilution(solu, lam_colD)
plot_degradation(solu)
plot_all_variables(solu, lam_colD)

scale_lam(csv, :rtca)

ω_ab_range = collect(0:0.0001:0.00125)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p = sweep_paramx2(rtc_model1!, params, :rtca, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p1 = sweep_paramx2(rtc_model1!, params, :rtcb, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p2 = sweep_paramx2(rtc_model1!, params, :rtcr, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p3 = sweep_paramx2(rtc_model1!, params, :rm_a, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p4 = sweep_paramx2(rtc_model1!, params, :rm_b, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p5 = sweep_paramx2(rtc_model1!, params, :rm_r, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p6 = sweep_paramx2(rtc_model1!, params, :rh, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p7 = sweep_paramx2(rtc_model1!, params, :rd, get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p8 = sweep_paramx2(rtc_model1!, params, :rt, get_ssval, :ω_r, :ω_ab, ω_ab_range)




params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p_std = sweep_paramx2(rtc_model1!, params, :rtca, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p1_std = sweep_paramx2(rtc_model1!, params, :rtcb, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p2_std = sweep_paramx2(rtc_model1!, params, :rtcr, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p3_std = sweep_paramx2(rtc_model1!, params, :rm_a, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p4_std = sweep_paramx2(rtc_model1!, params, :rm_b, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p5_std = sweep_paramx2(rtc_model1!, params, :rm_r, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p6_std = sweep_paramx2(rtc_model1!, params, :rh, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p7_std = sweep_paramx2(rtc_model1!, params, :rd, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, long_lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
p8_std = sweep_paramx2(rtc_model1!, params, :rt, check_get_ssval, :ω_r, :ω_ab, ω_ab_range)





