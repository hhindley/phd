using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

# load csv for growth rate and atp curves 
csv_gr = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
csv_gr = select!(csv_gr, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_atp = DataFrame(CSV.File("/home/holliehindley/phd/data/atp_for_rtcmodel.csv"))

# set atp and lam curves from files 
lam_t, new_df = extend_gr_curve(csv_gr)
lam_t[lam_t.< 0] .= 0 #zero(eltype(lam_colD))
p_lam = plot(scatter(x=new_df."t", y=lam_t), Layout(xaxis_type="log"))


lam_t = QuadraticInterpolation(csv_atp."lam",csv_atp."t")
atp_t = QuadraticInterpolation(csv_atp."atp",csv_atp."t")
p_atp = plot(scatter(x=csv_atp."t", y=atp_t), Layout(xaxis_type="log"))
p_lam = plot(scatter(x=csv_atp."t", y=lam_t), Layout(xaxis_type="log"))
[p_lam;p_atp]

# set kin curve  
rh = 11.29
kin_model = lam_t*rh/g_max
plot(scatter(x=new_df."t", y=kin_model), Layout(xaxis_type="log"))
kin_t = QuadraticInterpolation(kin_model, csv_atp."t") 
plot(scatter(x=new_df."t", y=kin_t), Layout(xaxis_type="log"))

# nothing varied over time 
tspan = (0, lam_t[end])
lam = 0.033; atp = 4000; kin = 0.054;
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu = sol(rtc_model, initial, tspan, params)
p = plotly_plot_sol(solu, "log", "", "nothing varied over time")


# just lam varied over time
atp = 4000; kin = 0.054;
params_lam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_lam = sol(rtc_model1!, initial, tspan, params_lam)
p = plotly_plot_sol(solu_lam, "log", "", "λ varied over time")
open("./model_output.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
p_all = plot_all_variables(solu_lam, lam_t)
open("./all_model_output.html", "w") do io
    PlotlyBase.to_html(io, p_all.plot)
end




rh = get_curve(solu_lam, :rh)
rd = get_curve(solu_lam, :rd); rt = get_curve(solu_lam, :rt)
rtot = @. rh+rd+rt
plot(scatter(x=solu_lam.t,y=rtot), Layout(xaxis_type="log"))
# open("./justlam.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end

# just vary kin over time 
atp = 4000; lam = 0.033;
params_kin = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_kin = sol(rtc_model_kin, initial, tspan, params_kin)

plotly_plot_sol(solu_kin, "log", "", "k_in varied over time")
rh = get_curve(solu_kin, :rh); rd = get_curve(solu_kin, :rd); rt = get_curve(solu_kin, :rt)
rtot = @. rh+rd+rt
plot(scatter(x=solu_lam.t,y=rtot), Layout(xaxis_type="log"))

# just atp varied over time 
lam = 0.033; kin = 0.054;
params_atp = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp_t, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_atp = sol(atp_t!, initial, tspan, params_atp)
plotly_plot_sol(solu_atp, "log", "", "atp varied over time")


# atp and lam varied over time 
kin = 0.054;
params_atp_lam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp_t, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_atp_lam = sol(atp_gr_t!, initial, tspan, params_atp_lam)
p = plotly_plot_sol(solu_atp_lam, "log", "", "atp and λ varied over time")

lam_model = get_lambda(solu_atp_lam, lam_t)
atp_model = get_atp(solu_atp_lam, atp_t)
plot(scatter(x=solu_atp_lam.t, y=lam_model), Layout(xaxis_type="log"))
plot(scatter(x=solu_atp_lam.t, y=atp_model), Layout(xaxis_type="log"))


# atp and kin varied over time 
lam = 0.033;
params_atp_kin = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp_t, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_atp_kin = sol(atp_kin_t!, initial, tspan, params_atp_kin)
plotly_plot_sol(solu_atp_kin, "log", "", "atp and k_in varied over time")

# kin and lam varied over time 
atp = 4000;
params_lam_kin = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_lam_kin = sol(lam_kin_t, initial, tspan, params_lam_kin)
p1 = plotly_plot_sol(solu_lam_kin, "log", "", "k_in and λ varied over time")
open("./lam_kin.html", "w") do io
    PlotlyBase.to_html(io, p1.plot)
end

rh = get_curve(solu_lam_kin, :rh); rd = get_curve(solu_lam_kin, :rd); rt = get_curve(solu_lam_kin, :rt)
rtot = @. rh+rd+rt
plot(scatter(x=solu_lam.t,y=rtot), Layout(xaxis_type="log"))