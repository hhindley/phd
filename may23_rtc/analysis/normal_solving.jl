push!(LOAD_PATH, "/home/holliehindley/phd/may23_rtc/Modules")
using CSV, DataFrames, DifferentialEquations, StaticArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics

using ModelsModule, SolvingModule, PlottingModule

include("/home/holliehindley/phd/may23_rtc/parameters/init.jl")
include("/home/holliehindley/phd/may23_rtc/parameters/params.jl")


# init = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];
# params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)



# nothing varied over time 
p_none = solvePlot_time(rtc_model, L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, 0.054, 4000, na, nb, nr, 0.033, "nothing varied over time", "")

open("./rtc_model_figs/nothing_varied.html", "w") do io
    PlotlyBase.to_html(io, p_none.plot)
end


# just lam varied over time
p_justlam = solvePlot_time(rtc_model_lam_t!, lam_t, 4000, 0.054, "just λ varied over time", "")

open("./rtc_model_figs/just_lam.html", "w") do io
    PlotlyBase.to_html(io, p_justlam.plot)
end


# just vary kin over time 
p_justkin = solvePlot_time(rtc_model_kin, 0.033, 4000, kin_t, "just kin varied over time", "")

open("./rtc_model_figs/justkin.html", "w") do io
    PlotlyBase.to_html(io, p_justkin.plot)
end


# just atp varied over time 
p_justatp = solvePlot_time(atp_t!, 0.033, atp_t, 0.054, "just ATP varied over time", "")

open("./rtc_model_figs/justatp.html", "w") do io
    PlotlyBase.to_html(io, p_justatp.plot)
end


# atp and lam varied over time 
p_atplam = solvePlot_time(atp_gr_t!, lam_t, atp_t, 0.054, "λ and ATP varied over time", "")

open("./rtc_model_figs/atplam.html", "w") do io
    PlotlyBase.to_html(io, p_atplam.plot)
end


# atp and kin varied over time 
p_atpkin = solvePlot_time(atp_kin_t!, 0.033, atp_t, kin_t, "ATP and kin varied over time", "")

open("./rtc_model_figs/atpkin.html", "w") do io
    PlotlyBase.to_html(io, p_atpkin.plot)
end


# kin and lam varied over time 
p_lamkin = solvePlot_time(lam_kin_t, lam_t, 4000, kin_t, "kin and λ varied over time", "")

open("./rtc_model_figs/kinlam.html", "w") do io
    PlotlyBase.to_html(io, p_lamkin.plot)
end


#all 
tspan = (0,1440)
p_all = solvePlot_time(rtc_all_t!, lam_t, atp_t, kin_t, "λ, ATP and kin varied over time", "")

savefig(p_all,"p_all.svg")


open("./rtc_model_figs/all.html", "w") do io
    PlotlyBase.to_html(io, p_all.plot)
end


sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_a, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))
sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_r, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))

# lam_model = get_lambda(solu_atplam, lam_t)
# atp_model = get_atp(solu_atplam, atp_t)
# plot(scatter(x=solu_atp_lam.t, y=lam_model), Layout(xaxis_type="log"))
# plot(scatter(x=solu_atp_lam.t, y=atp_model), Layout(xaxis_type="log"))
