using Plots, DifferentialEquations, DataFrames, Statistics, Measures
include("fullmodel_functions.jl")
include("fullmodel_params.jl")
include("models_new.jl") # use model RTC_SYSTEM2!

include("prototype_model_full.jl") # model here is the same as RTC_SYSTEM2!


# plotlyjs()

sol = simple_solve!(RTC_SYSTEM2!, init_rtc, tspan, params_new)
plot(sol[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_rtc, palette=:tab10)


sol1 = simple_solve!(rtc_model!, init_rtc, tspan, params_new)
plot(sol1[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_rtc, palette=:tab10)

# splitting species and looking at solution results
rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_reduced_model(species_rtc, sol)
plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_all)

plot(sol.t, rm_a, linecolor=:green3, annotations=([5e2],[rm_a[end]+72], "mRNA-RtcA/B"), annotationcolor=:green3, annotationfontsize=10, xaxis=(:log10, (1,1e3)), 
    yaxis=(:log10, (1,Inf)), title="Model solution", legend=false, xlabel="time", ylabel="amount of species")
plot!(sol.t, rtca, linecolor=:green, annotations=([6.9e2],[rtca[end]+10000], "RtcA/B"), annotationcolor=:green, annotationfontsize=10)
plot!(sol.t, rm_r, linecolor=:orange, annotations=([5.3e2],[rm_r[end]+8], "mRNA-RtcR"), annotationcolor=:orange, annotationfontsize=10)
plot!(sol.t, rtcr, linecolor=:darkorange2, annotations=([7.7e2],[rtcr[end]+1000], "RtcR"), annotationcolor=:darkorange2, annotationfontsize=10)
plot!(sol.t, ribo_h, linecolor=:dodgerblue2, annotations=([8.5e2],[ribo_h[end]+3], "Rh"), annotationcolor=:dodgerblue2, annotationfontsize=10)

plot(sol.t, ribo_h, linecolor=:dodgerblue2, annotations=([8e2],[ribo_h[end]+0.3], "Rh"), annotationcolor=:dodgerblue2, annotationfontsize=10, xaxis=(:log10, (1,1e3)), title="Model solution - ribosomes", legend=false, xlabel="time", ylabel="amount of species")
plot!(sol.t, ribo_d, linecolor=:blue, annotations=([8e2],[ribo_d[end]+1], "Rd"), annotationcolor=:blue, annotationfontsize=10)
plot!(sol.t, ribo_t, linecolor=:deepskyblue, annotations=([8e2],[ribo_t[end]+0.3], "Rt"), annotationcolor=:deepskyblue, annotationfontsize=10)

plot(sol.t, [rdrtca, rtrtcb], xaxis=(:log10, (1,Inf)))
a = sol.t[end]-50000
rm_a[end]


plot(sol.t, [rm_a, rm_b, rm_r], xaxis=(:log10, (1,Inf)))
plot(sol.t, [rm_r], xaxis=(:log10, (1,Inf)))
plot(sol.t, [rtca, rtcb, rtcr], xaxis=(:log10, (1,Inf)))
plot(sol.t, [ribo_h, ribo_d, ribo_t], xaxis=(:log10, (1,Inf)), labels=["Rh" "Rd" "Rt"])
plot(sol.t, [ribo_h], xaxis=(:log10, (1,Inf)))
plot(sol.t, [ribo_d], xaxis=(:log10, (1,Inf)))
plot(sol.t, [ribo_t], xaxis=(:log10, (1,Inf)))

p_sol = plot(sol.t, [ribo_d+ribo_h+ribo_t], xaxis=(:log10, (1,Inf)), margin=16mm)

p_rtc = plot(sol.t, xlabel="t", [rtca, rtcb, rtcr], labels=["RtcA" "RtcB" "RtcR"], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
p_mrna = plot(sol.t, xlabel="t", [rm_a, rm_b, rm_r], labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], xaxis=(:log10, (1,Inf)))
p_ribo = plot(sol.t, xlabel="t", [ribo_d, ribo_h, ribo_t], labels=["Rd" "Rh" "Rt"], xaxis=(:log10, (1,Inf)))#, title="Ribosomes")
p_int = plot(sol.t, xlabel="t", [rdrtca, rtrtcb], title="Intermediates", labels=["RdRtcA" "RtRtcB"], xaxis=(:log10, (1,Inf)))

plot(p_sol, p_rtc, p_mrna, p_ribo, p_int, size=(1000,1000))

plot(sol.t, ribo_d, xlabel="t", title="Damaged ribosomes", legend=false, xaxis=(:log10, (1,Inf)))
plot(sol.t, ribo_t, xlabel="t", title="Tagged ribosomes", legend=false, xaxis=(:log10, (1,Inf)))
plot(sol.t, ribo_h, xlabel="t", title="Healthy ribosomes", legend=false, xaxis=(:log10, (1,Inf)))


# varying parameters
rd_range = range(0, 10, length=6)
include("fullmodel_params.jl")
rd_p1, rd_p2 = logplots_reduced(rd_range, init_rtc, 8, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "Rd")
rd_p1
# rd_p2

rh_range = range(0, 10, length=6)
include("fullmodel_params.jl")
rh_p1, rh_p2 = logplots_reduced(rh_range, init_rtc, 7, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "Rh")
rh_p1
# rh_p2

atp_range = range(0, 100, length=6) 
include("fullmodel_params.jl")
atp_p1, atp_p2 = logplots_reduced(atp_range, params_new, 5, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "ATP")
atp_p1
# atp_p2

kr_range = range(0.01, 2, length=6)
include("fullmodel_params.jl") 
kr_p1, kr_p2 = logplots_reduced(kr_range, params_new, 1, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "kr")
kr_p1
# kr_p2

L_range = range(0, 100, length=6) 
include("fullmodel_params.jl")
L_p1, L_p2 = logplots_reduced(L_range, params_new, 2, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "L")
L_p1
# L_p2

c_range = range(0, 100, length=6)
include("fullmodel_params.jl")
c_p1, c_p2 = logplots_reduced(c_range, params_new, 3, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "c")
c_p1
# c_p2

sigma_range = range(0, 0.2, length=6)
include("fullmodel_params.jl")
sigma_p1, sigma_p2 = logplots_reduced(sigma_range, params_new, 4, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "sigma")
sigma_p1
# sigma_p2

k1_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1_p1, k1_p2 = logplots_reduced(k1_range, params_new, 11, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k1")
k1_p1
# k1_p2

k2_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2_p1, k2_p2 = logplots_reduced(k2_range, params_new, 12, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k2")
k2_p1
# k2_p2

k3_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k3_p1, k3_p2 = logplots_reduced(k3_range, params_new, 13, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k3")
k3_p1
# k3_p2

k4_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k4_p1, k4_p2 = logplots_reduced(k4_range, params_new, 14, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k4")
k4_p1
# k4_p2

k5_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k5_p1, k5_p2 = logplots_reduced(k5_range, params_new, 15, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k5")
k5_p1
# k5_p2

wrtc_range = range(0, 10, length=6)
include("fullmodel_params.jl")
wrtc_p1, wrtc_p2 = logplots_reduced(wrtc_range, params_new, 6, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "wrtc")
wrtc_p1
# wrtc_p2

thrtc_range = range(0, 20, length=6)
include("fullmodel_params.jl")
thrtc_p1, thrtc_p2 = logplots_reduced(thrtc_range, params_new, 7, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "thrtc")
thrtc_p1
# thrtc_p2

max_range = range(0, 10, length=6)
include("fullmodel_params.jl")
max_p1, max_p2 = logplots_reduced(max_range, params_new, 8, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "max")
max_p1
# max_p2

thr_range = range(0, 10, length=6)
include("fullmodel_params.jl")
thr_p1, thr_p2 = logplots_reduced(thr_range, params_new, 9, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "thr")
thr_p1
# thr_p2

# d_range = range(0, 5, length=6)
# include("fullmodel_params.jl")
# d_p1, d_p2 = logplots_reduced(d_range, params_new, 10, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "d")
# d_p1
# d_p2

k1a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1a_p1, k1a_p2 = logplots_reduced(k1a_range, params_new, 16, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k1a")
k1a_p1
# k1a_p2

k2a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2a_p1, k2a_p2 = logplots_reduced(k2a_range, params_new, 17, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k2a")
k2a_p1
# k2a_p2

k3a_range = range(0, 0.1, length=6)
include("fullmodel_params.jl")
k3a_p1, k3a_p2 = logplots_reduced(k3a_range, params_new, 18, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k3a")
k3a_p1
# k3a_p2

kab_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kab_p1, kab_p2 = logplots_reduced(kab_range, params_new, 19, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "kab")
kab_p1
# kab_p2

kbb_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kbb_p1, kbb_p2 = logplots_reduced(kbb_range, params_new, 20, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "kbb")
kbb_p1
# kbb_p2

kcb_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kcb_p1, kcb_p2 = logplots_reduced(kcb_range, params_new, 21, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "kcb")
kcb_p1
# kcb_p2

k_range = range(0, 1, length=6)
include("fullmodel_params.jl")
k_p1, k_p2 = logplots_reduced(k_range, params_new, 22, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "k")
k_p1
k_p2

rtca_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtca_p1, rtca_p2 = logplots_reduced(rtca_range, init_rtc, 2, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "rtca")
# rtca_p1
rtca_p2

rtcb_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtcb_p1, rtcb_p2 = logplots_reduced(rtcb_range, init_rtc, 4, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "rtcb")
rtcb_p1
# rtcb_p2

rt_range = range(0, 10, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in rtc model to ode and see if that works (probs will) 
include("fullmodel_params.jl")
rt_p1, rt_p2 = logplots_reduced(rt_range, init_rtc, 9, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, labels_rtc, "Rt")
rt_p1
# rt_p2

# plotting parameter against steady states 
atp_range_long = range(0,6000, length=10)
include("fullmodel_params.jl")
plot(atp_range, rm_a_sols, legend=false, annotations=([atp_range[end]],[rm_a_sols[end]], "mRNA-RtcA"))#text.("Rt", :bottom))
plot(atp_range, rtca_sols, label="RtcA")
plot(atp_range, rm_b_sols, label="mRNA-RtcB")
plot(atp_range, rtcb_sols, label="RtcB")
plot(atp_range, rm_r_sols, label="mRNA-RtcR")
plot(atp_range, rtcr_sols, label="RtcR")
plot(atp_range, rdrtca_sols, label="RdRtcA")
plot(atp_range, rtrtcb_sols, label="RtRtcB")
plot(atp_range, ribo_h_sols, label="Rh")
plot(atp_range, ribo_d_sols, label="Rd")
plot(atp_range, ribo_t_sols, label="Rt")

plot(atp_range, [ribo_h_sols, ribo_d_sols, ribo_t_sols])#, xlim=(0,6), ylim=(0,6))

atp_range_long = range(0,10000, length=1000)
include("fullmodel_params.jl")
rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = param_vs_ss(atp_range_long, params_new, 5, RTC_SYSTEM2!, init_rtc, params_new, species_rtc)
plot(atp_range_long, [rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols], yaxis=(:log10, (1,Inf)), labels=labels_rtc)

print(length(rtca_sols[1:end]))
plot(atp_range_long, rtca_sols[1:end], yaxis=:log10, xaxis=(:log10, (1, Inf)))
plot(atp_range_long, rtcr_sols)
plot(atp_range_long[atp_range_long.>0],[rtca_sols[rtca_sols.>0], rtcr_sols[rtcr_sols.>0]], yaxis=:log10, xaxis=:log10)

function plot_inset()
    p = plot(atp_range_long, [rtca_sols], title="Dependence on ATP", xlabel="ATP", 
    ylabel="amount of species at steady state", legend=false, linecolor=:green, 
    annotations=([atp_range_long[end]-2],[rtca_sols[end]], "RtcA/B"), annotationcolor=:green,
    xaxis=(:log10, (1, Inf)), yaxis=(:log10, (1, Inf)), 
    grid=false)
    plot!(atp_range_long, rtcr_sols, linecolor=:darkorange2, annotations=([atp_range_long[end]],
    [rtcr_sols[end]+5000], "RtcR"), annotationcolor=:darkorange2)
    # lens!([0,5],[0,1000], inset = (1, bbox(0.2, 0.1, 0.3, 0.3)), lc=:false)
    # plot!(p[2], atp_range, [rtca_sols], xlims=(0,5), ylims=(0,1000), legend=false, linecolor=:green, xlabel="ATP", ylabel="amount of species", xtickfontsize=7, ytickfontsize=7, xlabelfontsize=7, ylabelfontsize=7, annotations=([atp_range[15]],[rtca_sols[end]-102800], "RtcA/B"), annotationcolor=:green, annotationfontsize=7)
    # plot!(p[2], atp_range, [rtcr_sols], xlims=(0,5), ylims=(0,1000), legend=false, linecolor=:darkorange2, annotations=([atp_range[25]],[rtcr_sols[end]-3400], "RtcR"), annotationcolor=:darkorange2, annotationfontsize=7)
end

p1 = plot_inset()


# all param analysis 
atp_range_long = range(0,100, length=20)
include("fullmodel_params.jl")
atp_plot = plot_long_range_results(atp_range_long, params_new, 5, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "ATP")

kr_range_long = range(0.01,20, length=20)
include("fullmodel_params.jl")
kr_plot = plot_long_range_results(kr_range_long, params_new, 1, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "kr")

L_range_long = range(0,100, length=20)
include("fullmodel_params.jl")
L_plot = plot_long_range_results(L_range_long, params_new, 2, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "L")

c_range_long = range(0,100, length=20)
include("fullmodel_params.jl")
c_plot = plot_long_range_results(c_range_long, params_new, 3, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "c")

sigma_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
sigma_plot = plot_long_range_results(sigma_range_long, params_new, 4, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "sigma")

k1_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k1_plot = plot_long_range_results(k1_range_long, params_new, 11, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k1")

k2_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k2_plot = plot_long_range_results(k2_range_long, params_new, 12, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k2")

k3_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k3_plot = plot_long_range_results(k3_range_long, params_new, 13, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k3")

k4_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k4_plot = plot_long_range_results(k4_range_long, params_new, 14, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k4")

k5_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k5_plot = plot_long_range_results(k5_range_long, params_new, 15, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k5")

wrtc_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
wrtc_plot = plot_long_range_results(wrtc_range_long, params_new, 6, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "wrtc")

theta_range_long = range(0,20, length=20)
include("fullmodel_params.jl")
theta_plot = plot_long_range_results(theta_range_long, params_new, 7, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "theta")

max_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
max_plot = plot_long_range_results(max_range_long, params_new, 8, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "max")

thr_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
thr_plot = plot_long_range_results(thr_range_long, params_new, 9, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "thr")

k1a_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
k1a_plot = plot_long_range_results(k1a_range_long, params_new, 16, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k1a")

k2a_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
k2a_plot = plot_long_range_results(k2a_range_long, params_new, 17, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k2a")

k3a_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k3a_plot = plot_long_range_results(k3a_range_long, params_new, 18, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k3a")

kab_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
kab_plot = plot_long_range_results(kab_range_long, params_new, 19, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "kab")

kbb_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
kbb_plot = plot_long_range_results(kbb_range_long, params_new, 20, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "kbb")

kcb_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
kcb_plot = plot_long_range_results(kcb_range_long, params_new, 21, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "kcb")

k_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k_plot = plot_long_range_results(k_range_long, params_new, 22, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "k")














rd_range_long = range(0,10, length=11)
include("fullmodel_params.jl")
rd_plot = plot_long_range_results(rd_range_long, init_rtc, 8, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "Rd")

rh_range_long = range(0,10, length=11)
include("fullmodel_params.jl")
rh_plot = plot_long_range_results(rh_range_long, init_rtc, 7, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "Rh")

rt_range_long = range(0,10, length=11)
include("fullmodel_params.jl")
rt_plot = plot_long_range_results(rt_range_long, init_rtc, 9, RTC_SYSTEM2!, init_rtc, params_new, species_rtc, "Rt")
