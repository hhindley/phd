using PlotlyJS, DifferentialEquations, DataFrames, Statistics, Measures
include("fullmodel_functions.jl")
include("fullmodel_params.jl")
include("models_new.jl")
# plotlyjs()

sol = simple_solve!(new_ODEs!, init_new, tspan, params_new)
plot(sol[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_new, palette=:tab10)


# splitting species and looking at solution results
rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_full_model(species_full, sol)
plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_all)

p_sol = plot(sol.t, [ribo_d+ribo_h+ribo_t+rdrtca+rtrtcb], xaxis=(:log10, (1,Inf)), margin=16mm)

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
rd_p1, rd_p2 = logplots(rd_range, init_new, 10, new_ODEs!, init_new, params_new, species_full, labels_new, "Rd")
rd_p1
# rd_p2

rh_range = range(0, 10, length=6)
include("fullmodel_params.jl")
rh_p1, rh_p2 = logplots(rh_range, init_new, 9, new_ODEs!, init_new, params_new, species_full, labels_new, "Rh")
rh_p1
# rh_p2

atp_range = range(0, 2, length=6) 
include("fullmodel_params.jl")
atp_p1, atp_p2 = logplots(atp_range, params_new, 5, new_ODEs!, init_new, params_new, species_full, labels_new, "ATP")
atp_p1
# atp_p2

kr_range = range(0.01, 2, length=6)
include("fullmodel_params.jl") 
kr_p1, kr_p2 = logplots(kr_range, params_new, 1, new_ODEs!, init_new, params_new, species_full, labels_new, "kr")
kr_p1
# kr_p2

L_range = range(0, 100, length=6) 
include("fullmodel_params.jl")
L_p1, L_p2 = logplots(L_range, params_new, 2, new_ODEs!, init_new, params_new, species_full, labels_new, "L")
L_p1
# L_p2

c_range = range(0, 100, length=6)
include("fullmodel_params.jl")
c_p1, c_p2 = logplots(c_range, params_new, 3, new_ODEs!, init_new, params_new, species_full, labels_new, "c")
c_p1
# c_p2

sigma_range = range(0, 0.2, length=6)
include("fullmodel_params.jl")
sigma_p1, sigma_p2 = logplots(sigma_range, params_new, 4, new_ODEs!, init_new, params_new, species_full, labels_new, "sigma")
sigma_p1
# sigma_p2

k1_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1_p1, k1_p2 = logplots(k1_range, params_new, 11, new_ODEs!, init_new, params_new, species_full, labels_new, "k1")
k1_p1
# k1_p2

k2_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2_p1, k2_p2 = logplots(k2_range, params_new, 12, new_ODEs!, init_new, params_new, species_full, labels_new, "k2")
k2_p1
# k2_p2

k3_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k3_p1, k3_p2 = logplots(k3_range, params_new, 13, new_ODEs!, init_new, params_new, species_full, labels_new, "k3")
k3_p1
# k3_p2

k4_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k4_p1, k4_p2 = logplots(k4_range, params_new, 14, new_ODEs!, init_new, params_new, species_full, labels_new, "k4")
k4_p1
# k4_p2

k5_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k5_p1, k5_p2 = logplots(k5_range, params_new, 15, new_ODEs!, init_new, params_new, species_full, labels_new, "k5")
k5_p1
# k5_p2

wrtc_range = range(0, 10, length=6)
include("fullmodel_params.jl")
wrtc_p1, wrtc_p2 = logplots(wrtc_range, params_new, 6, new_ODEs!, init_new, params_new, species_full, labels_new, "wrtc")
wrtc_p1
# wrtc_p2

thrtc_range = range(0, 20, length=6)
include("fullmodel_params.jl")
thrtc_p1, thrtc_p2 = logplots(thrtc_range, params_new, 7, new_ODEs!, init_new, params_new, species_full, labels_new, "thrtc")
thrtc_p1
# thrtc_p2

max_range = range(0, 10, length=6)
include("fullmodel_params.jl")
max_p1, max_p2 = logplots(max_range, params_new, 8, new_ODEs!, init_new, params_new, species_full, labels_new, "max")
max_p1
# max_p2

thr_range = range(0, 10, length=6)
include("fullmodel_params.jl")
thr_p1, thr_p2 = logplots(thr_range, params_new, 9, new_ODEs!, init_new, params_new, species_full, labels_new, "thr")
thr_p1
# thr_p2

# d_range = range(0, 5, length=6)
# include("fullmodel_params.jl")
# d_p1, d_p2 = logplots(d_range, params_new, 10, new_ODEs!, init_new, params_new, species_full, labels_new, "d")
# d_p1
# d_p2

k1a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1a_p1, k1a_p2 = logplots(k1a_range, params_new, 16, new_ODEs!, init_new, params_new, species_full, labels_new, "k1a")
k1a_p1
# k1a_p2

k2a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2a_p1, k2a_p2 = logplots(k2a_range, params_new, 17, new_ODEs!, init_new, params_new, species_full, labels_new, "k2a")
k2a_p1
# k2a_p2

k3a_range = range(0, 0.1, length=6)
include("fullmodel_params.jl")
k3a_p1, k3a_p2 = logplots(k3a_range, params_new, 18, new_ODEs!, init_new, params_new, species_full, labels_new, "k3a")
k3a_p1
# k3a_p2

kab_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kab_p1, kab_p2 = logplots(kab_range, params_new, 19, new_ODEs!, init_new, params_new, species_full, labels_new, "kab")
kab_p1
# kab_p2

kbb_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kbb_p1, kbb_p2 = logplots(kbb_range, params_new, 20, new_ODEs!, init_new, params_new, species_full, labels_new, "kbb")
kbb_p1
# kbb_p2

kcb_range = range(0,2, length=6)
include("fullmodel_params.jl")
kcb_p1, kcb_p2 = logplots(kcb_range, params_new, 21, new_ODEs!, init_new, params_new, species_full, labels_new, "kcb")
kcb_p1
# kcb_p2

k_range = range(0, 1, length=6)
include("fullmodel_params.jl")
k_p1, k_p2 = logplots(k_range, params_new, 22, new_ODEs!, init_new, params_new, species_full, labels_new, "k")
k_p1
k_p2

rtca_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtca_p1, rtca_p2 = logplots(rtca_range, init_new, 2, new_ODEs!, init_new, params_new, species_full, labels_new, "rtca")
# rtca_p1
rtca_p2

rtcb_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtcb_p1, rtcb_p2 = logplots(rtcb_range, init_new, 4, new_ODEs!, init_new, params_new, species_full, labels_new, "rtcb")
rtcb_p1
# rtcb_p2

rt_range = range(0, 10, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in full model to ode and see if that works (probs will) 
include("fullmodel_params.jl")
rt_p1, rt_p2 = logplots(rt_range, init_new, 11, new_ODEs!, init_new, params_new, species_full, labels_new, "Rt")
rt_p1
# rt_p2












atp_range_long = range(0,20, length=20)
include("fullmodel_params.jl")
atp_plot = full_plot_long_range_results(atp_range_long, params_new, 5, new_ODEs!, init_new, params_new, species_full, "ATP")

kr_range_long = range(0.01,100, length=20)
include("fullmodel_params.jl")
kr_plot = full_plot_long_range_results(kr_range_long, params_new, 1, new_ODEs!, init_new, params_new, species_full, "kr")

L_range_long = range(0,100, length=20)
include("fullmodel_params.jl")
L_plot = full_plot_long_range_results(L_range_long, params_new, 2, new_ODEs!, init_new, params_new, species_full, "L")

c_range_long = range(0,100, length=20)
include("fullmodel_params.jl")
c_plot = full_plot_long_range_results(c_range_long, params_new, 3, new_ODEs!, init_new, params_new, species_full, "c")

sigma_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
sigma_plot = full_plot_long_range_results(sigma_range_long, params_new, 4, new_ODEs!, init_new, params_new, species_full, "sigma")

k1_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k1_plot = full_plot_long_range_results(k1_range_long, params_new, 11, new_ODEs!, init_new, params_new, species_full, "k1")

k2_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k2_plot = full_plot_long_range_results(k2_range_long, params_new, 12, new_ODEs!, init_new, params_new, species_full, "k2")

k3_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k3_plot = full_plot_long_range_results(k3_range_long, params_new, 13, new_ODEs!, init_new, params_new, species_full, "k3")

k4_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k4_plot = full_plot_long_range_results(k4_range_long, params_new, 14, new_ODEs!, init_new, params_new, species_full, "k4")

k5_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k5_plot = full_plot_long_range_results(k5_range_long, params_new, 15, new_ODEs!, init_new, params_new, species_full, "k5")

wrtc_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
wrtc_plot = full_plot_long_range_results(wrtc_range_long, params_new, 6, new_ODEs!, init_new, params_new, species_full, "wrtc")

theta_range_long = range(0,20, length=20)
include("fullmodel_params.jl")
theta_plot = full_plot_long_range_results(theta_range_long, params_new, 7, new_ODEs!, init_new, params_new, species_full, "theta")

max_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
max_plot = full_plot_long_range_results(max_range_long, params_new, 8, new_ODEs!, init_new, params_new, species_full, "max")

thr_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
thr_plot = full_plot_long_range_results(thr_range_long, params_new, 9, new_ODEs!, init_new, params_new, species_full, "thr")

k1a_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
k1a_plot = full_plot_long_range_results(k1a_range_long, params_new, 16, new_ODEs!, init_new, params_new, species_full, "k1a")

k2a_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
k2a_plot = full_plot_long_range_results(k2a_range_long, params_new, 17, new_ODEs!, init_new, params_new, species_full, "k2a")

k3a_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
k3a_plot = full_plot_long_range_results(k3a_range_long, params_new, 18, new_ODEs!, init_new, params_new, species_full, "k3a")

kab_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
kab_plot = full_plot_long_range_results(kab_range_long, params_new, 19, new_ODEs!, init_new, params_new, species_full, "kab")

kbb_range_long = range(0,2, length=20)
include("fullmodel_params.jl")
kbb_plot = full_plot_long_range_results(kbb_range_long, params_new, 20, new_ODEs!, init_new, params_new, species_full, "kbb")

kcb_range_long = range(0,1, length=20)
include("fullmodel_params.jl")
kcb_plot = full_plot_long_range_results(kcb_range_long, params_new, 21, new_ODEs!, init_new, params_new, species_full, "kcb")

k_range_long = range(0,10, length=20)
include("fullmodel_params.jl")
k_plot = full_plot_long_range_results(k_range_long, params_new, 22, new_ODEs!, init_new, params_new, species_full, "k")


