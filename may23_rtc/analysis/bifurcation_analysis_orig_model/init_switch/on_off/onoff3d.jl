using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, ProgressBars
include("/home/hollie_hindley/Documents/may23_rtc/functions/solving.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/set_ups.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/plotting.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/functions/sweep_params.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/rtc_orig.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/models/single_t.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/combinations_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

#branches = setup_ssvals_from_bfkit(1.)
#n=6000; l=10
#all_ranges = get_all_ranges(set_ss_range_zerotossval, branches, "ss_val_on", n, l)
#all, init_vals, unstable = get_rh_init_switch_all_ranges(all_ranges, branches.ss_val_on, :rh, l, params1)




params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# kdam_range = range(0.6359851,2.016996,length=500)
# logkdam_range = 10 .^ range(log10(2.016996),log10(0.6359851),length=500)
# kdam_range1 = vcat(range(0.6359851,1.479,length=10), range(1.48,2.016996,length=200))

tspan=(0,1e9)

# psm = deepcopy(params1)

# # upper to lower - double change - where we see change rtcr and rtcb
n1 = 1; l1= 50;
branches1 = setup_ssvals_from_bfkit(2.016996)
ps = deepcopy(params1)
ps.kdam = 2.016996
upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n1,l1)

sss, init1s, init2s, init3s = triple_init(branches1.ss_val_on, upper_ranges, branches1.ss_val_off, "rtcr", "rtcb", "rh", ps)
sss
binary=get_binary(branches1.ss_val_off[7], sss)


kdam_range = range(0.63575,2.0175,length=5)

tspan=(0,1e9)
sss_res=[]; branches=[]; xs=[]; ys=[]
rtcrs=[]; rtcbs=[]; rhs=[]
for kdam_val in ProgressBar(kdam_range)
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup)

    
    n = 1; l = 100;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)

    sss, init1s, init2s, init3s = triple_init(branches1.ss_val_on, upper_ranges, branches1.ss_val_off, "rtcr", "rtcb", "rh", ps)
    push!(sss_res, sss)
    push!(rtcrs, init1s)
    push!(rtcbs, init2s)
    push!(rhs, init3s)
    push!(branches, branches1.ss_val_on[7])
end
sss_res[1]
branches[1]
rtcrs[1]
rtcbs[1]
rhs[1]
# p1 = plot_binary(sss_res[1], set_shared_range_0ton(n,l), branches[1].ss_val_off[7], "rtcr", "rt", "kdam = $(kdam_range[1])")
# p2 = plot_binary(sss_res[2], set_shared_range_0ton(n,l), branches[2].ss_val_off[7], "rtcr", "rt", "kdam = $(kdam_range[2])")
# p3 = plot_binary(sss_res[3], set_shared_range_0ton(n,l), branches[3].ss_val_off[7], "rtcr", "rt", "kdam = $(kdam_range[3])")
# p4 = plot_binary(sss_res[4], set_shared_range_0ton(n,l), branches[4].ss_val_off[7], "rtcr", "rt", "kdam = $(kdam_range[4])")
# p5 = plot_binary(sss_res[5], set_shared_range_0ton(n,l), branches[5].ss_val_off[7], "rtcr", "rt", "kdam = $(kdam_range[5])")
# # p6 = plot_binary(sss_res[6], set_shared_range_0ton(n,l), branches[6].ss_val_off, "rtcr", "rt", "kdam = $(kdam_range[6])")
# # p7 = plot_binary(sss_res[7], set_shared_range_0ton(n,l), branches[7].ss_val_off, "rtcr", "rt", "kdam = $(kdam_range[7])")
# # p8 = plot_binary(sss_res[8], set_shared_range_0ton(n,l), branches[8].ss_val_off, "rtcr", "rt", "kdam = $(kdam_range[8])")
# # p9 = plot_binary(sss_res[9], set_shared_range_0ton(n,l), branches[9].ss_val_off, "rtcr", "rt", "kdam = $(kdam_range[9])")
# # p10 = plot_binary(sss_res[10], set_shared_range_0ton(n,l), branches[10].ss_val_off, "rtcr", "rt", "kdam = $(kdam_range[10])")

# p = [p1 p2 p3; p4 p5]

# open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/range_kdam_double.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end

xs1 = []
ys1 = []
for (i,j) in zip(xs,ys)
    push!(xs1, 100 .- i)
    push!(ys1, 100 .- j)
end

p6 = (scatter(x=xs1[1],y=ys1[1], name="$(kdam_range[1])"))
p7 = (scatter(x=xs1[2],y=ys1[2], name="$(kdam_range[2])"))
p8 = (scatter(x=xs1[3],y=ys1[3], name="$(kdam_range[3])"))
p9 = (scatter(x=xs1[4],y=ys1[4], name="$(kdam_range[4])"))
p10 = (scatter(x=xs1[5],y=ys1[5], name="$(kdam_range[5])"))
p11 = (scatter(x=xs1[6],y=ys1[6], name="$(kdam_range[6])"))
p12 = (scatter(x=xs1[7],y=ys1[7], name="$(kdam_range[7])"))
p13 = (scatter(x=xs1[8],y=ys1[8], name="$(kdam_range[8])"))
p14 = (scatter(x=xs1[9],y=ys1[9], name="$(kdam_range[9])"))
p15 = (scatter(x=xs1[10],y=ys1[10], name="$(kdam_range[10])"))


p = plot([p6,p7,p8,p9,p10,p11,p12,p13,p14,p15], Layout(xaxis_title="RtcB", yaxis_title="RtcR", title="line indicates turning point of on→off switch (% decrease)"))


open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/range_kdam_double.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end










params3 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 2, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
