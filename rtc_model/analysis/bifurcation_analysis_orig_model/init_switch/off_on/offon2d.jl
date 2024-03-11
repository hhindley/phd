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

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# kdam_range = range(0.6359851,2.016996,length=500)
# logkdam_range = 10 .^ range(log10(0.6359851),log10(2.016996),length=500)
# kdam_range1 = vcat(range(0.6359851,0.643,length=200), range(0.6451,2.016996,length=10))

tspan=(0,1e9)


# lower to upper - double change 
# n = 0.5; l = 100;
# branches1 = setup_ssvals_from_bfkit(0.6359851)
# ps = deepcopy(params1)
# ps.kdam = 0.6359851
# lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

# sss, rtcbs, rtcrs = double_init(branches1.ss_val_off, lower_ranges, branches1.ss_val_on, "rm_b", "rh", ps)
# # findall(x->x==branches1.ss_val_upper[7],sss)

# p = plot_binary(sss, set_shared_range_0ton(n,l), branches1.ss_val_off[7], "rtcr", "rt", "lower to upper - kdam = 1")

# # # open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper1.html", "w") do io
# # #     PlotlyBase.to_html(io, p.plot)
# # # end
# binary = get_binary(branches1.ss_val_off[7], sss)
# x, y = get_plot_vals(binary, l, n)
# plot(scatter(x=x,y=y), Layout(xaxis_range=(0,1000), yaxis_range=(0,1000)))




kdam_range = range(0.6359851,2.016996,length=10)
tspan=(0,1e9)
sss_res=[]; branches=[]; xs=[]; ys=[]

for kdam_val in ProgressBar(kdam_range)
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val)

    
    n = 25; l = 500;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

    sss, rtcbs, rtcrs = double_init(branches1.ss_val_off, lower_ranges, branches1.ss_val_on, "rtcr", "rt", ps)
    push!(sss_res, sss)
    push!(branches,branches1)
    binary = get_binary(branches1.ss_val_off[7], sss)
    x, y = get_plot_vals(binary, l, n)
    push!(xs, x)
    push!(ys, y)
end


data = DataFrame(sss=sss_res,xs=xs,ys=ys,branches=branches)
CSV.write("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/data_offon_2d.csv", data)


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

p6 = (scatter(x=xs[1],y=ys[1], name="$(kdam_range[1])"))
p7 = (scatter(x=xs[2],y=ys[2], name="$(kdam_range[2])"))
p8 = (scatter(x=xs[3],y=ys[3], name="$(kdam_range[3])"))
p9 = (scatter(x=xs[4],y=ys[4], name="$(kdam_range[4])"))
p10 = (scatter(x=xs[5],y=ys[5], name="$(kdam_range[5])"))
p11 = (scatter(x=xs[6],y=ys[6], name="$(kdam_range[6])"))
p12 = (scatter(x=xs[7],y=ys[7], name="$(kdam_range[7])"))
p13 = (scatter(x=xs[8],y=ys[8], name="$(kdam_range[8])"))
p14 = (scatter(x=xs[9],y=ys[9], name="$(kdam_range[9])"))
p15 = (scatter(x=xs[10],y=ys[10], name="$(kdam_range[10])"))

p1 = plot([p6,p7,p8,p9,p10,p11,p12,p13,p14,p15], Layout(xaxis_range=(0,2500), yaxis_range=(0,2500), xaxis_title="RtcR", yaxis_title="Rt", title="line indicates turning point of off→on switch (% increase)"))

open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/offon_double_rt_rtcr.html", "w") do io
    PlotlyBase.to_html(io, p1.plot)
end

# [p2 p3 p4; p5 p6 p7; p8 p9 p10]



# # lower to upper - triple change
# n = 5; l = 20;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches1, "ss_val_lower", n, l)

# sss, init1s, init2s, init3s = triple_init(branches1.ss_val_lower, lower_ranges, branches1.ss_val_upper, "rtcr", "rt", "rtcb")
# ind = findall(x->x==upper_branch.ss_val[7],sss)
# sss
# pars=[]
# for i in range(1,length(ind))
#     push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
# end
# pars










# params2 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.64, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

# branches2 = setup_ssvals(params2)

# # lower to upper - single change 
# n = 35; l = 1000;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)
# all = get_rh_init_switch_all_ranges(lower_ranges, branches2.ss_val_lower, :rh, l, params2)
# dfs = upper_or_lower(all, branches2.ss_val_lower, l)
# shared_range = set_shared_range(n,l)

# p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
# scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
# scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
# Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
# yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from lower to upper (rh) - kdam = 0.64"
# ))


# open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowerupper_064.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end

# # lower to upper - double change 
# n = 1; l = 100;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

# sss, rtcbs, rtcrs = double_init(branches2.ss_val_lower, lower_ranges, branches2.ss_val_upper, "rtcr", "rt", params2)
# findall(x->x==branches2.ss_val_upper[7],sss)

# p = plot_binary(sss, set_shared_range(n,l), branches2.ss_val_lower, "rt", "rtcr", "lower to upper - kdam = 0.64")

# open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end


# # lower to upper - triple change
# n = 5; l = 20;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

# sss, init1s, init2s, init3s = triple_init(lower_branch, lower_ranges, upper_branch, "rtcr", "rt", "rtcb")
# ind = findall(x->x==upper_branch.ss_val[7],sss)
# sss
# pars=[]
# for i in range(1,length(ind))
#     push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
# end
# pars









# params3 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 2, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

# branches3 = setup_ssvals(params3)

# # lower to upper - single change 
# n = 1000; l = 1000;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_lower", n, l)
# all = get_rh_init_switch_all_ranges(lower_ranges, branches3.ss_val_lower, :rh, l, params3)
# dfs = upper_or_lower(all, "lower", l)
# shared_range = set_shared_range(n,l)

# p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
# scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
# scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
# Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
# yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from lower to upper (rh) - kdam = 2"
# ))


# open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowerupper2.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end


# # lower to upper - double change 
# n = 1; l = 100;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

# sss, rtcbs, rtcrs = double_init(branches2.ss_val_lower, lower_ranges, branches2.ss_val_upper, "rtcr", "rt", params2)
# findall(x->x==branches2.ss_val_upper[7],sss)

# p = plot_binary(sss, set_shared_range(n,l), branches2.ss_val_lower, "rt", "rtcr", "lower to upper - kdam = 0.64")

# open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end


# # lower to upper - triple change
# n = 5; l = 20;
# lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

# sss, init1s, init2s, init3s = triple_init(lower_branch, lower_ranges, upper_branch, "rtcr", "rt", "rtcb")
# ind = findall(x->x==upper_branch.ss_val[7],sss)
# sss
# pars=[]
# for i in range(1,length(ind))
#     push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
# end
# pars

