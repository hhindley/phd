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

params_for_ssval_setup = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)
params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

tspan=(0,1e9)
# ps = deepcopy(params2)
ps.kdam = 0.3695
branches = setup_ssvals_from_bfkit(0.3695, params_for_ssval_setup2)
n=1; l=1000
all_ranges = get_all_ranges(set_ss_range_zerotossval, branches, "ss_val_on", n, l)
all, init_vals = get_rh_init_switch_all_ranges(all_ranges, branches.ss_val_on, :rh, l, ps)

full_find_differences_or_percs(all,get_percentages,init_vals,branches.ss_val_off[7],l,branches.ss_val_on,l,"on")


params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
kdam_range = range(0.63575,2.0175,length=100)
# logkdam_range = 10 .^ range(log10(2.016996),log10(0.6359851),length=500)
# kdam_range1 = vcat(range(0.6359851,1.479,length=10), range(1.48,2.016996,length=200))

tspan=(0,1e9)

all_res = []
init_res = []
switches = DataFrame(rtcbs=[],rhs=[],rtcrs=[],rts=[])
psm = deepcopy(params1)

for kdam_val in ProgressBar(kdam_range)
    psm = deepcopy(params1)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup)
    # @show psm
    
    n = 600; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[9]
    all, init_vals = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_on,:rh,l,psm)
    push!(all_res, all)
    push!(init_res, init_vals)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l)
    switch_ind = get_switch_ind(binary, l)
    switch_vals = get_switch_vals(switch_ind, init_vals)
    push!(switches.rtcbs, switch_vals[4])
    push!(switches.rhs, switch_vals[7])
    push!(switches.rtcrs, switch_vals[6])
    push!(switches.rts, switch_vals[9])
end

all_res
init_res
switches

rtcbs=[]
for i in switches
    push!(rtcbs, i[4])
end
rtcbs

switch_rtcb = scatter(x=kdam_range,y=switches.rtcbs, fill="tonexty")
switch_rtcr = scatter(x=kdam_range,y=switches.rtcrs)
switch_rh = scatter(x=kdam_range,y=switches.rhs, yaxis="y2")
switch_rt = scatter(x=kdam_range,y=switches.rts, yaxis="y2")

br = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)


colors=[:purple,:mediumpurple,:green,:plum,:purple,:green]
colors_r=[:gold,:darkorange,:red]
first,middle,last=dashed_lines_species(df, bf, colors)
r_first,r_middle,r_last = dashed_lines_ribosomes(df,bf,colors_r)
bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr, bf_rh, bf_rd, bf_rt = bf_scatter(bf, "darkblue")

rtcb_plot = plot([first[4],middle[4],last[4],bf_rtcb,
switch_rtcb],)
# Layout(legend=attr(x=0.75,y=1),yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)",
# yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))









rtcr_plot = plot([first[6],middle[6],last[6],bf_rtcr,
switch_rtcr],
Layout(legend=attr(x=0.75,y=1),yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))

rh_plot = plot([r_first[1],r_middle[1],r_last[1],bf_rh,
switch_rh],
Layout(legend=attr(x=0.75,y=1),yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))

rt_plot = plot([r_first[3],r_middle[3],r_last[3],bf_rt,
switch_rt],
Layout(legend=attr(x=0.75,y=1),yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))










