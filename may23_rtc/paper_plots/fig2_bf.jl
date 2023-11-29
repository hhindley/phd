using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");

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

tspan = (0,1e9)

params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

params_ = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)


br = get_br(rtc_mod, params2, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1, rtcb2, rtcb3, bf_rtcb = plot_rtc_bf(bf, df, kdam1, kdam2, :rtcb)
rtcr1, rtcr2, rtcr3, bf_rtcr = plot_rtc_bf(bf, df, kdam1, kdam2, :rtcr)
rtca1, rtca2, rtca3, bf_rtca = plot_rtc_bf(bf, df, kdam1, kdam2, :rtca)

p = plot([rtcb1, rtcb2, rtcb3, bf_rtcb, rtcr1, rtcr2, rtcr3, bf_rtcr, rtca1, rtca2, rtca3, bf_rtca],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))

rh1, rh2, rh3, bf_rh = plot_rtc_bf(bf, df, kdam1, kdam2, :rh)
rt1, rt2, rt3, bf_rt = plot_rtc_bf(bf, df, kdam1, kdam2, :rt)
rd1, rd2, rd3, bf_rd = plot_rtc_bf(bf, df, kdam1, kdam2, :rd)

p1 = plot([rh1, rh2, rh3, bf_rh, rt1, rt2, rt3, bf_rt, rd1, rd2, rd3, bf_rd],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


savefig(p, "/home/holliehindley/phd/may23_rtc/paper_plots/rtc_proteins.svg")
savefig(p1, "/home/holliehindley/phd/may23_rtc/paper_plots/ribosomes.svg")


p1 = plot([rh1, rh2, rh3, bf_rh],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
p2 = plot([rt1, rt2, rt3, bf_rt],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rt (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
p3 = plot([rtcb1, rtcb2, rtcb3, bf_rtcb],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
p4 = plot([rtcr1, rtcr2, rtcr3, bf_rtcr],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcR (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


savefig(p1, "/home/holliehindley/phd/may23_rtc/paper_plots/rh_bf.svg")
savefig(p2, "/home/holliehindley/phd/may23_rtc/paper_plots/rt_bf.svg")
savefig(p3, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcb_bf.svg")
savefig(p4, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcr_bf.svg")
