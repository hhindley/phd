using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");


br = get_br(rtc_mod, params_bf, rtc_init, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1, rtcb2, rtcb3 = plot_rtc_bf(df, kdam1, kdam2, :rtcb, "1", "b693ccff", "RtcB")
rtcr1, rtcr2, rtcr3 = plot_rtc_bf(df, kdam1, kdam2, :rtcr, "2", "4ca7a2ff", "RtcR")
rtca1, rtca2, rtca3 = plot_rtc_bf(df, kdam1, kdam2, :rtca, "3", "e48080ff", "RtcA")

p = plot([rtcb1, rtcb2, rtcb3, rtcr1, rtcr2, rtcr3, rtca1, rtca2, rtca3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rh1, rh2, rh3 = plot_rtc_bf(df, kdam1, kdam2, :rh, "1", "ffd30cff", "RtcB")
rt1, rt2, rt3 = plot_rtc_bf(df, kdam1, kdam2, :rt, "2", "e96100ff", "rt")
rd1, rd2, rd3 = plot_rtc_bf(df, kdam1, kdam2, :rd, "3", "ac0606ff", "rd")

p1 = plot([rh1, rh2, rh3, rt1, rt2, rt3, rd1, rd2, rd3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

[p1_sig p1]

savefig(p, "/home/holliehindley/phd/may23_rtc/paper_plots/rtc_proteins.svg")
savefig(p1, "/home/holliehindley/phd/may23_rtc/paper_plots/ribosomes.svg")

kdam_range1 = range(0,1.5,length=100)
kdam_range2 = range(1.5,0,length=100)
res = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rh, species_rtc, kdam_range1)
res1 = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rh, species_rtc, kdam_range2)
plot([scatter(x=kdam_range1,y=res), scatter(x=kdam_range2,y=res1)])

res

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
