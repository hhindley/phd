using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, LabelledArrays, DataFrames
# using Plots
using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");


wr_range = 10 .^ range(log10(1e-7),log10(2e-6),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wr_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_r=>i,))
    br = get_br(rtc_mod, params, rtc_init, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


wr_rtcb1, wr_rtcb2, wr_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))")
wr_rtcb1a, wr_rtcb2a, wr_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>r</sub> = $(round.(wr_range[2]; sigdigits=2))")
wr_rtcb1b, wr_rtcb2b, wr_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))")

wr_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))", legendgroup="3")

wr = plot([wr_rtcb1, wr_rtcb2, wr_rtcb3, wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_rtcb_nonbs],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

PlotlyJS.savefig(wr, "/home/holliehindley/phd/may23_rtc/paper_plots/wr_rtcb.svg")



wab_range = 10 .^range(log10(9e-5),log10(7e-4),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wab_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_ab=>i,))
    br = get_br(rtc_mod, params, initial, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wab_rtcb1, wab_rtcb2, wab_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))")
wab_rtcb1a, wab_rtcb2a, wab_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>ab</sub> = $(round.(wab_range[2]; sigdigits=2))")
wab_rtcb1b, wab_rtcb2b, wab_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))")

wab_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))", legendgroup="3")

wab = plot([wab_rtcb1, wab_rtcb2, wab_rtcb3,wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_rtcb1b, wab_rtcb2b, wab_rtcb3b],#wab_rtcb_nonbs],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


PlotlyJS.savefig(wab, "/home/holliehindley/phd/may23_rtc/paper_plots/wab_rtcb.svg")


atp_range = range(1000,stop=5000,length=3)

bfs=[]; dfs=[];
for i in ProgressBar(atp_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:atp=>i,))
    br = get_br(rtc_mod, params, rtc_init, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


atp_rtcb1, atp_rtcb2, atp_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ATP = $(round.(atp_range[1]; sigdigits=2))")
atp_rtcb1a, atp_rtcb2a, atp_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ATP = $(round.(atp_range[2]; sigdigits=2))")
atp_rtcb1b, atp_rtcb2b, atp_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ATP = $(round.(atp_range[3]; sigdigits=2))")

atp_p = plot([atp_rtcb1, atp_rtcb2, atp_rtcb3, atp_rtcb1a, atp_rtcb2a, atp_rtcb3a, atp_rtcb1b, atp_rtcb2b, atp_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


PlotlyJS.savefig(atp_p, "/home/holliehindley/phd/may23_rtc/paper_plots/atp_rtcb.svg")


lam_range = range(0.012,stop=0.02,length=3)

bfs=[]; dfs=[];
for i in ProgressBar(lam_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod, params, rtc_init, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


lam_rtcb1, lam_rtcb2, lam_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "λ = $(round.(lam_range[1]; sigdigits=2))")
lam_rtcb1a, lam_rtcb2a, lam_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "λ = $(round.(lam_range[2]; sigdigits=2))")
lam_rtcb1b, lam_rtcb2b, lam_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "λ = $(round.(lam_range[3]; sigdigits=2))")

lam_rtcb_nonbs = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="λ = $(round.(lam_range[1]; sigdigits=2))", legendgroup="1")

lam_p = plot([lam_rtcb_nonbs, lam_rtcb1a, lam_rtcb2a, lam_rtcb3a, lam_rtcb1b, lam_rtcb2b, lam_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

PlotlyJS.savefig(lam_p, "/home/holliehindley/phd/may23_rtc/paper_plots/lam_rtcb.svg")


