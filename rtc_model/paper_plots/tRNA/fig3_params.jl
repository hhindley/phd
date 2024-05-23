using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, PlotlyJS, ProgressBars, ModelingToolkit, Interpolations, QuadGK
# using Plots

PATH = "/home/holliehindley/phd"

include("$PATH/general_funcs/solving.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/rtc_model/models/rtc_trna_model.jl")
include("$PATH/rtc_model/functions/bf_funcs/bf_funcs.jl")

margins = attr(l=180,r=100,t=100,b=100)
axes_linewidth = 5.5
axes_title_fs = 34
tick_fs = attr(size=38)


kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)


lam_range = [0.01,0.013,0.015,0.02]#range(0.01,stop=0.02,length=4)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(lam_range)
    params_trna1[lam] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


# lam_rtcb1, lam_rtcb2, lam_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "λ = $(round.(lam_range[1]; sigdigits=2))","1")
lam_rtcb1a, lam_rtcb2a, lam_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "λ = $(round.(lam_range[2]; sigdigits=2))","1")
lam_rtcb1b, lam_rtcb2b, lam_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "λ = $(round.(lam_range[3]; sigdigits=2))","1")
lam_rtcb1c, lam_rtcb2c, lam_rtcb3c = plot_rtc_bf(dfs[4], findall(x->x==bfs[4].kdam[1],dfs[4].kdam)[1], findall(x->x==bfs[4].kdam[2],dfs[4].kdam)[1], :rtcb, "3", "800080ff", "λ = $(round.(lam_range[4]; sigdigits=2))","1")

lam_rtcb_nonbs = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="λ = $(round.(lam_range[4]; sigdigits=2))", legendgroup="3")
# lam_rtcb_nonbs1 = scatter(x=dfs[2].kdam,y=dfs[2].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="λ = $(round.(lam_range[4]; sigdigits=2))", legendgroup="3")

lam_p = plot([lam_rtcb_nonbs,lam_rtcb1a, lam_rtcb2a, lam_rtcb3a, lam_rtcb1b, lam_rtcb2b, lam_rtcb3b, lam_rtcb1c, lam_rtcb2c, lam_rtcb3c],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,4.2],tickvals=[0,1,2,3,4],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))

savefig(lam_p, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/lam_rtcb.svg")


atp_range = [1000,2500,3500,5000]#range(1000,stop=5000,length=4)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(atp_range)
    params_trna1[atp] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

atp_rtcb1, atp_rtcb2, atp_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ATP = $(round.(atp_range[1]; sigdigits=2))", "1")
atp_rtcb1a, atp_rtcb2a, atp_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ATP = $(round.(atp_range[2]; sigdigits=2))", "1")
atp_rtcb1b, atp_rtcb2b, atp_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ATP = $(round.(atp_range[3]; sigdigits=2))", "1")
atp_rtcb1c, atp_rtcb2c, atp_rtcb3c = plot_rtc_bf(dfs[4], findall(x->x==bfs[4].kdam[1],dfs[4].kdam)[1], findall(x->x==bfs[4].kdam[2],dfs[4].kdam)[1], :rtcb, "4", "800080ff", "ATP = $(round.(atp_range[4]; sigdigits=2))", "1")

atp_p = plot([atp_rtcb1, atp_rtcb2, atp_rtcb3, atp_rtcb1a, atp_rtcb2a, atp_rtcb3a, atp_rtcb1b, atp_rtcb2b, atp_rtcb3b, atp_rtcb1c, atp_rtcb2c, atp_rtcb3c],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,2.3],tickvals=[0,1,2,3,4,5],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))

savefig(atp_p, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/atp_rtcb.svg")


wab_range =  [2e-6,5e-6,1e-5,2e-5]#10 .^range(log10(2e-6),log10(2e-5),length=4)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(wab_range)
    params_trna1[ω_ab] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end
wab_rtcb1, wab_rtcb2, wab_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))", "1")
wab_rtcb1a, wab_rtcb2a, wab_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>ab</sub> = $(round.(wab_range[2]; sigdigits=2))", "1")
wab_rtcb1b, wab_rtcb2b, wab_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))", "1")
wab_rtcb1c, wab_rtcb2c, wab_rtcb3c = plot_rtc_bf(dfs[4], findall(x->x==bfs[4].kdam[1],dfs[4].kdam)[1], findall(x->x==bfs[4].kdam[2],dfs[4].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>ab</sub> = $(round.(wab_range[4]; sigdigits=2))", "1")

wab_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))", legendgroup="3")
wab_rtcb_nonbs = scatter(x=dfs[4].kdam,y=dfs[4].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[4]; sigdigits=2))", legendgroup="3")

wab = plot([wab_rtcb_nonbs, wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_rtcb1b, wab_rtcb2b, wab_rtcb3b, wab_rtcb_nonbs1], #wab_rtcb1, wab_rtcb2, wab_rtcb3
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,4.2],tickvals=[0,1,2,3,4,5],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))

savefig(wab, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/wab_rtcb.svg")


wr_range = [2e-7,5e-7,1e-6,2e-6]#10 .^ range(log10(2e-7),log10(2e-6),length=4)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(wr_range)
    params_trna1[ω_r] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wr_rtcb1, wr_rtcb2, wr_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))", "1")
wr_rtcb1a, wr_rtcb2a, wr_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>r</sub> = $(round.(wr_range[2]; sigdigits=2))", "1")
wr_rtcb1b, wr_rtcb2b, wr_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))", "1")

wr_rtcb_nonbs = scatter(x=dfs[4].kdam,y=dfs[4].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[4]; sigdigits=2))", legendgroup="3")
wr_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))", legendgroup="3")

wr = plot([wr_rtcb_nonbs1, wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_rtcb1b, wr_rtcb2b, wr_rtcb3b, wr_rtcb_nonbs],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,4],tickvals=[0,1,2,3,4,5],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))


savefig(wr, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/wr_rtcb.svg")


rh_range = range(20,stop=50,length=4)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(rh_range)
    params_trna1[rh] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

rh_rtcb1, rh_rtcb2, rh_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "Rh = $(round.(rh_range[4]; sigdigits=2))", "1")
rh_rtcb1a, rh_rtcb2a, rh_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "Rh = $(round.(rh_range[2]; sigdigits=2))", "1")
rh_rtcb1b, rh_rtcb2b, rh_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "Rh = $(round.(rh_range[3]; sigdigits=2))", "1")

rh_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="Rh = $(round.(rh_range[1]; sigdigits=2))", legendgroup="3")
rh_rtcb_nonbs2 = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="Rh = $(round.(rh_range[3]; sigdigits=2))", legendgroup="3")
rh_rtcb_nonbs3 = scatter(x=dfs[4].kdam,y=dfs[4].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="Rh = $(round.(rh_range[1]; sigdigits=2))", legendgroup="3")

rh_p = plot([rh_rtcb1, rh_rtcb2, rh_rtcb3, rh_rtcb1a, rh_rtcb2a, rh_rtcb3a, rh_rtcb1b, rh_rtcb2b, rh_rtcb3b, rh_rtcb_nonbs3],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,5.5],tickvals=[0,1,2,3,4,5],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))

savefig(rh_p, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/rh_rtcb.svg")



thr_range = range(0.28,stop=4,length=3)

bfs=[]; dfs=[];# res1=[]; res2=[];
params_trna1 = deepcopy(params_trna)
for i in ProgressBar(thr_range)
    params_trna1[thr_t] = i
    br = get_br(rtc_trna_model, ssvals_trna, params_trna1, 240.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

thr_rtcb1, thr_rtcb2, thr_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "θ_tRNA = $(round.(thr_range[1]; sigdigits=2))")
thr_rtcb1a, thr_rtcb2a, thr_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "θ_tRNA = $(round.(thr_range[2]; sigdigits=2))")
thr_rtcb1b, thr_rtcb2b, thr_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "θ_tRNA = $(round.(thr_range[3]; sigdigits=2))")

# rh_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="Rh = $(round.(rh_range[3]; sigdigits=2))", legendgroup="3")
# thr_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="θ_tRNA = $(round.(thr_range[1]; sigdigits=2))", legendgroup="3")
thr_p = plot([thr_rtcb1, thr_rtcb2, thr_rtcb3, thr_rtcb1a, thr_rtcb2a, thr_rtcb3a, thr_rtcb1b, thr_rtcb2b, thr_rtcb3b],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

savefig(thr_p, "$PATHrtc_model/paper_plots/tRNA/thr_rtcb.svg")




kin_range = [0.0001,0.0005,0.005,0.01]#0.01]#range(0.0001,stop=0.0003,length=4)
# lam_range = range(0.01,stop=0.03,length=4)

bfs=[]; dfs=[];
copyparams = deepcopy(params_trna)
for i in ProgressBar(kin_range)
    copyparams[kin] = i
    br = get_br(rtc_trna_model, ssvals_trna, copyparams, 30.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

kin_rtcb1, kin_rtcb2, kin_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "λ = $(round.(kin_range[1]; sigdigits=2))", "1")
kin_rtcb1a, kin_rtcb2a, kin_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "λ = $(round.(kin_range[2]; sigdigits=2))", "1")
kin_rtcb1b, kin_rtcb2b, kin_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "λ = $(round.(kin_range[3]; sigdigits=2))", "1")
kin_rtcb1c, kin_rtcb2c, kin_rtcb3c = plot_rtc_bf(dfs[4], findall(x->x==bfs[4].kdam[1],dfs[4].kdam)[1], findall(x->x==bfs[4].kdam[2],dfs[4].kdam)[1], :rtcb, "3", "800080ff", "λ = $(round.(kin_range[4]; sigdigits=2))", "1")

kin_rtcb_nonbs = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="λ = $(round.(kin_range[1]; sigdigits=2))", legendgroup="1")
kin_rtcb_nonbs1 = scatter(x=dfs[4].kdam,y=dfs[4].rtcb,showlegend=true,line=attr(width=9, color="#c1c1c1ff"), name="λ = $(round.(kin_range[1]; sigdigits=2))", legendgroup="1")

kin_p = plot([kin_rtcb1, kin_rtcb2, kin_rtcb3, kin_rtcb1a, kin_rtcb2a, kin_rtcb3a, kin_rtcb1b, kin_rtcb2b, kin_rtcb3b, kin_rtcb_nonbs1],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,2.3],tickvals=[0,0.5,1,1.5,2],tickfont=tick_fs,automargin=true),xaxis=attr(showline=true,linewidth=axes_linewidth,linecolor="black",range=[0,30],tickvals=[0,10,20,30],tickfont=tick_fs,automargin=true),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=axes_title_fs, color="black", family="sans-serif"), showlegend=false, margin=margins))

PlotlyJS.savefig(kin_p, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/plots/kin_rtcb.svg")


 
# creating the banana plots 
colours = ["#f4e5ffff","#e6c5ffff","#d7a1ffff","#c06affff","#a730ffff"]

atp_range = range(500,stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)

wab_range1 = 10 .^range(log10(1e-7),log10(1e-3),length=5)
results_wab, xvals = get_bs_region_results(rtc_trna_model, ssvals_trna, atp_range, atp, lam_range, lam, wab_range1, ω_ab, params_trna, 50.)

auc1 = area_under_curve(xvals[1],results_wab[1][1])
auc2 = area_under_curve(xvals[2],results_wab[2][1])
auc3 = area_under_curve(xvals[3],results_wab[3][1])
auc4 = area_under_curve(xvals[4],results_wab[4][1])
auc5 = area_under_curve(xvals[5],results_wab[5][1])

wab_banana = plot([
scatter(x=xvals[1],y=results_wab[1][2], line=attr(color=colours[1])), scatter(x=auc1.x,y=auc1.f1.(auc1.x), fill="tonexty", mode="lines", line=attr(color=colours[1])),
scatter(x=xvals[2],y=results_wab[2][2], line=attr(color=colours[2])), scatter(x=auc2.x,y=auc2.f1.(auc2.x), fill="tonexty", mode="lines", line=attr(color=colours[2])),
scatter(x=xvals[3],y=results_wab[3][2], line=attr(color=colours[3])), scatter(x=auc3.x,y=auc3.f1.(auc3.x), fill="tonexty", mode="lines", line=attr(color=colours[3])),
scatter(x=xvals[4],y=results_wab[4][2], line=attr(color=colours[4])), scatter(x=auc4.x,y=auc4.f1.(auc4.x), fill="tonexty", mode="lines", line=attr(color=colours[4])),
scatter(x=xvals[5],y=results_wab[5][2], line=attr(color=colours[5])), scatter(x=auc5.x,y=auc5.f1.(auc5.x), fill="tonexty", mode="lines", line=attr(color=colours[5]))],
Layout(xaxis_title="ATP (μM)", yaxis_title="λ (min<sup>-1</sup>)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black",range=(0.001,0.04)),xaxis=attr(showline=true,linewidth=3,linecolor="black", range=(500,5500)),# showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


savefig(wab_banana, "$PATHrtc_model/paper_plots/tRNA/wab_banana_final.svg")


wr_range1 = 10 .^ range(log10(1e-8),log10(1e-4),length=5)
results_wr, xvals_wr = get_bs_region_results(rtc_trna_model, ssvals_trna, atp_range, atp, lam_range, lam, wr_range1, ω_r, params_trna, 50.)


auc1_wr = area_under_curve(xvals_wr[1],results_wr[1][1])
auc2_wr = area_under_curve(xvals_wr[2],results_wr[2][1])
auc3_wr = area_under_curve(xvals_wr[3],results_wr[3][1])
auc4_wr = area_under_curve(xvals_wr[4],results_wr[4][1])
auc5_wr = area_under_curve(xvals_wr[5],results_wr[5][1])

wr_banana = plot([
scatter(x=xvals_wr[1],y=results_wr[1][2], line=attr(color=colours[1])), scatter(x=auc1_wr.x,y=auc1_wr.f1.(auc1_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[1])),
scatter(x=xvals_wr[2],y=results_wr[2][2], line=attr(color=colours[2])), scatter(x=auc2_wr.x,y=auc2_wr.f1.(auc2_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[2])),
scatter(x=xvals_wr[3],y=results_wr[3][2], line=attr(color=colours[3])), scatter(x=auc3_wr.x,y=auc3_wr.f1.(auc3_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[3])),
scatter(x=xvals_wr[4],y=results_wr[4][2], line=attr(color=colours[4])), scatter(x=auc4_wr.x,y=auc4_wr.f1.(auc4_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[4])),
scatter(x=xvals_wr[5],y=results_wr[5][2], line=attr(color=colours[5])), scatter(x=auc5_wr.x,y=auc5_wr.f1.(auc5_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[5]))],
Layout(xaxis_title="ATP (μM)", yaxis_title="λ (min<sup>-1</sup>)", showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black",range=(0.001,0.04)),xaxis=attr(showline=true,linewidth=3,linecolor="black", range=(500,5500)),# showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

savefig(wr_banana, "$PATHrtc_model/paper_plots/tRNA/wr_banana_final.svg")











