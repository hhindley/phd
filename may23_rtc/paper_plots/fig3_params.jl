using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, LabelledArrays, DataFrames
# using Plots
using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl")

wr_range = 10 .^ range(log10(2e-7),log10(2e-6),length=3)
wr_range = 10 .^ range(log10(1e-8),log10(1e-4),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wr_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_r=>i,))
    br = get_br(rtc_mod, params, init_rtc, 1.5)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


wr_rtcb1, wr_rtcb2, wr_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))")
wr_rtcb1a, wr_rtcb2a, wr_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>r</sub> = $(round.(wr_range[2]; sigdigits=2))")
wr_rtcb1b, wr_rtcb2b, wr_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))")

wr_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))", legendgroup="3")
wr_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))", legendgroup="3")
wr = plot([wr_rtcb_nonbs1, wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_rtcb_nonbs],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

wr = plot([wr_rtcb1, wr_rtcb2, wr_rtcb3, wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_rtcb_nonbs],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

PlotlyJS.savefig(wr, "/home/holliehindley/phd/may23_rtc/paper_plots/wr_rtcb.svg")



wab_range = 10 .^range(log10(2e-6),log10(2e-5),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wab_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_ab=>i,))
    br = get_br(rtc_mod, params, rtc_init, 1.5)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wab_rtcb1, wab_rtcb2, wab_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb, "1", "dda0ddff", "ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))")
wab_rtcb1a, wab_rtcb2a, wab_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb, "2", "ba55d3ff", "ω<sub>ab</sub> = $(round.(wab_range[2]; sigdigits=2))")
wab_rtcb1b, wab_rtcb2b, wab_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb, "3", "800080ff", "ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))")

wab_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))", legendgroup="3")

wab = plot([wab_rtcb1, wab_rtcb2, wab_rtcb3,wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_rtcb_nonbs],#wab_rtcb_nonbs],
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
    br = get_br(rtc_mod, params, rtc_init, 1.5)
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


lam_range = range(0.01,stop=0.02,length=3)

bfs=[]; dfs=[];
for i in ProgressBar(lam_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod, params, rtc_init, 1.5)
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





# looking at difference between varying wab and wr - they have the SAME effect 
wr_range = 10 .^ range(log10(1e-8),log10(1e-4),length=100)
wab_range = 10 .^range(log10(1e-7),log10(1e-3),length=100)
wr_range = range(1e-8,1e-4,length=50)
wab_range = range(1e-7,1e-3,length=50)
function change_w(param, w_range, params_rtc)
    rmr=[]
    rmb=[]
    rtcr=[]
    rtcb=[]
    rh=[]
    fas=[]
    for i in w_range
        copyparams = deepcopy(params_rtc)
        copyparams[param] = i
        solu = sol(rtc_model, init_rtc, tspan, copyparams)
        rm_r_ss = get_ssval(solu, :rm_r, species_rtc)
        rm_b_ss = get_ssval(solu, :rm_b, species_rtc)
        rtcr_ss = get_ssval(solu, :rtcr, species_rtc)
        rtcb_ss = get_ssval(solu, :rtcb, species_rtc)
        rh_ss = get_ssval(solu, :rh, species_rtc)
        rt = get_ssval(solu, :rt, species_rtc)
        alpha = rt/kr # unitless
        fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
        push!(rmr, rt)
        push!(rmb, rm_b_ss)
        push!(rtcr, rtcr_ss)
        push!(rtcb, rtcb_ss)
        push!(rh, rh_ss)
        push!(fas, fa)
    end
    return rmr, rmb, rtcr, rtcb, rh, fas
end

real_wab = (scatter(x=[ω_ab],y=[0], legendgroup="4", marker=attr(size=10, color=:red)))
real_wr = (scatter(x=[ω_r],y=[0], legendgroup="4", marker=attr(size=10, color=:red)))

params_rtc.kdam = 0
rmr, rmb, rtcr, rtcb, rh, fas = change_w(:ω_r, wr_range, params_rtc)
wr_rmr = (scatter(x=wr_range,y=rmr, name="kdam = 0", legendgroup="1"))#, Layout(xaxis_title="wr",yaxis_title="mRNA rtcR"))
wr_rmb = (scatter(x=wr_range,y=rmb, name="kdam = 0", legendgroup="1"))#, Layout(xaxis_title="wr",yaxis_title="mRNA rtcB"))
wr_rtcr = (scatter(x=wr_range,y=rtcr, name="kdam = 0", legendgroup="1"))#, Layout(xaxis_title="wr",yaxis_title="RtcR"))
wr_rtcb = (scatter(x=wr_range,y=rtcb, name="kdam = 0", legendgroup="1"))#, Layout(xaxis_title="wr",yaxis_title="RtcB"))
wr_rh = (scatter(x=wr_range, y=rh, name="kdam = 0", legendgroup="1"))
wr_fa = scatter(x=wr_range, y=fas, name="kdam = 0", legendgroup="1")


rmr_wab, rmb_wab, rtcr_wab, rtcb_wab, rh_wab, fas_wab = change_w(:ω_ab, wab_range, params_rtc)
wab_rmr = (scatter(x=wab_range,y=rmr_wab, name="kdam = 0", legendgroup="1"))
wab_rmb = (scatter(x=wab_range,y=rmb_wab, name="kdam = 0", legendgroup="1"))
wab_rtcr = (scatter(x=wab_range,y=rtcr_wab, name="kdam = 0", legendgroup="1"))
wab_rtcb = (scatter(x=wab_range,y=rtcb_wab, name="kdam = 0", legendgroup="1"))
wab_rh = (scatter(x=wab_range, y=rh_wab, name="kdam = 0", legendgroup="1"))
wab_fa = (scatter(x=wab_range, y=fas_wab, name="kdam = 0", legendgroup="1"))

# p_nodam = [wr_rmr wab_rmr; wr_rmb wab_rmb; wr_rtcr wab_rtcr; wr_rtcb wab_rtcb; wr_rh wab_rh]


params_rtc.kdam = 1
rmr1, rmb1, rtcr1, rtcb1, rh1, fas1 = change_w(:ω_r, wr_range, params_rtc)
wr_rmr1 = (scatter(x=wr_range,y=rmr1, name="kdam = 1", legendgroup="2"));
wr_rmb1 = (scatter(x=wr_range,y=rmb1, name="kdam = 1", legendgroup="2"));
wr_rtcr1 = (scatter(x=wr_range,y=rtcr1, name="kdam = 1", legendgroup="2"));
wr_rtcb1 = (scatter(x=wr_range,y=rtcb1, name="kdam = 1", legendgroup="2"));
wr_rh1 = (scatter(x=wr_range, y=rh1, name="kdam = 1", legendgroup="2"));
wr_fa1 = scatter(x=wr_range, y=fas1, name="kdam = 1", legendgroup="2")

rmr_wab1, rmb_wab1, rtcr_wab1, rtcb_wab1, rh_wab1, fas_wab1 = change_w(:ω_ab, wab_range, params_rtc)
wab_rmr1 = (scatter(x=wab_range,y=rmr_wab1, name="kdam = 1", legendgroup="2"))#, Layout(xaxis_title="wab",yaxis_title="mRNA rtcR"));
wab_rmb1 = (scatter(x=wab_range,y=rmb_wab1, name="kdam = 1", legendgroup="2"))#, Layout(xaxis_title="wab",yaxis_title="mRNA rtcB"));
wab_rtcr1 = (scatter(x=wab_range,y=rtcr_wab1, name="kdam = 1", legendgroup="2"))#, Layout(xaxis_title="wab",yaxis_title="RtcR"));
wab_rtcb1 = (scatter(x=wab_range,y=rtcb_wab1, name="kdam = 1", legendgroup="2"))#, Layout(xaxis_title="wab",yaxis_title="RtcB"));
wab_rh1 = (scatter(x=wab_range, y=rh_wab1, name="kdam = 1", legendgroup="2"))#);
wab_fas1 = (scatter(x=wab_range, y=fas_wab1, name="kdam = 1", legendgroup="2"))#);

# p_nodam = [wr_rmr1 wab_rmr1; wr_rmb1 wab_rmb1; wr_rtcr1 wab_rtcr1; wr_rtcb1 wab_rtcb1; wr_rh1 wab_rh1]

params_rtc.kdam = 5
rmr2, rmb2, rtcr2, rtcb2, rh2, fas2 = change_w(:ω_r, wr_range, params_rtc)
wr_rmr2 = (scatter(x=wr_range,y=rmr2, name="kdam = 5", legendgroup="3"));
wr_rmb2 = (scatter(x=wr_range,y=rmb2, name="kdam = 5", legendgroup="3"));
wr_rtcr2 = (scatter(x=wr_range,y=rtcr2, name="kdam = 5", legendgroup="3"));
wr_rtcb2 = (scatter(x=wr_range,y=rtcb2, name="kdam = 5", legendgroup="3"));
wr_rh2 = (scatter(x=wr_range, y=rh2, name="kdam = 5", legendgroup="3"));
wr_fa2 = scatter(x=wr_range, y=fas2, name="kdam = 5", legendgroup="3")

rmr_wab2, rmb_wab2, rtcr_wab2, rtcb_wab2, rh_wab2, fas_wab2 = change_w(:ω_ab, wab_range, params_rtc)
wab_rmr2 = (scatter(x=wab_range,y=rmr_wab2, name="kdam = 5", legendgroup="3"))#, Layout(xaxis_title="wab",yaxis_title="mRNA rtcR"));
wab_rmb2 = (scatter(x=wab_range,y=rmb_wab2, name="kdam = 5", legendgroup="3"))#, Layout(xaxis_title="wab",yaxis_title="mRNA rtcB"));
wab_rtcr2 = (scatter(x=wab_range,y=rtcr_wab2, name="kdam = 5", legendgroup="3"))#, Layout(xaxis_title="wab",yaxis_title="RtcR"));
wab_rtcb2 = (scatter(x=wab_range,y=rtcb_wab2, name="kdam = 5", legendgroup="3"))#, Layout(xaxis_title="wab",yaxis_title="RtcB"));
wab_rh2 = (scatter(x=wab_range, y=rh_wab2, name="kdam = 5", legendgroup="3"));
wab_fas2 = (scatter(x=wab_range, y=fas_wab2, name="kdam = 5", legendgroup="3"));

# p_nodam = [wr_rmr1 wab_rmr1; wr_rmb1 wab_rmb1; wr_rtcr1 wab_rtcr1; wr_rtcb1 wab_rtcb1; wr_rh1 wab_rh1]
[plot([wr_rmr, wr_rmr1, wr_rmr2,real_wr]) plot([wab_rmr, wab_rmr1, wab_rmr2,real_wab]);
plot([wr_rmb, wr_rmb1, wr_rmb2,real_wr]) plot([wab_rmb, wab_rmb1, wab_rmb2,real_wab]);
plot([wr_rtcr, wr_rtcr1, wr_rtcr2,real_wr]) plot([wab_rtcr, wab_rtcr1, wab_rtcr2,real_wab]);
plot([wr_rtcb, wr_rtcb1, wr_rtcb2,real_wr]) plot([wab_rtcb, wab_rtcb1, wab_rtcb2,real_wab]);
plot([wr_rh, wr_rh1, wr_rh2,real_wr]) plot([wab_rh, wab_rh1, wab_rh2,real_wab]);]

p = [plot([wr_rmr, wr_rmr1, wr_rmr2,real_wr], Layout(xaxis_type="log",title="wr",yaxis_title="rm_r",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_rmr, wab_rmr1, wab_rmr2,real_wab], Layout(xaxis_type="log", title="wab",xaxis_tickformat=".1e",xaxis_tickangle=75));
plot([wr_rmb, wr_rmb1, wr_rmb2,real_wr], Layout(xaxis_type="log",yaxis_title="rm_b",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_rmb, wab_rmb1, wab_rmb2,real_wab], Layout(xaxis_type="log",xaxis_tickformat=".1e",xaxis_tickangle=75));
plot([wr_rtcr, wr_rtcr1, wr_rtcr2,real_wr], Layout(xaxis_type="log",yaxis_title="RtcR",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_rtcr, wab_rtcr1, wab_rtcr2,real_wab], Layout(xaxis_type="log",xaxis_tickformat=".1e",xaxis_tickangle=75));
plot([wr_rtcb, wr_rtcb1, wr_rtcb2,real_wr], Layout(xaxis_type="log",yaxis_title="RtcB",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_rtcb, wab_rtcb1, wab_rtcb2,real_wab], Layout(xaxis_type="log",xaxis_tickformat=".1e",xaxis_tickangle=75));
plot([wr_rh, wr_rh1, wr_rh2,real_wr], Layout(xaxis_type="log", yaxis_title="Rh",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_rh, wab_rh1, wab_rh2,real_wab], Layout(xaxis_type="log",xaxis_tickformat=".1e",xaxis_tickangle=75));
plot([wr_fa, wr_fa1, wr_fa2,real_wr], Layout(xaxis_type="log", yaxis_title="fa",xaxis_tickformat=".1e",xaxis_tickangle=75)) plot([wab_fa, wab_fas1, wab_fas2,real_wab], Layout(xaxis_type="log",xaxis_tickformat=".1e",xaxis_tickangle=75));]



open("./example.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end

wr_range = 10 .^ range(log10(2e-7),log10(2e-6),length=3)
wr_range = 10 .^ range(log10(1e-8),log10(1e-4),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wr_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_r=>i,))
    br = get_br(rtc_mod, params, rtc_init, 1.5)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end


wr_rtcb1a, wr_rtcb2a, wr_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rm_r, "2", "ba55d3ff", "ω<sub>r</sub> = $(round.(wr_range[2]; sigdigits=2))")

wr_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rm_r,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[3]; sigdigits=2))", legendgroup="3")

wr_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rm_r,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>r</sub> = $(round.(wr_range[1]; sigdigits=2))", legendgroup="3")

wr = plot([wr_rtcb_nonbs1, wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_rtcb_nonbs],#wr_rtcb1b, wr_rtcb2b, wr_rtcb3b],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

PlotlyJS.savefig(wr, "/home/holliehindley/phd/may23_rtc/paper_plots/wr_rtcb.svg")



wab_range = 10 .^range(log10(2e-6),log10(2e-5),length=3)
wab_range = 10 .^range(log10(1e-7),log10(1e-3),length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wab_range)
    copyparams = deepcopy(params_bf)
    params = merge(copyparams, (:ω_ab=>i,))
    br = get_br(rtc_mod, params, rtc_init, 1.5)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wab_rtcb1, wab_rtcb2, wab_rtcb3 = plot_rtc_bf(dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rm_b, "1", "dda0ddff", "ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))")
wab_rtcb1a, wab_rtcb2a, wab_rtcb3a = plot_rtc_bf(dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rm_b, "2", "ba55d3ff", "ω<sub>ab</sub> = $(round.(wab_range[2]; sigdigits=2))")
wab_rtcb1b, wab_rtcb2b, wab_rtcb3b = plot_rtc_bf(dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rm_b, "3", "800080ff", "ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))")

wab_rtcb_nonbs = scatter(x=dfs[3].kdam,y=dfs[3].rm_b,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[3]; sigdigits=2))", legendgroup="3")

wab_rtcb_nonbs1 = scatter(x=dfs[1].kdam,y=dfs[1].rm_b,showlegend=true,line=attr(width=6.5, color="#c1c1c1ff"), name="ω<sub>ab</sub> = $(round.(wab_range[1]; sigdigits=2))", legendgroup="3")


wab = plot([wab_rtcb_nonbs1,wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_rtcb_nonbs],#wab_rtcb_nonbs],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


