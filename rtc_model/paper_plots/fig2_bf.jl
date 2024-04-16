using Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools
using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections
# using Plots
using PlotlyJS, ProgressBars

PATH = "/home/holliehindley/phd"

include("$PATH/general_funcs/solving.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/rtc_model/parameters/rtc_params_molecs.jl")
include("$PATH/rtc_model/functions/bf_funcs/bf_funcs.jl")


br = get_br(rtc_model, ssvals_rtc, params_rtc, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
plot(scatter(x=df.kdam, y=df.rtcb, name="RtcB"))

alpha = df.rt/kr_val # unitless
fa = @. (1+alpha)^6/(L_val*((1+c_val*alpha)^6)+(1+alpha)^6)
fa1 = scatter(x=df.kdam[1:kdam1], y=fa[1:kdam1], name="Fa", line=attr(width=6.5, color=:blue), showlegend=true, legendgroup=1, yaxis="y2")#, fill="tozeroy")
fa2 = scatter(x=df.kdam[kdam1:kdam2], y=fa[kdam1:kdam2], name="", mode="lines", line=attr(width=6.5,dash="dash", color=:blue),showlegend=false, legendgroup=1, yaxis="y2")
fa3 = scatter(x=df.kdam[kdam2:end], y=fa[kdam2:end], name="", line=attr(width=6.5, color=:blue),showlegend=false, legendgroup=1, yaxis="y2")

p = plot(scatter(x=df.kdam, y=fa, name="Fa", line=attr(width=6.5, color=:blue)),
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Fraction of active RtcR (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),yaxis2=attr(overlaying="y",showline=true,linewidth=3,linecolor="black",side="right"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))#, showlegend=false))

savefig(p, "$PATH/rtc_model/paper_plots/plots/fa.svg")


rtcb1, rtcb2, rtcb3 = plot_rtc_bf(df, kdam1, kdam2, :rtcb, "1", "b693ccff", "RtcB", "y1")
rtcr1, rtcr2, rtcr3 = plot_rtc_bf(df, kdam1, kdam2, :rtcr, "2", "4ca7a2ff", "RtcR", "y1")
rtca1, rtca2, rtca3 = plot_rtc_bf(df, kdam1, kdam2, :rtca, "3", "e48080ff", "RtcA", "y1")

p = plot([rtcb1, rtcb2, rtcb3, rtcr1, rtcr2, rtcr3, rtca1, rtca2, rtca3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),yaxis2=attr(overlaying="y",showline=true,linewidth=3,linecolor="black",side="right"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))#, showlegend=false))

rh1, rh2, rh3 = plot_rtc_bf(df, kdam1, kdam2, :rh, "1", "ffd30cff", "RtcB", "y1")
rt1, rt2, rt3 = plot_rtc_bf(df, kdam1, kdam2, :rt, "2", "e96100ff", "rt", "y1")
rd1, rd2, rd3 = plot_rtc_bf(df, kdam1, kdam2, :rd, "3", "ac0606ff", "rd", "y1")

p1 = plot([rh1, rh2, rh3, rt1, rt2, rt3, rd1, rd2, rd3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

[p1_sig p1]

savefig(p, "$PATH/rtc_model/paper_plots/plots/rtc_proteins_doubley.svg")
savefig(p1, "$PATH/may23_rtc/paper_plots/ribosomes.svg")



p = plot([rtcb1, rtcb2, rtcb3, rtca1, rtca2, rtca3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),yaxis2=attr(overlaying="y",showline=true,linewidth=3,linecolor="black",side="right"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif"), showlegend=false))#, showlegend=false))

p1 = plot([rtcr1, rtcr2, rtcr3, fa1, fa2, fa3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),yaxis2=attr(overlaying="y",showline=true,linewidth=3,linecolor="black",side="right"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif"), showlegend=false, yaxis_tickangle=77, yaxis2_tickangle=90))#, showlegend=false))

p2 = plot([rh1, rh2, rh3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif"), showlegend=false))

p3 = plot([rt1, rt2, rt3, rd1, rd2, rd3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif"), showlegend=false))


savefig(p, "/home/holliehindley/phd/rtc_model/paper_plots/plots/rtc_proteinBA.svg")
savefig(p1, "/home/holliehindley/phd/rtc_model/paper_plots/plots/rtcr_fa.svg")
savefig(p2, "/home/holliehindley/phd/rtc_model/paper_plots/plots/rh.svg")
savefig(p3, "/home/holliehindley/phd/rtc_model/paper_plots/plots/rt_rd.svg")



# kdam_range1 = range(0,1.5,length=100)
# kdam_range2 = range(1.5,0,length=100)
# res = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rh, species_rtc, kdam_range1)
# res1 = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rh, species_rtc, kdam_range2)
# plot([scatter(x=kdam_range1,y=res), scatter(x=kdam_range2,y=res1)])


# p1 = plot([rh1, rh2, rh3, bf_rh],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Rh (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
# p2 = plot([rt1, rt2, rt3, bf_rt],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Rt (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
# p3 = plot([rtcb1, rtcb2, rtcb3, bf_rtcb],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="RtcB (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))
# p4 = plot([rtcr1, rtcr2, rtcr3, bf_rtcr],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="RtcR (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


# savefig(p1, "$PATHmay23_rtc/paper_plots/rh_bf.svg")
# savefig(p2, "$PATHmay23_rtc/paper_plots/rt_bf.svg")
# savefig(p3, "$PATHmay23_rtc/paper_plots/rtcb_bf.svg")
# savefig(p4, "$PATHmay23_rtc/paper_plots/rtcr_bf.svg")




# kdam_range = range(0,1.5,length=1000)
# kdam_range2 = range(1.5,0,length=1000)

# res_trna1 = numerical_bistability_analysis(rtc_model, params_rtc, init_rtc, :rh, species_rtc, kdam_range)
# res_trna2 = numerical_bistability_analysis(rtc_model, params_rtc, init_rtc, :rh, species_rtc, kdam_range2)
# ptrna1 = scatter(x=kdam_range, y=res_trna1, name="Healthy tRNA", legendgroup=3, line=attr(color=:gold,linewidth=3))
# ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=3, showlegend=false, line=attr(color=:gold,linewidth=3))

# plot([ptrna1, ptrna2])

# res=[]
# ps = deepcopy(params_rtc)
# for i in kdam_range
#     ps.kdam = i
#     solu = sol(rtc_model, init_rtc, tspan, ps)
#     push!(res, get_all_ssvals(solu, species_rtc))
# end

# df_ssvals = DataFrame(vcat(transpose(res)...), :auto)
# rename!(df_ssvals, species_rtc)

# plot([scatter(x=kdam_range, y=col, name="$(names(df_ssvals)[i])") for (col,i) in zip(eachcol(df_ssvals), range(1,length(names(df_ssvals))))])
