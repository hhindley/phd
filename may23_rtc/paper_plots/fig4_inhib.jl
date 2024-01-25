using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars, QuadGK, Interpolations

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");
include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");


colours_rtcb = ["7e5c94ff", "c48fe7ff", "e4bbffff"]
colours_rtcr = ["28726dff", "46c6beff", "8ef8f1ff"]
colours_rtca = ["a1403fff", "e25a58ff", "ff9c9bff"]




rtcb_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcb, :rtcb, 10., colours_rtcb)
p_rtcb = plot([i for i in rtcb_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=24, color="black", family="sans-serif")))


rtca_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rtca, 10., colours_rtca)
p_rtca = plot([i for i in rtca_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcA (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rtcr, 10., colours_rtcr)
p_rtcr = plot([i for i in rtcr_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcR (μM)", yaxis_tickformat=".1e",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


# rtcb_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcb, :rh, 10., colours_rtcb)
# p_rtcb_rh = plot([i for i in rtcb_rh_traces],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Rh (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtca_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rh, 10., colours_rtca)
# p_rtca_rh = plot([i for i in rtca_rh_traces],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Rh (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtcr_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rh, 10., colours_rtcr)
# p_rtcr_rh = plot([i for i in rtcr_rh_traces],
# Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
# yaxis_title="Rh (μM)",
# yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
# xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


[p_rtcb p_rtcb_rh; p_rtca p_rtca_rh; p_rtcr p_rtcr_rh]


savefig(p_rtcb, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcb1.svg")
savefig(p_rtca, "/home/holliehindley/phd/may23_rtc/paper_plots/rtca.svg")
savefig(p_rtcr, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcr.svg")
savefig(p_rtcb_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcb_rh.svg")
savefig(p_rtca_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtca_rh.svg")
savefig(p_rtcr_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcr_rh.svg")



rtcb_auc = all_area_under_curve_rh(rtc_inhib_mod_rtcb, 10.)
rtcr_auc = all_area_under_curve_rh(rtc_inhib_mod_rtcr, 10.)
rtca_auc = all_area_under_curve_rh(rtc_inhib_mod_rtca, 10.)

a = 0.3
colours_rtcb_rgba = ["rgba(126,92,148,$a)", "rgba(196,143,231,$a)", "rgba(228,187,255,$a)"]
colours_rtcr_rgba = ["rgba(40,114,109,$a)", "rgba(70,198,190,$a)", "rgba(142,248,241,$a)"]
colours_rtca_rgba = ["rgba(161,64,63,$a)", "rgba(226,90,88,$a)", "rgba(255,156,155,$a)"]

rtcb_fill1 = scatter(x=rtcb_auc[1].x,y=rtcb_auc[1].f1.(rtcb_auc[1].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=("rgba(150, 150, 150, $a)"))
rtcb_fill2 = scatter(x=rtcb_auc[2].x,y=rtcb_auc[2].f1.(rtcb_auc[2].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcb_rgba[1])
rtcb_fill3 = scatter(x=rtcb_auc[3].x,y=rtcb_auc[3].f1.(rtcb_auc[3].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcb_rgba[2])
rtcb_fill4 = scatter(x=rtcb_auc[4].x,y=rtcb_auc[4].f1.(rtcb_auc[4].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcb_rgba[3])

rtca_fill1 = scatter(x=rtca_auc[1].x,y=rtca_auc[1].f1.(rtca_auc[1].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=("rgba(150, 150, 150, $a)"))
rtca_fill2 = scatter(x=rtca_auc[2].x,y=rtca_auc[2].f1.(rtca_auc[2].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtca_rgba[1])
rtca_fill3 = scatter(x=rtca_auc[3].x,y=rtca_auc[3].f1.(rtca_auc[3].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtca_rgba[2])
rtca_fill4 = scatter(x=rtca_auc[4].x,y=rtca_auc[4].f1.(rtca_auc[4].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtca_rgba[3])

rtcr_fill1 = scatter(x=rtcr_auc[1].x,y=rtcr_auc[1].f1.(rtcr_auc[1].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=("rgba(150, 150, 150, $a)"))
rtcr_fill2 = scatter(x=rtcr_auc[2].x,y=rtcr_auc[2].f1.(rtcr_auc[2].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcr_rgba[1])
rtcr_fill3 = scatter(x=rtcr_auc[3].x,y=rtcr_auc[3].f1.(rtcr_auc[3].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcr_rgba[2])
rtcr_fill4 = scatter(x=rtcr_auc[4].x,y=rtcr_auc[4].f1.(rtcr_auc[4].x), fill="tozeroy", showlegend=false, mode="none", fillcolor=colours_rtcr_rgba[3])



rtcb_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcb, :rh, 10., colours_rtcb)
push!(rtcb_rh_traces, rtcb_fill1, rtcb_fill2, rtcb_fill3, rtcb_fill4)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtca_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rh, 10., colours_rtca)
push!(rtca_rh_traces, rtca_fill1, rtca_fill2, rtca_fill3, rtca_fill4)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rh, 10., colours_rtcr)
push!(rtcr_rh_traces, rtcr_fill1, rtcr_fill2, rtcr_fill3, rtcr_fill4)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


savefig(p_rtcb_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcb_rh_plusfill.svg")
savefig(p_rtca_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtca_rh_plusfill.svg")
savefig(p_rtcr_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcr_rh_plusfill.svg")

percentage_size_rtcb = bf_size(rtc_inhib_mod_rtcb, 10.)
percentage_size_rtca = bf_size(rtc_inhib_mod_rtca, 10.)
percentage_size_rtcr = bf_size(rtc_inhib_mod_rtcr, 10.)

rtcb_dec = protein_decrease(rtc_inhib_mod_rtcb, :rtcb, 10.)
# rtcb_dec_rh = protein_decrease(rtc_inhib_mod_rtcb, :rh)
rtca_dec = protein_decrease(rtc_inhib_mod_rtca, :rtca, 10.)
# rtca_dec_rh = protein_decrease(rtc_inhib_mod_rtca, :rh)
rtcr_dec = protein_decrease(rtc_inhib_mod_rtcr, :rtcr, 10.)
# rtcr_dec_rh = protein_decrease(rtc_inhib_mod_rtcr, :rh)


av_dec_rtcb = [mean(rtcb_dec[i]) for i in range(1,3)]
# min_dec_rtcb_rh = [minimum(rtcb_dec_rh[i]) for i in range(1,3)]
# max_dec_rtcb_rh = [maximum(rtcb_dec_rh[i]) for i in range(1,3)]
# av_dec_rtcb_rh = [mean(rtcb_dec_rh[i]) for i in range(1,3)]

av_dec_rtca = [mean(rtca_dec[i]) for i in range(1,3)]
# min_dec_rtca_rh = [minimum(rtca_dec_rh[i]) for i in range(1,3)]
# max_dec_rtca_rh = [maximum(rtca_dec_rh[i]) for i in range(1,3)]
# av_dec_rtca_rh = [mean(rtca_dec_rh[i]) for i in range(1,3)]

av_dec_rtcr = [mean(rtcr_dec[i]) for i in range(1,3)]
# min_dec_rtcr_rh = [minimum(rtcr_dec_rh[i]) for i in range(1,3)]
# max_dec_rtcr_rh = [maximum(rtcr_dec_rh[i]) for i in range(1,3)]
# av_dec_rtcr_rh = [mean(rtcr_dec_rh[i]) for i in range(1,3)]

areas_rtcb = [rtcb_auc[i][:result] for i in range(1,4)]
areas_rtcr = [rtcr_auc[i][:result] for i in range(1,4)]
areas_rtca = [rtca_auc[i][:result] for i in range(1,4)]

percs_rtcb = [100*(areas_rtcb[i]/areas_rtcb[1]) for i in range(2,4)]
percs_rtcr = [100*(areas_rtcr[i]/areas_rtcr[1]) for i in range(2,4)]
percs_rtca = [100*(areas_rtca[i]/areas_rtca[1]) for i in range(2,4)]

rdf = DataFrame("data"=>["Rtc conc.","Bistability region","Growth Capacity"], 
"rtcb"=>[round(av_dec_rtcb[3],digits=2),round(percentage_size_rtcb[3],digits=2),round(percs_rtcb[3],digits=2)],
"rtcr"=>[round(av_dec_rtcr[3],digits=2),round(percentage_size_rtcr[3],digits=2),round(percs_rtcr[3],digits=2)],
"rtca"=>[round(av_dec_rtca[3],digits=2),round(percentage_size_rtca[3],digits=2),round(percs_rtca[3],digits=2)])


p = plot([bar(rdf, x=:data, y=:rtcb, text=:rtcb, textposition="auto", name=String(:rtcb), marker_color=["#7e5c94ff","#7e5c94ff","#7e5c94ff"]),
bar(rdf, x=:data, y=:rtcr, text=:rtcr, textposition="auto", name=String(:rtcr), marker_color=["#28726dff","#28726dff","#28726dff"]),
bar(rdf, x=:data, y=:rtca, text=:rtca, textposition="auto", name=String(:rtca), marker_color=["#a1403fff","#a1403fff","#a1403fff"])], 
Layout(yaxis_title="% of original",yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))


savefig(p,"/home/holliehindley/phd/may23_rtc/paper_plots/bar.svg")


