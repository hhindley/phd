using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars, QuadGK, Interpolations

PATH = "/home/holliehindley/phd"

include("$PATH/rtc_model/models/inhibition_models/rtc_inhibition_model.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/general_funcs/solving.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/rtc_model/functions/bf_funcs/bf_funcs.jl")

colours_rtcb = ["7e5c94ff", "c48fe7ff", "e4bbffff"]
colours_rtcr = ["28726dff", "46c6beff", "8ef8f1ff"]
colours_rtca = ["a1403fff", "e25a58ff", "ff9c9bff"]

solu_inhib_rtca = sol(rtca_inhib_model, init_inhib_rtca, tspan, params_inhib)
ssvals_rtca = get_all_ssvals(solu_inhib_rtca, species_inhib)

solu_inhib_rtcb = sol(rtcb_inhib_model, init_inhib_rtcb, tspan, params_inhib)
ssvals_rtcb = get_all_ssvals(solu_inhib_rtcb, species_inhib)

solu_inhib_rtcr = sol(rtcr_inhib_model, init_inhib_rtcr, tspan, params_inhib)
ssvals_rtcr = get_all_ssvals(solu_inhib_rtcr, species_inhib)


rtcb_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtcb_inhib_model, ssvals_rtcb, params_inhib, :rtcb, 1.5, colours_rtcb, k_inhib_vals)
p_rtcb = plot([i for i in rtcb_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=24, color="black", family="sans-serif")))


rtca_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtca_inhib_model, ssvals_rtca, params_inhib, :rtca, 1.5, colours_rtca, k_inhib_vals)
p_rtca = plot([i for i in rtca_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtcr_inhib_model, ssvals_rtcr, params_inhib, :rtcr, 1.5, colours_rtcr, k_inhib_vals)
p_rtcr = plot([i for i in rtcr_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcR (μM)", showlegend=false,#yaxis_tickformat=".1e",
yaxis=attr(showline=true,linewidth=3,linecolor="black",tickangle=77),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


rtcb_rh_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtcb_inhib_model, ssvals_rtcb, params_inhib, :rh, 1.5, colours_rtcb, k_inhib_vals)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtca_rh_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtca_inhib_model, ssvals_rtca, params_inhib, :rh, 1.5, colours_rtca, k_inhib_vals)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_rh_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtcr_inhib_model, ssvals_rtcr, params_inhib, :rh, 1.5, colours_rtcr, k_inhib_vals)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


[p_rtcb p_rtcb_rh; p_rtca p_rtca_rh; p_rtcr p_rtcr_rh]


savefig(p_rtcb, "$PATHmay23_rtc/paper_plots/rtcb1.svg")
savefig(p_rtca, "$PATHmay23_rtc/paper_plots/rtca1.svg")
savefig(p_rtcr, "$PATHmay23_rtc/paper_plots/rtcr2.svg")
savefig(p_rtcb_rh, "$PATHmay23_rtc/paper_plots/rtcb_rh.svg")
savefig(p_rtca_rh, "$PATHmay23_rtc/paper_plots/rtca_rh.svg")
savefig(p_rtcr_rh, "$PATHmay23_rtc/paper_plots/rtcr_rh.svg")



rtcb_auc = all_area_under_curve_rh(rtcb_inhib_model, params_inhib, ssvals_rtcb, 1.5, k_inhib_vals)
rtcr_auc = all_area_under_curve_rh(rtcr_inhib_model, params_inhib, ssvals_rtcr, 1.5, k_inhib_vals)
rtca_auc = all_area_under_curve_rh(rtca_inhib_model, params_inhib, ssvals_rtca, 1.5, k_inhib_vals)

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



# rtcb_rh_traces = creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtc_inhib_mod_rtcb, :rh, 1.5, colours_rtcb)
push!(rtcb_rh_traces, rtcb_fill1, rtcb_fill2, rtcb_fill3, rtcb_fill4)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtca_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rh, 1.5, colours_rtca)
push!(rtca_rh_traces, rtca_fill1, rtca_fill2, rtca_fill3, rtca_fill4)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtcr_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rh, 1.5, colours_rtcr)
push!(rtcr_rh_traces, rtcr_fill1, rtcr_fill2, rtcr_fill3, rtcr_fill4)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


savefig(p_rtcb_rh, "$PATHmay23_rtc/paper_plots/rtcb_rh_plusfill.svg")
savefig(p_rtca_rh, "$PATHmay23_rtc/paper_plots/rtca_rh_plusfill.svg")
savefig(p_rtcr_rh, "$PATHmay23_rtc/paper_plots/rtcr_rh_plusfill.svg")

percentage_size_rtcb = bf_size(rtcb_inhib_model, ssvals_rtcb, 10., k_inhib_vals)
percentage_size_rtca = bf_size(rtca_inhib_model, ssvals_rtca, 10., k_inhib_vals)
percentage_size_rtcr = bf_size(rtcr_inhib_model, ssvals_rtcr, 10., k_inhib_vals)

rtcb_dec = protein_decrease(rtcb_inhib_model, ssvals_rtcb, :rtcb, 10., k_inhib_vals)
# rtcb_dec_rh = protein_decrease(rtc_inhib_mod_rtcb, :rh)
rtca_dec = protein_decrease(rtca_inhib_model, ssvals_rtca, :rtca, 10., k_inhib_vals)
# rtca_dec_rh = protein_decrease(rtc_inhib_mod_rtca, :rh)
rtcr_dec = protein_decrease(rtcr_inhib_model, ssvals_rtcr, :rtcr, 10., k_inhib_vals)
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
Layout(xaxis_tickangle=0,yaxis_title="% of original",yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",showlegend=false,font=attr(size=24, color="black", family="sans-serif")))


savefig(p,"$PATHmay23_rtc/paper_plots/bar.svg")


