using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, ProgressBars, LabelledArrays, DataFrames, PlotlyJS, ModelingToolkit, Interpolations, QuadGK

PATH = "/home/holliehindley/phd"

include("$PATH/general_funcs/solving.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/rtc_model/models/rtc_trna_model.jl")
include("$PATH/rtc_model/functions/bf_funcs/bf_funcs.jl")
include("$PATH/rtc_model/models/inhibition_models/trna_inhib_models.jl")

colours_rtcb = ["7e5c94ff", "c48fe7ff", "e4bbffff"]
colours_rtcr = ["28726dff", "46c6beff", "8ef8f1ff"]
colours_rtca = ["a1403fff", "e25a58ff", "ff9c9bff"]

rtcb_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcb_trna_inhib_model, ssvals_trna_rtcb, params_trna_inhib, :rtcb, 20., colours_rtcb, k_trna_inhib_vals)
p_rtcb = plot([i for i in rtcb_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=24, color="black", family="sans-serif")))


rtca_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtca_trna_inhib_model, ssvals_trna_rtca, params_trna_inhib, :rtca, 20., colours_rtca, k_trna_inhib_vals)
p_rtca = plot([i for i in rtca_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcr_trna_inhib_model, ssvals_trna_rtcr, params_trna_inhib, :rtcr, 20., colours_rtcr, k_trna_inhib_vals)
p_rtcr = plot([i for i in rtcr_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcR (μM)", showlegend=false,#yaxis_tickformat=".1e",
yaxis=attr(showline=true,linewidth=3,linecolor="black",tickangle=77),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


rtcb_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcb_trna_inhib_model, ssvals_trna_rtcb, params_trna_inhib, :rh, 20., colours_rtcb, k_trna_inhib_vals)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtca_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtca_trna_inhib_model, ssvals_trna_rtca, params_trna_inhib, :rh, 20., colours_rtca, k_trna_inhib_vals)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rtcr_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcr_trna_inhib_model, ssvals_trna_rtcr, params_trna_inhib, :rh, 20., colours_rtcr, k_trna_inhib_vals)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


[p_rtcb p_rtcb_rh; p_rtca p_rtca_rh; p_rtcr p_rtcr_rh]


savefig(p_rtcb, "$PATHrtc_model/paper_plots/tRNA/rtcb1.svg")
savefig(p_rtca, "$PATHrtc_model/paper_plots/tRNA/rtca1.svg")
savefig(p_rtcr, "$PATHrtc_model/paper_plots/tRNA/rtcr2.svg")
savefig(p_rtcb_rh, "$PATHrtc_model/paper_plots/tRNA/rtcb_rh.svg")
savefig(p_rtca_rh, "$PATHrtc_model/paper_plots/tRNA/rtca_rh.svg")
savefig(p_rtcr_rh, "$PATHrtc_model/paper_plots/tRNA/rtcr_rh.svg")



br = get_br(rtc_trna_model, ssvals_trna, params_trna, 20.)
bf0 = bf_point_df(br)
df0 = create_br_df(br)

s0 = bf0.kdam[1]-bf0.kdam[2]


bf, df, kdam1, kdam2 = different_levels_inhibition(rtcr_trna_inhib_model, ssvals_trna_rtcr_inhib, params_trna_inhib, k_trna_inhib_vals[1], 20.)

s1 = bf.kdam[1]-bf.kdam[2]
percentage_of_original_size = 100*(s1/s0) 


plot(scatter(x=df.kdam,y=df.rtcr))

rtcb_auc = all_area_under_curve_rh(rtcb_trna_inhib_model, params_trna_inhib, ssvals_trna_rtcb_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
rtcr_auc = all_area_under_curve_rh(rtcr_trna_inhib_model, params_trna_inhib, ssvals_trna_rtcr_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
rtca_auc = all_area_under_curve_rh(rtca_trna_inhib_model, params_trna_inhib, ssvals_trna_rtca_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)

a = 0.3
colours_rtcb_rgba = ["rgba(228,187,255,$a)", "rgba(196,143,231,$a)", "rgba(126,92,148,$a)"]
colours_rtcr_rgba = ["rgba(142,248,241,$a)", "rgba(70,198,190,$a)", "rgba(40,114,109,$a)"]
colours_rtca_rgba = ["rgba(255,156,155,$a)", "rgba(226,90,88,$a)", "rgba(161,64,63,$a)"]

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


rtcb_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcb_trna_inhib_model, ssvals_trna_rtcb_inhib, params_trna_inhib, :rh, 20., colours_rtcb, k_trna_inhib_vals)
rtca_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtca_trna_inhib_model, ssvals_trna_rtca_inhib, params_trna_inhib, :rh, 20., colours_rtca, k_trna_inhib_vals)
rtcr_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcr_trna_inhib_model, ssvals_trna_rtcr_inhib, params_trna_inhib, :rh, 20., colours_rtcr, k_trna_inhib_vals)

# rtcb_rh_traces = creating_rtc_inhib_plot(rtc_trna_model, ssvals_trna, params_trna, rtcb_trna_inhib_model, :rh, 20., colours_rtcb, k_trna_inhib_vals)
push!(rtcb_rh_traces, rtcb_fill1, rtcb_fill2, rtcb_fill3, rtcb_fill4)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtca_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rh, 1.5, colours_rtca)
push!(rtca_rh_traces, rtca_fill1, rtca_fill2, rtca_fill3, rtca_fill4)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

# rtcr_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rh, 1.5, colours_rtcr)
push!(rtcr_rh_traces, rtcr_fill1, rtcr_fill2, rtcr_fill3, rtcr_fill4)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Healthy tRNA (μM)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


savefig(p_rtcb_rh, "$PATHrtc_model/paper_plots/tRNA/rtcb_rh_plusfill.svg")
savefig(p_rtca_rh, "$PATHrtc_model/paper_plots/tRNA/rtca_rh_plusfill.svg")
savefig(p_rtcr_rh, "$PATHrtc_model/paper_plots/tRNA/rtcr_rh_plusfill.svg")

percentage_size_rtcb = bf_size(rtcb_trna_inhib_model, ssvals_trna_rtcb_inhib, params_trna_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
percentage_size_rtca = bf_size(rtca_trna_inhib_model, ssvals_trna_rtca_inhib, params_trna_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
percentage_size_rtcr = bf_size(rtcr_trna_inhib_model, ssvals_trna_rtcr_inhib, params_trna_inhib, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)

rtcb_dec = protein_decrease(rtcb_trna_inhib_model, ssvals_trna_rtcb, params_trna_inhib, :rtcb, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
# rtcb_dec_rh = protein_decrease(rtc_inhib_mod_rtcb, :rh)
rtca_dec = protein_decrease(rtca_trna_inhib_model, ssvals_trna_rtca, params_trna_inhib, :rtca, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
# rtca_dec_rh = protein_decrease(rtc_inhib_mod_rtca, :rh)
rtcr_dec = protein_decrease(rtcr_trna_inhib_model, ssvals_trna_rtcr, params_trna_inhib, :rtcr, 20., k_trna_inhib_vals, rtc_trna_model, ssvals_trna, params_trna)
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


savefig(p,"$PATHrtc_model/paper_plots/tRNA/bar.svg")


