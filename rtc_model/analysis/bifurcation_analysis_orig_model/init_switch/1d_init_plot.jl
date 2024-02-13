using PlotlyJS, CSV, DataFrames

# offon_diffs = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/diffsNEW.csv"))
# offon_percs = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/percsNEW.csv"))
# offon_fold = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/foldNEW.csv"))

# onoff_diffs = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/diffsNEW.csv"))
onoff_percs = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/percs_newkdamrange.csv"))
# onoff_fold = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/foldNEW.csv"))

offon_percs = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/percs_to100%_newkdamrange.csv"))

offon_percs.kdam[200]

function plot_traces(offon_diffs, onoff_diffs, ytitle, y2title, title)
    diff_rma_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rm_a[1:200], name="RtcA mRNA", line=attr(color=:purple), legendgroup="RtcA mRNA")
    diff_rtca_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rtca[1:200], name="RtcA", line=attr(color=:mediumpurple), legendgroup="RtcA")
    diff_rmb_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rm_b[1:200], name="RtcB mRNA", line=attr(color=:purple), legendgroup="RtcB mRNA")
    diff_rtcb_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rtcb[1:200], name="RtcB", line=attr(color=:plum), legendgroup="RtcB")
    diff_rmr_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rm_r[1:200], name="RtcR mRNA", line=attr(color=:lightgreen), legendgroup="RtcR mRNA")
    diff_rtcr_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rtcr[1:200], name="RtcR", line=attr(color=:green), legendgroup="RtcR")
    diff_rh_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rh[1:200], name="Rh", line=attr(color=:gold), legendgroup="Rh")
    diff_rd_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rd[1:200], name="Rd", line=attr(color=:darkorange), legendgroup="Rd")
    diff_rt_offon = scatter(x=offon_diffs.kdam[1:200], y=offon_diffs.rt[1:200], name="Rt", line=attr(color=:red), legendgroup="Rt")

    diff_rma_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rm_a[9:end], yaxis="y2", name="RtcA mRNA", line=attr(color=:purple, dash="dash"), legendgroup="RtcA mRNA")
    diff_rtca_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rtca[9:end], yaxis="y2", name="RtcA", line=attr(color=:mediumpurple, dash="dash"), legendgroup="RtcA")
    diff_rmb_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rm_b[9:end], yaxis="y2", name="RtcB mRNA", line=attr(color=:blue, dash="dash"), legendgroup="RtcB mRNA")
    diff_rmr_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rm_r[9:end], yaxis="y2", name="RtcR mRNA", line=attr(color=:purple, dash="dash"), legendgroup="RtcR mRNA")
    diff_rtcr_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rtcr[9:end], yaxis="y2", name="RtcR", line=attr(color=:green, dash="dash"), legendgroup="RtcR", showlegend=false,yaxis2=attr(overlaying="y",side="right", showline=true, mirror=true))
    diff_rh_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rh[9:end], yaxis="y2", name="Rh", line=attr(color=:gold, dash="dash"), legendgroup="Rh", showlegend=false,yaxis2=attr(overlaying="y",side="right", showline=true, mirror=true))
    diff_rd_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rd[9:end], yaxis="y2", name="Rd", line=attr(color=:darkorange, dash="dash"), legendgroup="Rd")
    diff_rt_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rt[9:end], yaxis="y2", name="Rt", line=attr(color=:red, dash="dash"), legendgroup="Rt", showlegend=false,yaxis2=attr(overlaying="y",side="right", showline=true, mirror=true))
    diff_rtcb_onoff = scatter(x=onoff_diffs.kdam[9:end], y=onoff_diffs.rtcb[9:end], yaxis="y2", name="on → off", line=attr(color=:plum, dash="dash"), legendgroup="RtcB",yaxis2=attr(overlaying="y",side="right", showline=true, mirror=true))#, showlegend=false)

    return diff_rma_offon, diff_rtca_offon, diff_rmb_offon, diff_rtcb_offon,  diff_rmr_offon, diff_rtcr_offon, diff_rh_offon, diff_rd_offon, diff_rt_offon,
    diff_rma_onoff, diff_rtca_onoff, diff_rmb_onoff, diff_rtcb_onoff, diff_rmr_onoff, diff_rtcr_onoff, diff_rh_onoff, diff_rd_onoff, diff_rt_onoff

    # return plot([diff_rmb_offon,diff_rtcb_offon,diff_rmr_offon,diff_rtcr_offon,diff_rh_offon,diff_rd_offon,diff_rt_offon,
    # diff_rmb_onoff,diff_rtcb_onoff,diff_rmr_onoff,diff_rtcr_onoff,diff_rh_onoff,diff_rd_onoff,diff_rt_onoff],
    # Layout(yaxis2=attr(overlaying="y",side="right"), yaxis_title=ytitle, yaxis2_title=y2title, xaxis_title="Damage rate (min<sup>-1</sup>)",
    # yaxis2_type="log", yaxis_type="log", 
    # yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"), xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))

    # return plot([diff_rmb_offon,diff_rtcb_offon,diff_rmr_offon,diff_rtcr_offon,diff_rh_offon,diff_rd_offon,diff_rt_offon,
    # diff_rmb_onoff,diff_rtcb_onoff,diff_rmr_onoff,diff_rtcr_onoff,diff_rh_onoff,diff_rd_onoff,diff_rt_onoff],
    # Layout(yaxis2=attr(overlaying="y",side="right"), yaxis_title=ytitle, yaxis2_title=y2title, xaxis_title="Damage rate (min<sup>-1</sup>)",
    # # yaxis2_type="log", yaxis_type="log", 
    # yaxis_range=(0,100),
    # yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"), xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))
end

# diff_plot = plot_traces(offon_diffs, onoff_diffs, "Difference from ssval (μM) - off to on", "Difference from ssval (μM) - on to off", "")
perc_plot = plot_traces(offon_percs100, onoff_percs, "% increase from ssval - off to on", "% decrease from ssval - on to off", "")
# fold_plot = plot_traces(offon_fold, onoff_fold, "Fold-change from ssval - off to on", "Fold-change from ssval - on to off", "")

diff_rma_offon, diff_rtca_offon, diff_rmb_offon, diff_rtcb_offon,  diff_rmr_offon, diff_rtcr_offon, diff_rh_offon, diff_rd_offon, diff_rt_offon,
diff_rma_onoff, diff_rtca_onoff, diff_rmb_onoff, diff_rtcb_onoff, diff_rmr_onoff, diff_rtcr_onoff, diff_rh_onoff, diff_rd_onoff, diff_rt_onoff = plot_traces(offon_percs, onoff_percs, "", "", "")



p = make_subplots(rows=1,cols=2,shared_yaxes=true, horizontal_spacing=0.02);

add_trace!(p, diff_rt_offon, row=1, col=1)
add_trace!(p, diff_rtcr_offon, row=1, col=1)
add_trace!(p, diff_rtcb_offon, row=1, col=1)

add_trace!(p, diff_rtcr_onoff, row=1, col=2)
add_trace!(p, diff_rh_onoff, row=1, col=2)
add_trace!(p, diff_rt_onoff, row=1, col=2)
add_trace!(p, diff_rtcb_onoff, row=1, col=2)


relayout!(p, xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="% ↑/↓ from steady state value for off→on/on→off",
yaxis2=attr(overlaying="y",side="right", showline=true, mirror=true),
yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"), xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",
xaxis2=attr(showline=true,linewidth=1,linecolor="black"),
)
p

# yaxis2=attr(showline=false,linewidth=1,linecolor="black", mirror=false)




savefig(p,"/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/init_switch_perc_broken.svg")

open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/init_switch_broken.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end



open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/init_switch_fold.html", "w") do io
    PlotlyBase.to_html(io, fold_plot.plot)
end
open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/init_switch_perc.html", "w") do io
    PlotlyBase.to_html(io, perc_plot.plot)
end
open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/init_switch_diff.html", "w") do io
    PlotlyBase.to_html(io, diff_plot.plot)
end
