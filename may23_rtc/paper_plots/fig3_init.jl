

svals_onoff = DataFrame(CSV.File("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis_orig_model/init_switch/on_off/data/PAPERswitch_vals_NEW.csv"))

br = get_br(rtc_mod, params2, initial, 3.)
df = create_br_df(br)
bf = bf_point_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
kdam_range_onoff = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=100)
rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")

rtcb_01, rtcb_02, rtcb_03, bf_rtcb0 = plot_rtc_bf(bf, df, kdam1, kdam2, :rtcb)

p_rtcb = plot([rtcb_01, rtcb_02, rtcb_03, bf_rtcb0, rtcb_onoff],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))


rtcr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcr, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")



rh_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rh, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")
rt_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rt, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")

rh1, rh2, rh3, bf_rh = plot_rtc_bf(bf0, df0, kdam01, kdam02, :rh)
rt1, rt2, rt3, bf_rt = plot_rtc_bf(bf0, df0, kdam01, kdam02, :rt)

rh_traces = [rh1, rh2, rh3, bf_rh, rh_onoff]
rt_traces = [rt1, rt2, rt3, bf_rt, rt_onoff]

p_rh = plot([i for i in rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))

p_rt = plot([i for i in rt_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rt (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

savefig(p_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rh.svg")
savefig(p_rt, "/home/holliehindley/phd/may23_rtc/paper_plots/rt.svg")







bf0 = bf_point_df(br)
df0 = create_br_df(br)
kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1, k_inhib2, inhib)

bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1a, k_inhib2, inhib)

bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1b, k_inhib2, inhib)

# result0 = trapezoidal_rule(df0.kdam, df0.rh)
# result = trapezoidal_rule(df.kdam, df.rh)
# resulta = trapezoidal_rule(dfa.kdam, dfa.rh)
# resultb = trapezoidal_rule(dfb.kdam, dfb.rh)
0.6357115
using QuadGK, Interpolations


function area_under_curve(df,kdam1,specie)
    df=Float64.(df)
    x = df.kdam[54:kdam1]
    y = df[!,specie][54:kdam1]
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end
function auc_initswitch(x, y)
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)
    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end

rtcb_auc = area_under_curve(df, kdam1, :rtcb)
rtcb_fill1 = scatter(x=rtcb_auc.x,y=rtcb_auc.f1.(rtcb_auc.x), fill="tozeroy", showlegend=false, mode="none")

kdam_range_onoff = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=100)
rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")

x_rtcb = kdam_range_onoff[62:end]
y_rtcb = svals_onoff.rtcb[62:end]

auc_initswitch_rtcb = auc_initswitch(x_rtcb, y_rtcb)

rtcb_is_fill1 = scatter(x=auc_initswitch_rtcb.x,y=auc_initswitch_rtcb.f1.(auc_initswitch_rtcb.x), fill="tozeroy", showlegend=false, mode="none")

p_rtcb = plot([rtcb_01, rtcb_02, rtcb_03, bf_rtcb0, rtcb_onoff, rtcb_fill1, rtcb_is_fill1],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))

100*(auc_initswitch_rtcb.result/rtcb_auc.result)
100-(100*(auc_initswitch_rtcb.result/rtcb_auc.result))




rtcr_01, rtcr_02, rtcr_03, bf_rtcr0 = plot_rtc_bf(bf, df, kdam1, kdam2, :rtcr)

rtcr_auc = area_under_curve(df, kdam1, :rtcr)
rtcr_fill1 = scatter(x=rtcr_auc.x,y=rtcr_auc.f1.(rtcr_auc.x), fill="tozeroy", showlegend=false, mode="none")

p_rtcr = plot([rtcr_01, rtcr_02, rtcr_03, bf_rtcr0, rtcr_onoff, rtcr_fill1],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))

kdam_range_onoff = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=100)
rtcr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcr, name="switch point", showlegend=false, line=attr(width=3, color="#9f9f9fff", dash="dot"))#, fill="tozeroy")

x_rtcr = kdam_range_onoff[96:end]
y_rtcr = svals_onoff.rtcr[96:end]
auc_initswitch_rtcr = auc_initswitch(x_rtcr, y_rtcr)

rtcr_is_fill1 = scatter(x=auc_initswitch_rtcr.x,y=auc_initswitch_rtcr.f1.(auc_initswitch_rtcr.x), fill="tozeroy", showlegend=false, mode="none")

p_rtcr = plot([rtcr_01, rtcr_02, rtcr_03, bf_rtcr0, rtcr_onoff, rtcr_is_fill1, rtcr_fill1],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))

100*(auc_initswitch_rtcr.result/rtcr_auc.result)

100-(100*(auc_initswitch_rtcr.result/rtcr_auc.result))
