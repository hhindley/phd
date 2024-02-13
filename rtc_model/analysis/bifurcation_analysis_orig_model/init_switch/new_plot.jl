using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, ProgressBars
include("/home/hollie_hindley/Documents/may23_rtc/functions/solving.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/set_ups.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/plotting.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/functions/sweep_params.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/rtc_orig.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/models/single_t.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/combinations_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

params_for_ssval_setup = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)
params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

kdam_range_onoff = range(0.63575,2.0175,length=250)
kdam_range_offon = range(0.636,2.0175,length=250)


svals_onoff = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/switch_vals.csv"))
svals_offon = DataFrame(CSV.File("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/switch_vals.csv"))

# # # svals
br = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)

kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rtcr1 = scatter(x=df.kdam[1:kdam1], y=df.rtcr[1:kdam1], name="RtcR", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcr2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcr[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcr3 = scatter(x=df.kdam[kdam2:end], y=df.rtcr[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rh1 = scatter(x=df.kdam[1:kdam1], y=df.rh[1:kdam1], name="Rh", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rh2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rh[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rh3 = scatter(x=df.kdam[kdam2:end], y=df.rh[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rt1 = scatter(x=df.kdam[1:kdam1], y=df.rt[1:kdam1], name="Rt", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rt2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rt[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rt3 = scatter(x=df.kdam[kdam2:end], y=df.rt[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rd1 = scatter(x=df.kdam[1:kdam1], y=df.rd[1:kdam1], name="Rd", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tonexty")
rd2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rd[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rd3 = scatter(x=df.kdam[kdam2:end], y=df.rd[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rtca1 = scatter(x=df.kdam[1:kdam1], y=df.rtca[1:kdam1], name="RtcA", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tonexty")
rtca2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtca[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtca3 = scatter(x=df.kdam[kdam2:end], y=df.rtca[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rma1 = scatter(x=df.kdam[1:kdam1], y=df.rm_a[1:kdam1], name="mRNA RtcA", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tonexty")
rma2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_a[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rma3 = scatter(x=df.kdam[kdam2:end], y=df.rm_a[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rmr1 = scatter(x=df.kdam[1:kdam1], y=df.rm_r[1:kdam1], name="mRNA RtcR", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tonexty")
rmr2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_r[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rmr3 = scatter(x=df.kdam[kdam2:end], y=df.rm_r[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

rmb1 = scatter(x=df.kdam[1:kdam1], y=df.rm_b[1:kdam1], name="mRNA RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tonexty")
rmb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_b[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rmb3 = scatter(x=df.kdam[kdam2:end], y=df.rm_b[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rtca = scatter(x=bf.kdam, y=bf.rtca, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rtcr = scatter(x=bf.kdam, y=bf.rtcr, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rma = scatter(x=bf.kdam, y=bf.rm_a, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rmb = scatter(x=bf.kdam, y=bf.rm_b, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rmr = scatter(x=bf.kdam, y=bf.rm_r, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rh = scatter(x=bf.kdam, y=bf.rh, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rd = scatter(x=bf.kdam, y=bf.rd, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
bf_rt = scatter(x=bf.kdam, y=bf.rt, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

kdam_range_onoff = range(0.63575,2.0175,length=250)
kdam_range_offon = range(0.636,2.0175,length=250)

rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rtcb_offon = scatter(x=kdam_range_offon, y=svals_offon.rtcb, name="switch point", showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rtcr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcr, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rtcr_offon = scatter(x=kdam_range_offon, y=svals_offon.rtcr, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rh_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rh, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rh_offon = scatter(x=kdam_range_offon, y=svals_offon.rh, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rt_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rt, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rt_offon = scatter(x=kdam_range_offon, y=svals_offon.rt, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rd_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rd, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rd_offon = scatter(x=kdam_range_offon, y=svals_offon.rd, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rtca_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtca, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rtca_offon = scatter(x=kdam_range_offon, y=svals_offon.rtca, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rma_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_a, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rma_offon = scatter(x=kdam_range_offon, y=svals_offon.rm_a, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rmb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_b, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rmb_offon = scatter(x=kdam_range_offon, y=svals_offon.rm_b, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

rmr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_r, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rmr_offon = scatter(x=kdam_range_offon, y=svals_offon.rm_r, name="switch point",showlegend=false, line=attr(color=:green, dash="dot"))#, fill="tonexty")

# on→off plots
rtcb_onoff_plot = plot([rtcb_onoff, rtcb1, rtcb2, rtcb3, bf_rtcb], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcB - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rtcr_onoff_plot = plot([rtcr_onoff, rtcr1, rtcr2, rtcr3, bf_rtcr], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcR - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rh_onoff_plot = plot([rh_onoff, rh1, rh2, rh3, bf_rh], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rh - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rt_onoff_plot = plot([rt_onoff, rt1, rt2, rt3, bf_rt], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rt - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))

onoff_all = [rtcb_onoff_plot rtcr_onoff_plot; rh_onoff_plot rt_onoff_plot]

open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/onoff_all.html", "w") do io
    PlotlyBase.to_html(io, onoff_all.plot)
end

# off→on plots
rtcb_offon_plot = plot([rtcb_offon, rtcb1, rtcb2, rtcb3, bf_rtcb], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcB - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rtcr_offon_plot = plot([rtcr_offon, rtcr1, rtcr2, rtcr3, bf_rtcr], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcR - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rh_offon_plot = plot([rh_offon, rh1, rh2, rh3, bf_rh], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rh - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rt_offon_plot = plot([rt_offon, rt1, rt2, rt3, bf_rt], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rt - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rd_offon_plot = plot([rd_offon, rd1, rd2, rd3, bf_rd], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rd - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rmb_offon_plot = plot([rmb_offon, rmb1, rmb2, rmb3, bf_rmb], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="mRNA RtcB - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))
rmr_offon_plot = plot([rmr_offon, rmr1, rmr2, rmr3, bf_rmr], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="mRNA RtcR - off → on", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))

offon_all = [rtcb_offon_plot rtcr_offon_plot; rmb_offon_plot rmr_offon_plot; rh_offon_plot rt_offon_plot rd_offon_plot]

open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/offon_all.html", "w") do io
    PlotlyBase.to_html(io, offon_all.plot)
end

# savefig(rtcb_offon_plot, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/rtcb_plot.svg")
# savefig(rtcr_onoff_plot, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/rtcr_plot.svg")
# savefig(rh_onoff_plot, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/rh_plot.svg")
# savefig(rt_onoff_plot, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/rt_plot.svg")


# offon normalised plot 
ons_offon=[]
offs_offon=[]
for i in kdam_range_offon
    df = setup_ssvals_from_bfkit(i, params_for_ssval_setup)
    a = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df,1))])
    select!(a, Not(:column))
    a = rename(a, :1=>:rm_a, :2=>:rtca, :3=>:rm_b, :4=>:rtcb, :5=>:rm_r, :6=>:rtcr, :7=>:rh, :8=>:rd, :9=>:rt)
    delete!(a, 1)
    on = DataFrame(a[1,:])
    off = DataFrame(a[2,:])
    push!(ons_offon, on)
    push!(offs_offon, off)
end

on_offon_df = reduce(vcat, [i for i in ons_offon])
off_offon_df = reduce(vcat, [i for i in offs_offon])

norm_res_offon=DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[])
for (on_col, off_col, sval_col, new_col) in zip(eachcol(on_offon_df[:,1:7]), eachcol(off_offon_df[:,1:7]), eachcol(svals_offon[:,1:7]), eachcol(norm_res_offon))
    for (on, off, sval) in zip(on_col, off_col, sval_col)
        push!(new_col, sval/(on-off))
    end
end

norm_res2_offon = DataFrame(rd=[],rt=[])
for (on_col, off_col, sval_col, new_col) in zip(eachcol(on_offon_df[:,8:9]), eachcol(off_offon_df[:,8:9]), eachcol(svals_offon[:,8:9]), eachcol(norm_res2_offon))
    for (on, off, sval) in zip(on_col, off_col, sval_col)
        push!(new_col, sval/(off-on))
    end
end

norm_init_offon = hcat(norm_res_offon, norm_res2_offon)

rmb_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rm_b, name="RtcB mRNA", line=attr(color=:purple), legendgroup="mRNA RtcB")
rtcb_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rtcb, name="RtcB", line=attr(color=:plum), legendgroup="RtcB")
rtcr_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rtcr, name="RtcR", line=attr(color=:green), legendgroup="RtcR")
rh_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rh, name="Rh", line=attr(color=:gold), legendgroup="Rh")
rd_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rd, name="Rd", line=attr(color=:darkorange), legendgroup="Rd")
rt_norm_offon = scatter(x=kdam_range_offon, y=norm_init_offon.rt, name="Rt", line=attr(color=:red), legendgroup="Rt")

plot([rmb_norm_offon, rtcb_norm_offon, rtcr_norm_offon, rh_norm_offon, rd_norm_offon, rt_norm_offon], Layout(yaxis_type="log"))


# onoff normalised plot 
ons_onoff=[]
offs_onoff=[]
for i in kdam_range_onoff
    df = setup_ssvals_from_bfkit(i, params_for_ssval_setup)
    a = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df,1))])
    select!(a, Not(:column))
    a = rename(a, :1=>:rm_a, :2=>:rtca, :3=>:rm_b, :4=>:rtcb, :5=>:rm_r, :6=>:rtcr, :7=>:rh, :8=>:rd, :9=>:rt)
    delete!(a, 1)
    on = DataFrame(a[1,:])
    off = DataFrame(a[2,:])
    push!(ons_onoff, on)
    push!(offs_onoff, off)
end

on_onoff_df = reduce(vcat, [i for i in ons_onoff])
off_onoff_df = reduce(vcat, [i for i in offs_onoff])

norm_res_onoff=DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[])
for (on_col, off_col, sval_col, new_col) in zip(eachcol(on_onoff_df[:,1:7]), eachcol(off_onoff_df[:,1:7]), eachcol(svals_onoff[:,1:7]), eachcol(norm_res_onoff))
    for (on, off, sval) in zip(on_col, off_col, sval_col)
        push!(new_col, sval/(on-off))
    end
end

norm_res2_onoff = DataFrame(rd=[],rt=[])
for (on_col, off_col, sval_col, new_col) in zip(eachcol(on_onoff_df[:,8:9]), eachcol(off_onoff_df[:,8:9]), eachcol(svals_onoff[:,8:9]), eachcol(norm_res2_onoff))
    for (on, off, sval) in zip(on_col, off_col, sval_col)
        push!(new_col, sval/(off-on))
    end
end

norm_init_onoff = hcat(norm_res_onoff, norm_res2_onoff)

rtcb_norm_onoff = scatter(x=kdam_range_onoff, y=norm_init_onoff.rtcb, name="RtcB", line=attr(dash="dot", color=:plum), legendgroup="RtcB")
rtcr_norm_onoff = scatter(x=kdam_range_onoff, y=norm_init_onoff.rtcr, name="RtcR", line=attr(dash="dot", color=:green), legendgroup="RtcR")
rh_norm_onoff = scatter(x=kdam_range_onoff, y=norm_init_onoff.rh, name="Rh", line=attr(dash="dot", color=:gold), legendgroup="Rh")
rt_norm_onoff = scatter(x=kdam_range_onoff, y=norm_init_onoff.rt, name="Rt", line=attr(dash="dot", color=:red), legendgroup="Rt")

plot([rtcb_norm_onoff, rtcr_norm_onoff, rh_norm_onoff, rt_norm_onoff], Layout(yaxis_type="log"))



p = plot([rmb_norm_offon, rtcb_norm_offon, rtcr_norm_offon, rh_norm_offon, rd_norm_offon, rt_norm_offon, rtcb_norm_onoff, rtcr_norm_onoff, rh_norm_onoff, rt_norm_onoff], 
Layout(yaxis_type="log", xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="Normalised difference from ssval",
yaxis=attr(showline=true, linewidth=1, linecolor=:black, mirror=true), xaxis=attr(showline=true, linewidth=1,linecolor=:black), 
xaxis_showgrid=false, yaxis_showgrid=false, plot_bgcolor=:white, title="switching from off→on (solid) and on→off (dotted)"))

open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/all_norm_diff.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
savefig(p, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/all_norm_diff.svg")
