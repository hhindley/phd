using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, Printf
include("$PATHmay23_rtc/functions/solving.jl"); include("$PATHmay23_rtc/functions/set_ups.jl"); include("$PATHmay23_rtc/functions/plotting.jl"); 
include("$PATHmay23_rtc/functions/sweep_params.jl"); include("$PATHmay23_rtc/models/rtc_orig.jl"); include("$PATHmay23_rtc/models/atp_lam_kin_t.jl"); 
include("$PATHmay23_rtc/models/single_t.jl"); include("$PATHmay23_rtc/models/combinations_t.jl"); include("$PATHmay23_rtc/functions/bf_funcs/bf_funcs.jl");

include("$PATHmay23_rtc/functions/bf_funcs/init_switch_funcs.jl");


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

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
kdam_range = range(0.6359851,2.016996,length=200)
tspan=(0,1e9)

ps = deepcopy(params1)

all_percs_onoff=[]

for kdam_val in kdam_range
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params_bf, initial)
    @show ps
    
    n = 6000; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)

    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, upper_ranges, branches1.ss_val_on,:rh,l,ps,9, all_species)
    
    push!(all_percs_onoff,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_on,l,"on"))
end

df_percs_onoff = create_resdf(all_percs_onoff,kdam_range)


p_perc_onoff = plot([scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rm_a,name="rm_a"),scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rtca,name="rtca"),scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rm_b,name="rm_b"),
scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rtcb,name="rtcb"),scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rm_r,name="rm_r"),scatter(x=[df_percs_onoff.kdam[end]],y=[df_percs_onoff.rtcr[end]],name="rtcr"),
scatter(x=[df_percs_onoff.kdam[end]],y=[df_percs_onoff.rh[end]],name="rh"),scatter(x=df_percs_onoff.kdam,y=df_percs_onoff.rd,name="rd"),scatter(x=[df_percs_onoff.kdam[end]],y=[df_percs_onoff.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from on to off",))



all_percs_offon=[]

for kdam_val in kdam_range
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params_bf, initial)
    @show ps
    
    n = 6000; l = 1000;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, lower_ranges, branches1.ss_val_off, :rh, l, ps, 9, all_species)
    
    push!(all_percs_offon,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
end

df_percs_offon = create_resdf(all_percs_offon,kdam_range)


p_perc_offon = plot([scatter(x=df_percs_offon.kdam,y=df_percs_offon.rm_a,name="rm_a"),scatter(x=df_percs_offon.kdam,y=df_percs_offon.rtca,name="rtca"),scatter(x=df_percs_offon.kdam,y=df_percs_offon.rm_b,name="rm_b"),
scatter(x=df_percs_offon.kdam,y=df_percs_offon.rtcb,name="rtcb"),scatter(x=df_percs_offon.kdam,y=df_percs_offon.rm_r,name="rm_r"),scatter(x=[df_percs_offon.kdam[end]],y=[df_percs_offon.rtcr[end]],name="rtcr"),
scatter(x=[df_percs_offon.kdam[end]],y=[df_percs_offon.rh[end]],name="rh"),scatter(x=df_percs_offon.kdam,y=df_percs_offon.rd,name="rd"),scatter(x=[df_percs_offon.kdam[end]],y=[df_percs_offon.rt[end]],name="rt")],
Layout(yaxis_type="log",xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from off to on",))


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
end

# diff_plot = plot_traces(offon_diffs, onoff_diffs, "Difference from ssval (μM) - off to on", "Difference from ssval (μM) - on to off", "")
perc_plot = plot_traces(df_percs_offon, df_percs_onoff, "% increase from ssval - off to on", "% decrease from ssval - on to off", "")
# fold_plot = plot_traces(offon_fold, onoff_fold, "Fold-change from ssval - off to on", "Fold-change from ssval - on to off", "")

diff_rma_offon, diff_rtca_offon, diff_rmb_offon, diff_rtcb_offon,  diff_rmr_offon, diff_rtcr_offon, diff_rh_offon, diff_rd_offon, diff_rt_offon,
diff_rma_onoff, diff_rtca_onoff, diff_rmb_onoff, diff_rtcb_onoff, diff_rmr_onoff, diff_rtcr_onoff, diff_rh_onoff, diff_rd_onoff, diff_rt_onoff = plot_traces(offon_percs, onoff_percs, "", "", "")





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




#tRNA
rh = 11.29 #75 # conc of ribosomes in exponential phase 
thr_t = 5#30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 1

init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
kdam_range_trna = range(0.,294.39,length=20)

params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

all_percs_onoff_trna=[]

for kdam_val in kdam_range_trna
    ps = deepcopy(params_trna)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod_trna, kdam_val, params_trna2, init_trna)
    @show ps
    
    n = 6000; l = 10000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)

    all, init_vals = get_rh_init_switch_all_ranges(rtc_model_trna, upper_ranges, branches1.ss_val_on,:rh,l,ps,9, all_species)
    
    push!(all_percs_onoff_trna,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_on,l,"on"))
end

df_percs_onoff_trna = create_resdf(all_percs_onoff_trna,kdam_range)


p_perc_onoff_trna = plot([scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rm_a,name="rm_a"),scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rtca,name="rtca"),scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rm_b,name="rm_b"),
scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rtcb,name="rtcb"),scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rm_r,name="rm_r"),scatter(x=[df_percs_onoff_trna.kdam[end]],y=[df_percs_onoff_trna.rtcr[end]],name="rtcr"),
scatter(x=[df_percs_onoff_trna.kdam[end]],y=[df_percs_onoff_trna.rh[end]],name="rh"),scatter(x=df_percs_onoff_trna.kdam,y=df_percs_onoff_trna.rd,name="rd"),scatter(x=[df_percs_onoff_trna.kdam[end]],y=[df_percs_onoff_trna.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from on to off",))



all_percs_offon_trna=[]

for kdam_val in kdam_range_trna
    ps = deepcopy(params_trna)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod_trna, kdam_val, params_trna2, init_trna)
    @show ps
    
    n = 6000; l = 10000;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

    all, init_vals = get_rh_init_switch_all_ranges(rtc_model_trna, lower_ranges, branches1.ss_val_off, :rh, l, ps, 9, all_species)
    
    push!(all_percs_offon_trna,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
end

df_percs_offon_trna = create_resdf(all_percs_offon_trna,kdam_range)


p_perc_offon_trna = plot([scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rm_a,name="rm_a"),scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rtca,name="rtca"),scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rm_b,name="rm_b"),
scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rtcb,name="rtcb"),scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rm_r,name="rm_r"),scatter(x=[df_percs_offon_trna.kdam[end]],y=[df_percs_offon_trna.rtcr[end]],name="rtcr"),
scatter(x=[df_percs_offon_trna.kdam[end]],y=[df_percs_offon_trna.rh[end]],name="rh"),scatter(x=df_percs_offon_trna.kdam,y=df_percs_offon_trna.rd,name="rd"),scatter(x=[df_percs_offon_trna.kdam[end]],y=[df_percs_offon_trna.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from on to off",))
