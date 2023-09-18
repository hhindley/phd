using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf

using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


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

br = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

kdam_range_onoff = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=10)
tspan=(0,1e9)

branches1 = setup_ssvals_from_bfkit(rtc_mod, 1, params_for_ssval_setup)
branches1.ss_val_off
initial


svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

for kdam_val in ProgressBar(kdam_range_onoff)
    psm = deepcopy(params1)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params_for_ssval_setup)
    # @show psm
    
    n = 600; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, upper_ranges, branches1.ss_val_on,:rh,l,psm,9)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l, 9)
    inds = get_switch_ind(binary, l)
    vals = get_switch_vals(inds, init_vals)
    push!(svals_onoff.rm_a, vals[1])
    push!(svals_onoff.rtca, vals[2])
    push!(svals_onoff.rm_b, vals[3])
    push!(svals_onoff.rtcb, vals[4])
    push!(svals_onoff.rm_r, vals[5])
    push!(svals_onoff.rtcr, vals[6])
    push!(svals_onoff.rh, vals[7])
    push!(svals_onoff.rd, vals[8])
    push!(svals_onoff.rt, vals[9])
end

svals_onoff

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")

plot([rtcb1,rtcb2,rtcb3,bf_rtcb,rtcb_onoff])

svals_onoff.rtcb
((collect(kdam_range_onoff)[7]-df.kdam[kdam2])/(df.kdam[kdam1] - df.kdam[kdam2]))*100

65.3032


atp_range = range(2500,stop=4000,length=3)
all_res = []
kdam_ranges = []
on_states = []
off_states = []
mids = []
bfs = []
for i in ProgressBar(atp_range)
    # svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    svals_onoff = []
    psm = deepcopy(params1)
    psm.atp = i
    param_br = deepcopy(params_for_ssval_setup)
    param_br = merge(param_br, (:atp=>i,))
    br = get_br(rtc_mod, param_br, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    kdam_range = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=20)
    push!(kdam_ranges, kdam_range)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    push!(on_states, rtcb1)
    push!(off_states, rtcb3)
    push!(mids, rtcb2)
    push!(bfs, bf_rtcb)
    for kdam_val in (kdam_range)
        # psm = deepcopy(params1)
        psm.kdam = kdam_val
        branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup)
        # @show psm
        
        n = 600; l = 1000;
        upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
        all, init_vals = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_on,:rh,l,psm)
        binary = upper_or_lower(all, branches1.ss_val_off[7], l)
        inds = get_switch_ind(binary, l)
        vals = get_switch_vals(inds, init_vals)
        push!(svals_onoff, vals[4])

    end
    push!(all_res, svals_onoff)
end
all_res
kdam_ranges[5]


a1 = scatter(x=kdam_ranges[1], y=all_res[1])
a2 = scatter(x=kdam_ranges[2], y=all_res[2])
a3 = scatter(x=kdam_ranges[3], y=all_res[3])
a4 = scatter(x=kdam_ranges[4], y=all_res[4])
a5 = scatter(x=kdam_ranges[5], y=all_res[5])

plot([a1,a2,a3,a4,a5])
plot([a1,a2,a3,a4,a5, on_states[1], on_states[2], on_states[3], on_states[4], on_states[5], off_states[1], off_states[2], off_states[3], off_states[4], off_states[5], mids[1], mids[2], mids[3], mids[4], mids[5], bfs[1], bfs[2], bfs[3], bfs[4],bfs[5]])

plot([a1,a2,a3, on_states[1], on_states[2], on_states[3], off_states[1], off_states[2], off_states[3], mids[1], mids[2], mids[3], bfs[1], bfs[2], bfs[3]])

kdam_ranges[5]



br = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

plot([rtcb1,rtcb2,rtcb3,bf_rtcb])

function plot_all_bf(bf, df)
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

    
    return rtcb_onoff, rtcr_onoff, rh_onoff, rt_onoff, rd_onoff, rtca_onoff, rma_onoff, rmb_onoff, rmr_onoff
end

rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rtcr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcr, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rh_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rh, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rt_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rt, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")
rd_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rd, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rtca_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtca, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rma_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_a, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rmb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_b, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rmr_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rm_r, name="switch point",showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tonexty")
rtcb_onoff, rtcr_onoff, rh_onoff, rt_onoff, rd_onoff, rtca_onoff, rma_onoff, rmb_onoff, rmr_onoff = plot_all_bf()
# on→off plots
rtcb_onoff_plot = plot([rtcb_onoff, rtcb1, rtcb2, rtcb3, bf_rtcb], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcB - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rtcr_onoff_plot = plot([rtcr_onoff, rtcr1, rtcr2, rtcr3, bf_rtcr], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcR - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rh_onoff_plot = plot([rh_onoff, rh1, rh2, rh3, bf_rh], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rh - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))
rt_onoff_plot = plot([rt_onoff, rt1, rt2, rt3, bf_rt], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="Rt - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))

onoff_all = [rtcb_onoff_plot rtcr_onoff_plot; rh_onoff_plot rt_onoff_plot]

# open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/onoff_all.html", "w") do io
#     PlotlyBase.to_html(io, onoff_all.plot)
# end




psm = deepcopy(params1)
psm.atp = 4500
param_br = deepcopy(params_for_ssval_setup)
param_br = merge(param_br, (:atp=>4500,))


br = get_br(rtc_mod, param_br, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

kdam_range = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=20)


svals_onoff2 = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

for kdam_val in ProgressBar(kdam_range)
    psm = deepcopy(params1)
    psm.atp = 4500
    param_br = deepcopy(params_for_ssval_setup)
    param_br = merge(param_br, (:atp=>4500,))

    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, param_br)
    # @show psm
    
    n = 600; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_on,:rh,l,psm)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l)
    inds = get_switch_ind(binary, l)
    vals = get_switch_vals(inds, init_vals)
    push!(svals_onoff2.rm_a, vals[1])
    push!(svals_onoff2.rtca, vals[2])
    push!(svals_onoff2.rm_b, vals[3])
    push!(svals_onoff2.rtcb, vals[4])
    push!(svals_onoff2.rm_r, vals[5])
    push!(svals_onoff2.rtcr, vals[6])
    push!(svals_onoff2.rh, vals[7])
    push!(svals_onoff2.rd, vals[8])
    push!(svals_onoff2.rt, vals[9])
end

rtcb_onoff2 = scatter(x=kdam_range, y=svals_onoff2.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")

rtcb_onoff_plot = plot([rtcb_onoff2, rtcb1, rtcb2, rtcb3, bf_rtcb], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration (μM)",title="RtcB - on → off", yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white"))#, xaxis_range=(df.kdam[kdam2],df.kdam[kdam1])))




plot(scatter(x=kdam_range, y=svals_onoff2.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot")))#, fill="tozeroy")
svals_onoff2













kin_range = range(0,stop=0.1,length=5)
atp_range = range(500,stop=5000,length=5)
lam_range = range(0.01,stop=0.04,length=5)
wr_range = range(1e-7,stop=0.01,length=5)
wab_range = range(0.001, stop=3, length=5)
all_res = []
kdam_ranges = []
on_states = []
off_states = []
mids = []
bfs = []
for i in ProgressBar(atp_range)
    # svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    svals_onoff = []
    psm = deepcopy(params1)
    psm.atp = i
    param_br = deepcopy(params_for_ssval_setup)
    param_br = merge(param_br, (:atp=>i,))
    br = get_br(rtc_mod, param_br, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    kdam_range = range(df.kdam[kdam2]+0.001*df.kdam[kdam2], df.kdam[kdam1]-0.001*df.kdam[kdam1], length=100)
    push!(kdam_ranges, kdam_range)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    push!(on_states, rtcb1)
    push!(off_states, rtcb3)
    push!(mids, rtcb2)
    push!(bfs, bf_rtcb)
    for kdam_val in (kdam_range)
        # psm = deepcopy(params1)
        psm.kdam = kdam_val
        branches1 = setup_ssvals_from_bfkit(kdam_val, param_br)
        # @show psm
        
        n = 600; l = 1000;
        upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
        all, init_vals = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_on,:rh,l,psm)
        binary = upper_or_lower(all, branches1.ss_val_off[7], l)
        inds = get_switch_ind(binary, l)
        vals = get_switch_vals(inds, init_vals)
        push!(svals_onoff, vals[4])

    end
    push!(all_res, svals_onoff)
end


a1 = scatter(x=kdam_ranges[1], y=all_res[1], name="ATP = $(atp_range[1])")
a2 = scatter(x=kdam_ranges[2], y=all_res[2], name="ATP = $(atp_range[2])")
a3 = scatter(x=kdam_ranges[3], y=all_res[3], name="ATP = $(atp_range[3])")
a4 = scatter(x=kdam_ranges[4], y=all_res[4], name="ATP = $(atp_range[4])")
a5 = scatter(x=kdam_ranges[5], y=all_res[5], name="ATP = $(atp_range[5])")

plot([a1,a2,a3,a4,a5])
p = plot([a1,a2,a3,a4,a5, on_states[1], on_states[2], on_states[3], on_states[4], on_states[5], off_states[1], off_states[2], off_states[3], off_states[4], off_states[5], mids[1], mids[2], mids[3], mids[4], mids[5], bfs[1], bfs[2], bfs[3], bfs[4],bfs[5]],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",legend=attr(x=0.75,y=1)))

savefig(p, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/atp_bf.svg")
open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/atp_bf.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end










kin_range = range(0,stop=0.1,length=5)
atp_range = range(500,stop=5000,length=5)
lam_range = range(0.014,stop=0.02,length=5)
wr_range = range(1e-7,stop=0.01,length=5)
wab_range = range(0.001, stop=3, length=5)
all_res = []
kdam_ranges = []
on_states = []
off_states = []
mids = []
bfs = []
for i in ProgressBar(lam_range)
    # svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    svals_onoff = []
    psm = deepcopy(params1)
    psm.lam = i
    param_br = deepcopy(params_for_ssval_setup)
    param_br = merge(param_br, (:lam=>i,))
    br = get_br(rtc_mod, param_br, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    kdam_range = range(df.kdam[kdam2]+0.001*df.kdam[kdam2], df.kdam[kdam1]-0.001*df.kdam[kdam1], length=20)
    push!(kdam_ranges, kdam_range)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    push!(on_states, rtcb1)
    push!(off_states, rtcb3)
    push!(mids, rtcb2)
    push!(bfs, bf_rtcb)
    for kdam_val in (kdam_range)
        # psm = deepcopy(params1)
        psm.kdam = kdam_val
        branches1 = setup_ssvals_from_bfkit(kdam_val, param_br)
        # @show psm
        
        n = 600; l = 100;
        upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
        all, init_vals = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_on,:rh,l,psm)
        binary = upper_or_lower(all, branches1.ss_val_off[7], l)
        inds = get_switch_ind(binary, l)
        vals = get_switch_vals(inds, init_vals)
        push!(svals_onoff, vals[4])

    end
    push!(all_res, svals_onoff)
end


a1 = scatter(x=kdam_ranges[1], y=all_res[1], name="λ = $(lam_range[1])")
a2 = scatter(x=kdam_ranges[2], y=all_res[2], name="λ = $(lam_range[2])")
a3 = scatter(x=kdam_ranges[3], y=all_res[3], name="λ = $(lam_range[3])")
a4 = scatter(x=kdam_ranges[4], y=all_res[4], name="λ = $(lam_range[4])")
a5 = scatter(x=kdam_ranges[5], y=all_res[5], name="λ = $(lam_range[5])")

plot([a1,a2,a3,a4,a5])
p = plot([a1,a2,a3,a4,a5, on_states[1], on_states[2], on_states[3], on_states[4], on_states[5], off_states[1], off_states[2], off_states[3], off_states[4], off_states[5], mids[1], mids[2], mids[3], mids[4], mids[5], bfs[1], bfs[2], bfs[3], bfs[4],bfs[5]],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",legend=attr(x=0.75,y=1)))

savefig(p, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/lam_bf.svg")
open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/lam_bf.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end








kdam_range_onoff = range(0.63575,2.0175,length=250)

res=[]
ss=[]
for kdam_val in ProgressBar(kdam_range_onoff)
    psm = deepcopy(params1)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup)
    # @show psm
    n = 600; l = 1000;
    initial = deepcopy(branches1.ss_val_on)
    initial[4] = 0
    initial[6] = 0
    initial[7] = branches1.ss_val_on[7]-(0.9*branches1.ss_val_on[7])
    solu = sol(rtc_model, initial, tspan, psm)
    push!(res, get_ssval(solu, :rtcb))
    push!(ss, branches1.ss_val_on[4])

end

ss
res

br = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

points = scatter(x=kdam_range_onoff, y=res, mode="markers")

plot([rtcb1, rtcb2, rtcb3, bf_rtcb, points])

