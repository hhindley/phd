using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf

using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/holliehindley/phd/may23_rtc/models/rtc_inhibition_model.jl");

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


tspan = (0,1e9)
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

colors1 = [:yellow, :pink, :orange, :blue, :green, :purple, :lightblue, :darkgreen, :magenta]
params_for_ssval_setup = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)
params = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
br0 = get_br(rtc_mod, params_for_ssval_setup, initial, 3.)
bf0 = bf_point_df(br0)
df0 = create_br_df(br0)
kdam10 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam20 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

bf_rma0, bf_rtca0, bf_rmb0, bf_rtcb0, bf_rmr0, bf_rtcr0, bf_rh0, bf_rd0, bf_rt0 = bf_scatter(bf0, :black)
rma_p0, rtca_p0, rmb_p0, rtcb_p0, rmr_p0, rtcr_p0, rh_p0, rd_p0, rt_p0 = full_lines(df0, 3., colors1)

p0_p = plot([rma_p0, rtca_p0, rmb_p0, rtcb_p0, rmr_p0, rtcr_p0], Layout(title="original"))
p0_r = plot([rh_p0, rd_p0, rt_p0], Layout(title="original"))
rtca10 = scatter(x=df0.kdam[1:kdam10], y=df0.rtca[1:kdam10], name="RtcA", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtca20 = scatter(x=df0.kdam[kdam10:kdam20], y=df0.rtca[kdam10:kdam20], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtca30 = scatter(x=df0.kdam[kdam20:end], y=df0.rtca[kdam20:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtca0 = scatter(x=bf0.kdam, y=bf0.rtca, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
plot([rtca10, rtca20, rtca30, bf_rtca0])






k_inhib1 = 1
k_inhib2 = 0.0025
inhib = 0.1
# ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002
initial_i = [0,0,0,0,0,0,11.29,0,0,0]
params_for_ssval_setup_inhib = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib)

br = get_br(rtc_inhib_mod_rtca, params_for_ssval_setup_inhib, initial_i, 3.)
bf = bf_point_df_inhib(br)
df = create_br_df_inhib(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr, bf_rh, bf_rd, bf_rt = bf_scatter(bf, :black)
rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, rh_p, rd_p, rt_p = full_lines(df, 3., colors1)
plot([rtca_p, bf_rtca, rtca_p0, bf_rtca0])

rtca1 = scatter(x=df.kdam[1:kdam1], y=df.rtca[1:kdam1], name="RtcA inhib", line=attr(width=3, color=:purple), legendgroup="1")#, fill="tozeroy")
rtca2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtca[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtca3 = scatter(x=df.kdam[kdam2:end], y=df.rtca[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")

bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod_rtca, k_inhib1, k_inhib2, inhib)
rtca1, rtca2, rtca3, bf_rtca = plot_rtca_bf(bf, df, kdam1, kdam2)

k_inhib1a = 0.3
k_inhib2a = 0.0025
inhiba = 0.1
bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod_rtca, k_inhib1a, k_inhib2a, inhiba)
rtca1a, rtca2a, rtca3a, bf_rtcaa = plot_rtca_bf(bfa, dfa, kdam1a, kdam2a)

k_inhib1b = 0.5
k_inhib2b = 0.0025
inhibb = 0.1
bfa1, dfa1, kdam1a1, kdam2a1 = different_levels_inhibition(rtc_inhib_mod_rtca, k_inhib1a1, k_inhib2a1, inhiba1)
rtca1a1, rtca2a1, rtca3a1, bf_rtcaa1 = plot_rtca_bf(bfb, dfb, kdam1b, kdam2b)

# rtca_inhib = plot([rtca1, rtca2, rtca3, rtca10, rtca20, rtca30])
# savefig(rtca_inhib, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/rtca_inhib.png")

p = plot([rtca10, rtca20, rtca30, bf_rtca0, rtca1, rtca2, rtca3, bf_rtca, rtca1a, rtca2a, rtca3a, bf_rtcaa, rtca1a1, rtca2a1, rtca3a1, bf_rtcaa1],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black"),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))

savefig(p, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/on_off/rtca_inhib.svg")


plot(scatter(x=df.kdam,y=df.rtcb))

bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr, bf_rh, bf_rd, bf_rt = bf_scatter(bf, :black)
rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, rh_p, rd_p, rt_p = full_lines(df, 3, colors1)
p_rtcb_i = scatter(x=df.kdam, y=df.rtcb_i, name="RtcB_i", line=attr(width=3,color=:lightgreen))
p_p = plot([rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, p_rtcb_i, bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr], Layout(title="inhib"))
p_r = plot([rh_p, rd_p, rt_p], Layout(title="inhib"))

[p0_p p0_r; p_p p_r]
plot([scatter(x=df.kdam, y=df.rtcb, name="inhib"), scatter(x=df0.kdam, y=df0.rtcb, name="orig")])



kdam_range = range(df.kdam[kdam2]+0.01*df.kdam[kdam2], df.kdam[kdam1]-0.01*df.kdam[kdam1], length=10)

svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
for kdam_val in ProgressBar(kdam_range)
    psm = deepcopy(params_inhib)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit_inhib(rtc_inhib_mod, kdam_val, params_for_ssval_setup_inhib)
    # @show psm
    
    n = 600; l = 500;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_inhib_model_rtcb, upper_ranges, branches1.ss_val_on,:rh,l,psm,10)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l,10)
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

rtcb_onoff = scatter(x=kdam_range, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")

plot([rtcb1,rtcb2,rtcb3,bf_rtcb,rtcb_onoff])



rtcr1 = scatter(x=df.kdam[1:kdam1], y=df.rtcr[1:kdam1], name="RtcR", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcr2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcr[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcr3 = scatter(x=df.kdam[kdam2:end], y=df.rtcr[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rtcr = scatter(x=bf.kdam, y=bf.rtcr, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

rtcr_onoff = scatter(x=kdam_range, y=svals_onoff.rtcr, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")

plot([rtcr1,rtcr2,rtcr3,bf_rtcr,rtcr_onoff])


rh1 = scatter(x=df.kdam[1:kdam1], y=df.rh[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
rh2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rh[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rh3 = scatter(x=df.kdam[kdam2:end], y=df.rh[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
bf_rh = scatter(x=bf.kdam, y=bf.rh, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

rh_onoff = scatter(x=kdam_range, y=svals_onoff.rh, name="switch point", showlegend=false, line=attr(color=:red, dash="dot"))#, fill="tozeroy")

plot([rh1,rh2,rh3,bf_rh,rh_onoff])


kdam_range = range(0,3,length=100)
init_rtca = [0,0,0,0,0,0,0,0,11.29,0]
species_rtca = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtca_i]
k_inhib1 = 1
k_inhib2 = 0.025
inhib = 0.1
params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014, k_inhib1, k_inhib2, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :k_inhib1, :k_inhib2, :inhib)
res = change_param(kdam_range, :kdam, rtc_inhib_model_rtca, init_rtca, get_ssval, species_rtca, params_inhib)

plot(scatter(x=kdam_range, y=res[:rtca]))