using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/init_switch_funcs.jl")

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





params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

br = get_br(rtc_mod, params2, initial, 3.)



bf0 = bf_point_df(br)
df0 = create_br_df(br)
kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]




params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

kdam_range_onoff = range(df0.kdam[kdam02]+0.0001*df0.kdam[kdam02], df0.kdam[kdam01]-0.00008*df0.kdam[kdam01], length=100)
tspan=(0,1e9)

svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

for kdam_val in ProgressBar(kdam_range_onoff)
    psm = deepcopy(params1)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params2, initial)
    # @show psm
    
    n = 600; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, upper_ranges, branches1.ss_val_on,:rh,l,psm,9, all_species)
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

CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis_orig_model/init_switch/on_off/data/PAPERswitch_vals_NEW.csv", svals_onoff)

rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color="#005356ff"), showlegend=false, legendgroup="1")#, fill="tozeroy")
rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color="#f04e53ff"),showlegend=false, legendgroup="1")
bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")

rtcb_onoff = scatter(x=kdam_range_onoff, y=svals_onoff.rtcb, name="switch point", showlegend=false, line=attr(color="#9f9f9fff", dash="dot"))#, fill="tozeroy")

plot([rtcb1,rtcb2,rtcb3,bf_rtcb,rtcb_onoff])
