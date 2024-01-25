using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/init_switch_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl");

br = get_br(rtc_mod, params_bf, rtc_init, 10.)
bf0 = bf_point_df(br)
df0 = create_br_df(br)
kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

kdam_range_onoff = range(df0.kdam[kdam02]+0.0001*df0.kdam[kdam02], df0.kdam[kdam01]-0.00008*df0.kdam[kdam01], length=100)

svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

for kdam_val in ProgressBar(kdam_range_onoff)
    psm = deepcopy(params_rtc)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params_bf, rtc_init, 10.)
    # @show psm
    
    n = 600; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, upper_ranges, branches1.ss_val_on,:rh,l,psm,9, species_rtc)
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
