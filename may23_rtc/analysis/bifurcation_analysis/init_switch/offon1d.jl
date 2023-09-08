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


params_for_ssval_setup2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014)
params2 = @LArray [L, c, kr, Vmax_init, Km_init, 1, 0.0001, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
kdam_range = range(0.6359851,2.016996,length=500)
logkdam_range = 10 .^ range(log10(0.6359851),log10(2.016996),length=500)
kdam_range1 = vcat(range(0.6359851,0.643,length=200), range(0.6451,2.016996,length=10))

tspan=(0,1e9)
all_diffs=[]
all_percs=[]
all_multiples=[]
psm = deepcopy(params1)
# instab=[]
kdam_range = range(0.251777,0.3706,length=200)

for kdam_val in ProgressBar(kdam_range)
    psm = deepcopy(params2)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup2)
    # @show psm
    
    #n = 6000; l = 7500;
    n = 1; l = 1000;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)
    # @show lower_ranges[9]
    all, init_vals = get_rh_init_switch_all_ranges(lower_ranges, branches1.ss_val_off,:rh,l,psm)
    # push!(instab,unstable)
    # push!(all_diffs, upper_or_lower(all, branches1.ss_val_off, l))

    # push!(all_diffs,full_find_differences_or_percs(all,get_diffs,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
    push!(all_percs,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
    # push!(all_multiples,full_find_differences_or_percs(all,get_multiples,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
end

# df_res = create_resdf(all_diffs,kdam_range)
df_percs = create_resdf(all_percs,kdam_range)
# df_multiples = create_resdf(all_multiples,kdam_range)


# CSV.write("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/diffsNEW.csv", df_res)
CSV.write("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/percs_to100%_newkdamrange.csv", df_percs)
# CSV.write("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/foldNEW.csv", df_multiples)


# p_diff = plot([scatter(x=df_res.kdam,y=df_res.rm_a,name="rm_a"),scatter(x=df_res.kdam,y=df_res.rtca,name="rtca"),scatter(x=df_res.kdam,y=df_res.rm_b,name="rm_b"),
# scatter(x=df_res.kdam,y=df_res.rtcb,name="rtcb"),scatter(x=df_res.kdam,y=df_res.rm_r,name="rm_r"),scatter(x=df_res.kdam,y=df_res.rtcr,name="rtcr"),
# scatter(x=df_res.kdam,y=df_res.rh,name="rh"),scatter(x=df_res.kdam,y=df_res.rd,name="rd"),scatter(x=df_res.kdam,y=df_res.rt,name="rt")],
# Layout(xaxis_title="kdam",yaxis_title="difference from ssval (μM)", title="switching from off to on",
# yaxis_type="log"))#,xaxis_type="log"))

p_perc = plot([scatter(x=df_percs.kdam,y=df_percs.rm_a,name="rm_a"),scatter(x=df_percs.kdam,y=df_percs.rtca,name="rtca"),scatter(x=df_percs.kdam,y=df_percs.rm_b,name="rm_b"),
scatter(x=df_percs.kdam,y=df_percs.rtcb,name="rtcb"),scatter(x=df_percs.kdam,y=df_percs.rm_r,name="rm_r"),scatter(x=df_percs.kdam,y=df_percs.rtcr,name="rtcr"),
scatter(x=df_percs.kdam,y=df_percs.rh,name="rh"),scatter(x=df_percs.kdam,y=df_percs.rd,name="rd"),scatter(x=df_percs.kdam,y=df_percs.rt,name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from off to on",
xaxis_type="log"))

# p_fold = plot([scatter(x=df_multiples.kdam,y=df_multiples.rm_a,name="rm_a"),scatter(x=df_multiples.kdam,y=df_multiples.rtca,name="rtca"),scatter(x=df_multiples.kdam,y=df_res.rm_b,name="rm_b"),
# scatter(x=df_multiples.kdam,y=df_multiples.rtcb,name="rtcb"),scatter(x=df_multiples.kdam,y=df_multiples.rm_r,name="rm_r"),scatter(x=df_multiples.kdam,y=df_res.rtcr,name="rtcr"),
# scatter(x=df_multiples.kdam,y=df_multiples.rh,name="rh"),scatter(x=df_multiples.kdam,y=df_res.rd,name="rd"),scatter(x=df_multiples.kdam,y=df_res.rt,name="rt")],
# Layout(xaxis_title="kdam",yaxis_title="fold-change from ssval", title="switching from off to on",
# yaxis_type="log"))

# open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_fold_offon.html", "w") do io
#     PlotlyBase.to_html(io, p_fold.plot)
# end
# open("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/rh_comp.html", "w") do io
#     PlotlyBase.to_html(io, p_perc.plot)
# end
# open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_diff_offon.html", "w") do io
#     PlotlyBase.to_html(io, p_diff.plot)
# end

