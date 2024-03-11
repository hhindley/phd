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
kdam_range = range(0.6359851,2.016996,length=20)
tspan=(0,1e9)
all_diffs=[]
all_percs=[]
all_multiples=[]

ps = deepcopy(params1)
allres=[]
all2=[]
instab=[]
all_diffs=[]
all_percs=[]
all_multiples=[]
for kdam_val in kdam_range
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(rtc_mod, kdam_val, params_bf)
    @show ps
    
    n = 6000; l = 1000;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[9]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_model, upper_ranges, branches1.ss_val_on,:rh,l,ps,9, all_species)
    # push!(instab,unstable)
    # push!(all2, all)
    # push!(all2, get_switch_vals(get_switch_ind(upper_or_lower(all, branches1.ss_val_off, l),l),init_vals))
    # push!(all2,upper_or_lower(all,branches1.ss_val_off,l))

    # push!(all_diffs,full_find_differences_or_percs(all,get_diffs,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_on,l,"on"))
    push!(all_percs,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_on,l,"on"))
    # push!(all_multiples,full_find_differences_or_percs(all,get_multiples,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_on,l,"on"))
end


df_res = create_resdf(all_diffs,kdam_range)
df_percs = create_resdf(all_percs,kdam_range)
df_multiples = create_resdf(all_multiples,kdam_range)

CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/on_off/diffs.csv", df_res)
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/on_off/percs.csv", df_percs)
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/on_off/fold.csv", df_multiples)


p_diff = plot([scatter(x=df_res.kdam,y=df_res.rm_a,name="rm_a"),scatter(x=df_res.kdam,y=df_res.rtca,name="rtca"),scatter(x=df_res.kdam,y=df_res.rm_b,name="rm_b"),
scatter(x=df_res.kdam,y=df_res.rtcb,name="rtcb"),scatter(x=df_res.kdam,y=df_res.rm_r,name="rm_r"),scatter(x=[df_res.kdam[end]],y=[df_res.rtcr[end]],name="rtcr"),
scatter(x=[df_res.kdam[end]],y=[df_res.rh[end]],name="rh"),scatter(x=df_res.kdam,y=df_res.rd,name="rd"),scatter(x=[df_res.kdam[end]],y=[df_res.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (μM)", title="switching from on to off",))
# yaxis_type="log"))#,xaxis_type="log"))


p_perc = plot([scatter(x=df_percs.kdam,y=df_percs.rm_a,name="rm_a"),scatter(x=df_percs.kdam,y=df_percs.rtca,name="rtca"),scatter(x=df_percs.kdam,y=df_percs.rm_b,name="rm_b"),
scatter(x=df_percs.kdam,y=df_percs.rtcb,name="rtcb"),scatter(x=df_percs.kdam,y=df_percs.rm_r,name="rm_r"),scatter(x=[df_percs.kdam[end]],y=[df_percs.rtcr[end]],name="rtcr"),
scatter(x=[df_percs.kdam[end]],y=[df_percs.rh[end]],name="rh"),scatter(x=df_percs.kdam,y=df_percs.rd,name="rd"),scatter(x=[df_percs.kdam[end]],y=[df_percs.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from on to off",))
# yaxis_type="log"))

p_fold = plot([scatter(x=df_multiples.kdam,y=df_multiples.rm_a,name="rm_a"),scatter(x=df_multiples.kdam,y=df_multiples.rtca,name="rtca"),scatter(x=df_multiples.kdam,y=df_multiples.rm_b,name="rm_b"),
scatter(x=df_multiples.kdam,y=df_multiples.rtcb,name="rtcb"),scatter(x=df_multiples.kdam,y=df_multiples.rm_r,name="rm_r"),scatter(x=[df_multiples.kdam[end]],y=[df_multiples.rtcr[end]],name="rtcr"),
scatter(x=[df_multiples.kdam[end]],y=[df_multiples.rh[end]],name="rh"),scatter(x=df_multiples.kdam,y=df_multiples.rd,name="rd"),scatter(x=[df_multiples.kdam[end]],y=[df_multiples.rt[end]],name="rt")],
Layout(xaxis_title="kdam",yaxis_title="fold-change from ssval", title="switching from on to off",))
# yaxis_type="log"))

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_fold_onoff.html", "w") do io
    PlotlyBase.to_html(io, p_fold.plot)
end
open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_perc_onoff.html", "w") do io
    PlotlyBase.to_html(io, p_perc.plot)
end
open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_diff_onoff.html", "w") do io
    PlotlyBase.to_html(io, p_diff.plot)
end




# upper to lower - single change - don't see anything
n1 = 1000; l1= 1000;
upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_upper, :rh, l1, params1)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh)"
))




# upper to lower - double change - where we see change rtcr and rtcb
n1 = 0; l1= 100;
upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches1.ss_val_upper, upper_ranges, branches1.ss_val_lower, "rtcb", "rtcr", params1)
p = plot_binary(sss, range(-100,0,length=l1), branches1.ss_val_lower, "rtcb", "rtcr", "upper to lower - rh, kdam = 1")


open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


params2 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.64, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches2 = setup_ssvals(params2)

# upper to lower - single change - don't see anything
n1 = 1000; l1= 100;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches2.ss_val_upper, :rh, l1, params2)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh)"
))


# upper to lower - double change - where we see change rtcr and rtcb
n1 = 10000; l1= 50;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches2.ss_val_upper, upper_ranges, branches2.ss_val_lower, "rtcr", "rt", params2)
findall(x->x==branches2.ss_val_lower[7],sss)
sss
p = plot_binary(sss, range(-100,0,length=l1), branches2.ss_val_lower, "rt", "rtcr", "upper to lower - rh")
sss
round(sss[1];digits=6)

branches2
branches2.ss_val_upper[7]

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end









params3 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 2, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches3 = setup_ssvals(params3)

# upper to lower - single change - don't see anything
n1 = 1000; l1= 1000;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches3.ss_val_upper, :rh, l1, params3)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh) - kdam = 2"
))

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_upperlower2.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# upper to lower - double change - where we see change rtcr and rtcb
n1 = 30; l1= 50;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches3.ss_val_upper, upper_ranges, branches3.ss_val_lower, "rtcr", "rt", params3)
findall(x->x==branches3.ss_val_lower[7],sss)
shared_range1 = set_shared_range(n1,l1)
p = plot_binary(sss, shared_range1, branches3.ss_val_lower, "rt", "rtcr", "upper to lower - rh, kdam = 2")
sss
round(sss[1];digits=6)

branches2
branches2.ss_val_upper[7]

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
