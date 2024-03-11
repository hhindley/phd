using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS
include("$PATHmay23_rtc/functions/solving.jl"); include("$PATHmay23_rtc/functions/set_ups.jl"); include("$PATHmay23_rtc/functions/plotting.jl"); 
include("$PATHmay23_rtc/functions/sweep_params.jl"); include("$PATHmay23_rtc/models/rtc_orig.jl"); include("$PATHmay23_rtc/models/atp_lam_kin_t.jl"); 
include("$PATHmay23_rtc/models/single_t.jl"); include("$PATHmay23_rtc/models/combinations_t.jl"); 
include("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


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

# branches = setup_ssvals_from_bfkit(0.6357116)
# n=600; l=1000
# all_ranges = get_all_ranges(set_ss_range_Nssval, branches, "ss_val_off", n, l)

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
kdam_range = range(0.6359851,2.016996,length=10)
tspan=(0,1e9)
all_diffs=[]
all_percs=[]
all_multiples=[]
ps = deepcopy(params1)
instab=[]
for kdam_val in kdam_range
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val)
    @show ps
    
    n = 6000; l = 7500;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)
    # @show lower_ranges[9]
    all, init_vals, unstable = get_rh_init_switch_all_ranges(lower_ranges, branches1.ss_val_off,:rh,l,ps)
    push!(instab,unstable)
    push!(all_diffs,full_find_differences_or_percs(all,get_diffs,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
    push!(all_percs,full_find_differences_or_percs(all,get_percentages,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
    push!(all_multiples,full_find_differences_or_percs(all,get_multiples,init_vals,branches1.ss_val_off[7],l,branches1.ss_val_off,0,"off"))
end
all_diffs

df_res = create_resdf(all_diffs,kdam_range)
df_percs = create_resdf(all_percs,kdam_range)
df_multiples = create_resdf(all_multiples,kdam_range)

CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/off_on/diffs.csv", df_res)
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/off_on/percs.csv", df_percs)
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/init_switch/off_on/fold.csv", df_multiples)


p_diff = plot([scatter(x=df_res.kdam,y=df_res.rm_a,name="rm_a"),scatter(x=df_res.kdam,y=df_res.rtca,name="rtca"),scatter(x=df_res.kdam,y=df_res.rm_b,name="rm_b"),
scatter(x=df_res.kdam,y=df_res.rtcb,name="rtcb"),scatter(x=df_res.kdam,y=df_res.rm_r,name="rm_r"),scatter(x=df_res.kdam,y=df_res.rtcr,name="rtcr"),
scatter(x=df_res.kdam,y=df_res.rh,name="rh"),scatter(x=df_res.kdam,y=df_res.rd,name="rd"),scatter(x=df_res.kdam,y=df_res.rt,name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (μM)", title="switching from off to on",
yaxis_type="log"))#,xaxis_type="log"))

p_perc = plot([scatter(x=df_percs.kdam,y=df_percs.rm_a,name="rm_a"),scatter(x=df_percs.kdam,y=df_percs.rtca,name="rtca"),scatter(x=df_percs.kdam,y=df_res.rm_b,name="rm_b"),
scatter(x=df_percs.kdam,y=df_percs.rtcb,name="rtcb"),scatter(x=df_percs.kdam,y=df_percs.rm_r,name="rm_r"),scatter(x=df_percs.kdam,y=df_res.rtcr,name="rtcr"),
scatter(x=df_percs.kdam,y=df_percs.rh,name="rh"),scatter(x=df_percs.kdam,y=df_res.rd,name="rd"),scatter(x=df_percs.kdam,y=df_res.rt,name="rt")],
Layout(xaxis_title="kdam",yaxis_title="difference from ssval (%)", title="switching from off to on",
yaxis_type="log"))

p_fold = plot([scatter(x=df_multiples.kdam,y=df_multiples.rm_a,name="rm_a"),scatter(x=df_multiples.kdam,y=df_multiples.rtca,name="rtca"),scatter(x=df_multiples.kdam,y=df_res.rm_b,name="rm_b"),
scatter(x=df_multiples.kdam,y=df_multiples.rtcb,name="rtcb"),scatter(x=df_multiples.kdam,y=df_multiples.rm_r,name="rm_r"),scatter(x=df_multiples.kdam,y=df_res.rtcr,name="rtcr"),
scatter(x=df_multiples.kdam,y=df_multiples.rh,name="rh"),scatter(x=df_multiples.kdam,y=df_res.rd,name="rd"),scatter(x=df_multiples.kdam,y=df_res.rt,name="rt")],
Layout(xaxis_title="kdam",yaxis_title="fold-change from ssval", title="switching from off to on",
yaxis_type="log"))

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_fold_offon.html", "w") do io
    PlotlyBase.to_html(io, p_fold.plot)
end
open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_perc_offon.html", "w") do io
    PlotlyBase.to_html(io, p_perc.plot)
end
open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_diff_offon.html", "w") do io
    PlotlyBase.to_html(io, p_diff.plot)
end


#find difference from ss val that it takes to swap from off branch to on branch 
tspan=(0,1e9)
branches1 = setup_ssvals_from_bfkit(kdam_range[4])
params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam_range[4], ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# lower to upper - single change 
n = 300; l = 1000;
# branches1 = setup_ssvals(params1)
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches1, "ss_val_off", n, l)
all, init_vals = get_rh_init_switch_all_ranges(lower_ranges, branches1.ss_val_off,:rh,l,params1)
dfs = upper_or_lower(all, branches1.ss_val_off, l)
shared_range = set_shared_range(n,l)
diffs=full_find_differences_or_percs(all,get_diffs,init_vals,branches1)


df = DataFrame(species=all_species,diff=diffs)
p = plot(df[3:end,:], x=:species, y=:diff, kind="bar", Layout(yaxis_type="log",yaxis_title="difference from ss_val to make switch (μM)", title="from OFF to ON with kdam = 1"))



p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["OFF","ON"], title="switch from lower to upper (rh) - kdam = 1"
))


p = plot([scatter(x=init_vals.rm_a,y=dfs.rm_a, name="rm_a"),scatter(x=init_vals.rtca,y=dfs.rtca, name="rtca"),scatter(x=init_vals.rm_b,y=dfs.rm_b, name="rm_b"),
scatter(x=init_vals.rtcb,y=dfs.rtcb,name="rtcb"),scatter(x=init_vals.rm_r,y=dfs.rm_r,name="rm_r"),scatter(x=init_vals.rtcr,y=dfs.rtcr,name="rtcr"),
scatter(x=init_vals.rh,y=dfs.rh,name="rh"),scatter(x=init_vals.rd,y=dfs.rd,name="rd"),scatter(x=init_vals.rt,y=dfs.rt,name="rt"),
scatter(x=[branches1.ss_val_lower[1]],y=[0],name="rm_a0",mode="markers",line=attr(color="#636efa")),scatter(x=[branches1.ss_val_lower[2]],y=[0],name="rtca0",mode="markers",line=attr(color="#ef553b")),scatter(x=[branches1.ss_val_lower[3]],y=[0],name="rm_b0",mode="markers",line=attr(color="#00cc96")),
scatter(x=[branches1.ss_val_lower[4]],y=[0],name="rtcb0",mode="markers",line=attr(color="#ab63fa")),scatter(x=[branches1.ss_val_lower[5]],y=[0],name="rm_r0",mode="markers",line=attr(color="#ffa15a")),scatter(x=[branches1.ss_val_lower[6]],y=[0],name="rtcr0",mode="markers",line=attr(color="#19d3f3")),
scatter(x=[branches1.ss_val_lower[7]],y=[0],name="rh_0",mode="markers",line=attr(color="#ff6692")),scatter(x=[branches1.ss_val_lower[8]],y=[0],name="rd0",mode="markers",line=attr(color="#b6e880")),scatter(x=[branches1.ss_val_lower[9]],y=[0],name="rt0",mode="markers",line=attr(color="#ff97ff"))],
Layout(xaxis_title="concentration (μM)", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["OFF","ON"], title="switch from lower to upper (rh) - kdam = 1",
xaxis_type="log"
))


open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_bar.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# lower to upper - double change 
n = 7; l = 100;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches1, "ss_val_lower", n, l)

sss, rtcbs, rtcrs = double_init(branches1.ss_val_lower, lower_ranges, branches1.ss_val_upper, "rtcr", "rt", params1)
# findall(x->x==branches1.ss_val_upper[7],sss)
p = plot_binary(sss, set_shared_range(n,l), branches1.ss_val_lower, "rtcr", "rt", "lower to upper - kdam = 1")

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper1.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# lower to upper - triple change
n = 5; l = 20;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches1, "ss_val_lower", n, l)

sss, init1s, init2s, init3s = triple_init(branches1.ss_val_lower, lower_ranges, branches1.ss_val_upper, "rtcr", "rt", "rtcb")
ind = findall(x->x==upper_branch.ss_val[7],sss)
sss
pars=[]
for i in range(1,length(ind))
    push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
end
pars










params2 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.64, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches2 = setup_ssvals(params2)

# lower to upper - single change 
n = 35; l = 1000;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)
all = get_rh_init_switch_all_ranges(lower_ranges, branches2.ss_val_lower, :rh, l, params2)
dfs = upper_or_lower(all, branches2.ss_val_lower, l)
shared_range = set_shared_range(n,l)

p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from lower to upper (rh) - kdam = 0.64"
))


open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowerupper_064.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end

# lower to upper - double change 
n = 1; l = 100;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

sss, rtcbs, rtcrs = double_init(branches2.ss_val_lower, lower_ranges, branches2.ss_val_upper, "rtcr", "rt", params2)
findall(x->x==branches2.ss_val_upper[7],sss)

p = plot_binary(sss, set_shared_range(n,l), branches2.ss_val_lower, "rt", "rtcr", "lower to upper - kdam = 0.64")

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# lower to upper - triple change
n = 5; l = 20;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

sss, init1s, init2s, init3s = triple_init(lower_branch, lower_ranges, upper_branch, "rtcr", "rt", "rtcb")
ind = findall(x->x==upper_branch.ss_val[7],sss)
sss
pars=[]
for i in range(1,length(ind))
    push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
end
pars









params3 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 2, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches3 = setup_ssvals(params3)

# lower to upper - single change 
n = 1000; l = 1000;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_lower", n, l)
all = get_rh_init_switch_all_ranges(lower_ranges, branches3.ss_val_lower, :rh, l, params3)
dfs = upper_or_lower(all, "lower", l)
shared_range = set_shared_range(n,l)

p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from lower to upper (rh) - kdam = 2"
))


open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowerupper2.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# lower to upper - double change 
n = 1; l = 100;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

sss, rtcbs, rtcrs = double_init(branches2.ss_val_lower, lower_ranges, branches2.ss_val_upper, "rtcr", "rt", params2)
findall(x->x==branches2.ss_val_upper[7],sss)

p = plot_binary(sss, set_shared_range(n,l), branches2.ss_val_lower, "rt", "rtcr", "lower to upper - kdam = 0.64")

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# lower to upper - triple change
n = 5; l = 20;
lower_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_lower", n, l)

sss, init1s, init2s, init3s = triple_init(lower_branch, lower_ranges, upper_branch, "rtcr", "rt", "rtcb")
ind = findall(x->x==upper_branch.ss_val[7],sss)
sss
pars=[]
for i in range(1,length(ind))
    push!(pars, (init1s[ind[i]],init2s[ind[i]],init3s[ind[i]]))
end
pars