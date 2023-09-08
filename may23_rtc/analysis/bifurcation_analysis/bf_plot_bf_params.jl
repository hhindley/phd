using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames, ProgressBars

include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")

kin_range = range(0,stop=0.1,length=10)
atp_range = range(500,stop=5000,length=3)
lam_range = range(0.001,stop=0.04,length=3)
wr_range = (range(1e-7,stop=0.01,length=10))
wab_range = range(0.001, stop=3, length=10)

function plot_different_bp_lines(param_range, params1, param, specie, kdam_max)
    p = plot();
    copyparams = deepcopy(params1)
    for i in param_range
        params = merge(copyparams, (param=>i,))
        br = get_br(rtc_mod, params, initial, kdam_max)
        plot!(p, br, vars = (:param, specie), legend=false)#, label="$param = $i", putspecialptlegend=false)
    end
    return p
end

params_for_ssval_setup = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

p1_rma = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rm_a, 3.);
p2_rma = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rm_a, 3.);
p3_rma = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rm_a, 3.);
p4_rma = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rm_a, 3.);
p5_rma = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rm_a, 3.);

# l = @layout [a b; c d; e];
# all_rma = plot(p1_rma, p2_rma, p3_rma, p4_rma, p5_rma, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rma, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rma.svg")

p1_rmr = plot_different_bp_lines(atp_range, params1, :atp, :rm_r, 3.);
p2_rmr = plot_different_bp_lines(kin_range, params1, :kin, :rm_r, 3.);
p3_rmr = plot_different_bp_lines(lam_range, params1, :lam, :rm_r, 3.);
p4_rmr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_r, 3.);
p5_rmr = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_r, 3.);

# all_rmr = plot(p1_rmr, p2_rmr, p3_rmr, p4_rmr, p5_rmr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rmr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rmr.svg")

p1_rtca = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rtca, 3.);
p2_rtca = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rtca, 3.);
p3_rtca = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rtca, 3.);
p4_rtca = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rtca, 3.);
p5_rtca = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rtca, 3.);

# all_rtca = plot(p1_rtca, p2_rtca, p3_rtca, p4_rtca, p5_rtca, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rtca, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtca.svg")

p1_rtcr = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rtcr, 3.);
p2_rtcr = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rtcr, 3.);
p3_rtcr = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rtcr, 3.);
p4_rtcr = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rtcr, 3.);
p5_rtcr = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rtcr, 3.);

# all_rtcr = plot(p1_rtcr, p2_rtcr, p3_rtcr, p4_rtcr, p5_rtcr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rtcr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtcr.svg")

p1_rh = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rh, 3.);
p2_rh = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rh, 3.);
p3_rh = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rh, 3.);
p4_rh = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rh, 3.);
p5_rh = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rh, 3.);

# all_rh = plot(p1_rh, p2_rh, p3_rh, p4_rh, p5_rh, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rh.svg")

p1_rt = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rt, 3.);
p2_rt = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rt, 3.);
p3_rt = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rt, 3.);
p4_rt = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rt, 3.);
p5_rt = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rt, 3.);

# all_rt = plot(p1_rt, p2_rt, p3_rt, p4_rt, p5_rt, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rt, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rt.svg")

p1_rd = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rd, 3.);
p2_rd = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rd, 3.);
p3_rd = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rd, 3.);
p4_rd = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rd, 3.);
p5_rd = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rd, 3.);

# all_rd = plot(p1_rd, p2_rd, p3_rd, p4_rd, p5_rd, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rd, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rd.svg")

p1_rtcb = plot_different_bp_lines(atp_range, params_for_ssval_setup, :atp, :rtcb, 3.);
p2_rtcb = plot_different_bp_lines(kin_range, params_for_ssval_setup, :kin, :rtcb, 3.);
p3_rtcb = plot_different_bp_lines(lam_range, params_for_ssval_setup, :lam, :rtcb, 3.);
p4_rtcb = plot_different_bp_lines(wab_range, params_for_ssval_setup, :ω_ab, :rtcb, 3.);
p5_rtcb = plot_different_bp_lines(wr_range, params_for_ssval_setup, :ω_r, :rtcb, 3.);

all_rtca = plot(p1_rtcb, p2_rtcb, p3_rtcb, p4_rtcb, p5_rtcb, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")


l = @layout [a b; c d; e f; h i]
atp_plots = plot(p1_rma, p1_rtca, p1_rtcb, p1_rmr, p1_rtcr, p1_rh, p1_rt, p1_rd, layout=l, size=(1200,800), plot_title="ATP = $([round.(i; digits=4) for i in atp_range])", titlefontsize=6, left_margin=4Plots.mm)
kin_plots = plot(p2_rma, p2_rtca, p1_rtcb, p2_rmr, p2_rtcr, p2_rh, p2_rt, p2_rd, layout=l, size=(1200,800), plot_title="kin = $([round.(i; digits=4) for i in kin_range])", titlefontsize=6, left_margin=4Plots.mm)
lam_plots = plot(p3_rma, p3_rtca, p1_rtcb, p3_rmr, p3_rtcr, p3_rh, p3_rt, p3_rd, layout=l, size=(1200,800), plot_title="λ = $([round.(i; digits=4) for i in lam_range])", titlefontsize=6, left_margin=4Plots.mm)
wab_plots = plot(p4_rma, p4_rtca, p1_rtcb, p4_rmr, p4_rtcr, p4_rh, p4_rt, p4_rd, layout=l, size=(1200,800), plot_title="ω_ab = $([round.(i; digits=4) for i in wab_range])", titlefontsize=6, left_margin=4Plots.mm)
wr_plots = plot(p5_rma, p5_rtca, p1_rtcb, p5_rmr, p5_rtcr, p5_rh, p5_rt, p5_rd, layout=l, size=(1200,800), plot_title="ω_r = $([round.(i; digits=4) for i in wr_range])", titlefontsize=6, left_margin=4Plots.mm)

# savefig(atp_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/atp_plots.svg")
# savefig(kin_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/kin_plots.svg")
# savefig(lam_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/lam_plots.svg")
# savefig(wab_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wab_plots.svg")
# savefig(wr_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wr_plots.svg")



kin_range = range(0,stop=0.1,length=100)
atp_range = range(300,stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)
wr_range = (range(1e-7,stop=0.01,length=100))
wab_range = range(0.001, stop=3, length=100)

function bp_points_plot(params1, kdam_max)
    df = DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    copyparams = deepcopy(params1)
    param_ranges = [atp_range, lam_range, kin_range, wab_range, wr_range]
    param = [:atp, :lam, :kin, :ω_ab, :ω_r]
    range_labels = ["ATP = 300-5500", "λ = 0.001-0.04", "kin = 0-0.1", "ω_ab = 0.001-3", "ω_r = 1e-7-0.01"]

    dfs = []
    for (range, par) in zip(param_ranges, param)
        df = DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
        for i in range
            params = merge(copyparams, (par=>i,))
            br = get_br(rtc_mod, params, initial, kdam_max)
            if length(br.specialpoint) == 2
                push!(df, (i, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end
        push!(dfs, df)
    end
    bsps = [df[df.bs .== :bp, :] for df in dfs]
    
    diffs = [(bsp.kdam1 - bsp.kdam2) for bsp in bsps]

    p=plot()
    p=[plot!(range(1, stop=100, length=length(bsps[i].param)), diffs[i], label="$(range_labels[i])", c=palette(:tab10)[i], linewidth=1.5, xformatter=Returns(""), titlefontsize=10, title="size of bistability region on increasing parameter values, kdam=$kdam_max", ylabel="size of bistability region", xlabel="param range") for i in range(1,length(range_labels))]

    p = plot!(p[1])
end


p1 = bp_points_plot(params1, 1.)
p5 = bp_points_plot(params1, 3.5)

savefig(p1, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p1_bs_area.svg")
savefig(p5, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p5_bs_area.svg")



function bp_points_plot(params1, kdam_max)
    df = DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    copyparams = deepcopy(params1)
    param_ranges = [atp_range, lam_range, kin_range, wab_range, wr_range]
    param = [:atp, :lam, :kin, :ω_ab, :ω_r]
    range_labels = ["ATP = 300-5500", "λ = 0.001-0.04", "kin = 0-0.1", "ω_ab = 0.001-3", "ω_r = 1e-7-0.01"]

    dfs = []
    for (range, par) in zip(param_ranges, param)
        df = DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
        for i in range
            params = merge(copyparams, (par=>i,))
            br = get_br(rtc_mod, params, initial, kdam_max)
            if length(br.specialpoint) == 2
                push!(df, (i, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end
        push!(dfs, df)
    end
    bsps = [df[df.bs .== :bp, :] for df in dfs]
    
    diffs = [(bsp.kdam1 - bsp.kdam2) for bsp in bsps]

    # p=plot()
    # p=[plot!(range(1, stop=100, length=length(bsps[i].param)), diffs[i], label="$(range_labels[i])", xformatter=Returns(""), titlefontsize=10, title="size of bistability region on increasing parameter values, kdam=$kdam_max", ylabel="size of bistability region", xlabel="param range") for i in range(1,length(range_labels))]

    # p = plot!(p[1])
    return bsps, diffs
end

range_labels = ["ATP = 300-5500", "λ = 0.001-0.04", "kin = 0-0.1", "ω_ab = 0.001-3", "ω_r = 1e-7-0.01"]
kdam_max = 1.

bsps1, diffs1 = bp_points_plot(params1, 1.)
bsps3, diffs3 = bp_points_plot(params1, 3.)

p=plot()
p=[plot!(range(1, stop=100, length=length(bsps1[i].param)), diffs1[i], label="$(range_labels[i])", linewidth=1.5, c=palette(:tab10)[i], xformatter=Returns(""), titlefontsize=10, title="size of bistability region on increasing parameter values, kdam=$kdam_max", ylabel="size of bistability region", xlabel="param range") for i in range(1,length(range_labels))]
p=[plot!(range(1, stop=100, length=length(bsps3[i].param)), diffs3[i], c=palette(:tab10)[i], linewidth=1.5, linestyle=:dash, xformatter=Returns(""), titlefontsize=10, label="", title="dashed: kdam = 3, solid: kdam = 1", ylabel="size of bistability region", xlabel="param range") for i in range(1,length(range_labels))]

p = plot!(p[1])

savefig(p, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/pboth_bs_area.svg")






function param_change(param_range, params1, param, specie, kdam_max)
    p = plot();
    copyparams = deepcopy(params1)
    for i in param_range
        params = merge(copyparams, (param=>i,))
        br = get_br(rtc_mod, params, initial, kdam_max)
        # plot!(p, br, vars = (:param, specie), legend=false)#, label="$param = $i", putspecialptlegend=false)
    end
    return p
end

using PlotlyJS

atp_range = range(1000,stop=5000,length=3)
lam_range = range(0.01,stop=0.04,length=3)

copyparams = deepcopy(params_for_ssval_setup)

bfs=[]; dfs=[];
for i in ProgressBar(atp_range)
    copyparams = deepcopy(params_for_ssval_setup)
    params = merge(copyparams, (:atp=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

function plot_all_bf(bf, df, color, param_val, param)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="$param = $param_val", line=attr(width=3, color=color), showlegend=true, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rtcr1 = scatter(x=df.kdam[1:kdam1], y=df.rtcr[1:kdam1], name="RtcR", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcr2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcr[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rtcr3 = scatter(x=df.kdam[kdam2:end], y=df.rtcr[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rh1 = scatter(x=df.kdam[1:kdam1], y=df.rh[1:kdam1], name="Rh", yaxis="y2", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rh2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rh[kdam1:kdam2], name="", yaxis="y2", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rh3 = scatter(x=df.kdam[kdam2:end], y=df.rh[kdam2:end], name="", yaxis="y2", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rt1 = scatter(x=df.kdam[1:kdam1], y=df.rt[1:kdam1], name="Rt", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rt2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rt[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rt3 = scatter(x=df.kdam[kdam2:end], y=df.rt[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rd1 = scatter(x=df.kdam[1:kdam1], y=df.rd[1:kdam1], name="Rd", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tonexty")
    rd2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rd[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rd3 = scatter(x=df.kdam[kdam2:end], y=df.rd[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rtca1 = scatter(x=df.kdam[1:kdam1], y=df.rtca[1:kdam1], name="RtcA", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tonexty")
    rtca2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtca[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rtca3 = scatter(x=df.kdam[kdam2:end], y=df.rtca[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rma1 = scatter(x=df.kdam[1:kdam1], y=df.rm_a[1:kdam1], name="mRNA RtcA", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tonexty")
    rma2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_a[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rma3 = scatter(x=df.kdam[kdam2:end], y=df.rm_a[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rmr1 = scatter(x=df.kdam[1:kdam1], y=df.rm_r[1:kdam1], name="mRNA RtcR", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tonexty")
    rmr2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_r[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rmr3 = scatter(x=df.kdam[kdam2:end], y=df.rm_r[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    rmb1 = scatter(x=df.kdam[1:kdam1], y=df.rm_b[1:kdam1], name="mRNA RtcB", line=attr(width=3, color=color), showlegend=false, legendgroup="1")#, fill="tonexty")
    rmb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rm_b[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=color),showlegend=false, legendgroup="1")
    rmb3 = scatter(x=df.kdam[kdam2:end], y=df.rm_b[kdam2:end], name="", line=attr(width=3, color=color),showlegend=false, legendgroup="1")

    bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rtca = scatter(x=bf.kdam, y=bf.rtca, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rtcr = scatter(x=bf.kdam, y=bf.rtcr, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rma = scatter(x=bf.kdam, y=bf.rm_a, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rmb = scatter(x=bf.kdam, y=bf.rm_b, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rmr = scatter(x=bf.kdam, y=bf.rm_r, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rh = scatter(x=bf.kdam, y=bf.rh, mode="markers", name="Bifurcation point", yaxis="y2", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rd = scatter(x=bf.kdam, y=bf.rd, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")
    bf_rt = scatter(x=bf.kdam, y=bf.rt, mode="markers", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="1")

    
    return rtcb1, rtcb2, rtcb3, rtcr1, rtcr2, rtcr3, rh1, rh2, rh3, rt1, rt2, rt3, rd1, rd2, rd3, rtca1, rtca2, rtca3, rma1, rma2, rma3, rmr1, rmr2, rmr3, rmb1, rmb2, rmb3, bf_rtcb, bf_rtca, bf_rtcr, bf_rma, bf_rmb, bf_rmr, bf_rh, bf_rd, bf_rt

end


rtcb1_1, rtcb2_1, rtcb3_1, rtcr1_1, rtcr2_1, rtcr3_1, rh1_1, rh2_1, rh3_1, rt1_1, rt2_1, rt3_1, rd1_1, rd2_1, rd3_1, rtca1_1, rtca2_1, rtca3_1, rma1_1, rma2_1, rma3_1, rmr1_1, rmr2_1, rmr3_1, rmb1_1, rmb2_1, rmb3_1, bf_rtcb_1, bf_rtca_1, bf_rtcr_1, bf_rma_1, bf_rmb_1, bf_rmr_1, bf_rh_1, bf_rd_1, bf_rt_1 = plot_all_bf(bfs[1],dfs[1], :plum, atp_range[1], "ATP")
rtcb1_2, rtcb2_2, rtcb3_2, rtcr1_2, rtcr2_2, rtcr3_2, rh1_2, rh2_2, rh3_2, rt1_2, rt2_2, rt3_2, rd1_2, rd2_2, rd3_2, rtca1_2, rtca2_2, rtca3_2, rma1_2, rma2_2, rma3_2, rmr1_2, rmr2_2, rmr3_2, rmb1_2, rmb2_2, rmb3_2, bf_rtcb_2, bf_rtca_2, bf_rtcr_2, bf_rma_2, bf_rmb_2, bf_rmr_2, bf_rh_2, bf_rd_2, bf_rt_2 = plot_all_bf(bfs[2],dfs[2], :mediumorchid, atp_range[2], "ATP")
rtcb1_3, rtcb2_3, rtcb3_3, rtcr1_3, rtcr2_3, rtcr3_3, rh1_3, rh2_3, rh3_3, rt1_3, rt2_3, rt3_3, rd1_3, rd2_3, rd3_3, rtca1_3, rtca2_3, rtca3_3, rma1_3, rma2_3, rma3_3, rmr1_3, rmr2_3, rmr3_3, rmb1_3, rmb2_3, rmb3_3, bf_rtcb_3, bf_rtca_3, bf_rtcr_3, bf_rma_3, bf_rmb_3, bf_rmr_3, bf_rh_3, bf_rd_3, bf_rt_3 = plot_all_bf(bfs[3],dfs[3], :purple, atp_range[3], "ATP")


rtcbp = plot([rtcb1_1, rtcb2_1, rtcb3_1, bf_rtcb_1, rtcb1_2, rtcb2_2, rtcb3_2, bf_rtcb_2, rtcb1_3, rtcb2_3, rtcb3_3, bf_rtcb_3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",legend=attr(x=0.75,y=1)))

savefig(rtcbp, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/rtcb_atp.svg")


plot([rh1_1, rh2_1, rh3_1, bf_rh_1, rh1_2, rh2_2, rh3_2, bf_rh_2, rh1_3, rh2_3, rh3_3, bf_rh_3])


lam_range = range(0.014,stop=0.04,length=3)

lam_range = [0.01, 0.014, 0.02]
bfs=[]; dfs=[];
for i in (lam_range)
    copyparams = deepcopy(params_for_ssval_setup)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

rtcb1_2, rtcb2_2, rtcb3_2, rtcr1_2, rtcr2_2, rtcr3_2, rh1_2, rh2_2, rh3_2, rt1_2, rt2_2, rt3_2, rd1_2, rd2_2, rd3_2, rtca1_2, rtca2_2, rtca3_2, rma1_2, rma2_2, rma3_2, rmr1_2, rmr2_2, rmr3_2, rmb1_2, rmb2_2, rmb3_2, bf_rtcb_2, bf_rtca_2, bf_rtcr_2, bf_rma_2, bf_rmb_2, bf_rmr_2, bf_rh_2, bf_rd_2, bf_rt_2 = plot_all_bf(bfs[2],dfs[2], :mediumorchid, lam_range[2], "λ")
rtcb1_3, rtcb2_3, rtcb3_3, rtcr1_3, rtcr2_3, rtcr3_3, rh1_3, rh2_3, rh3_3, rt1_3, rt2_3, rt3_3, rd1_3, rd2_3, rd3_3, rtca1_3, rtca2_3, rtca3_3, rma1_3, rma2_3, rma3_3, rmr1_3, rmr2_3, rmr3_3, rmb1_3, rmb2_3, rmb3_3, bf_rtcb_3, bf_rtca_3, bf_rtcr_3, bf_rma_3, bf_rmb_3, bf_rmr_3, bf_rh_3, bf_rd_3, bf_rt_3 = plot_all_bf(bfs[3],dfs[3], :purple, lam_range[3], "λ")

a = (scatter(x=dfs[1].kdam, y=dfs[1].rtcb, line=attr(width=3, color=:plum), name="λ = $(lam_range[1])"))

rtcb_lam = plot([a, rtcb1_2, rtcb2_2, rtcb3_2, bf_rtcb_2, rtcb1_3, rtcb2_3, rtcb3_3, bf_rtcb_3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",legend=attr(x=0.75,y=1)))

savefig(rtcb_lam, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/rtcb_lam.svg")
