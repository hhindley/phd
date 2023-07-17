using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")

kin_range = range(0,stop=0.1,length=10)
atp_range = range(300,stop=5500,length=10)
lam_range = range(0.001,stop=0.04,length=10)
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



p1_rma = plot_different_bp_lines(atp_range, params1, :atp, :rm_a, 3.);
p2_rma = plot_different_bp_lines(kin_range, params1, :kin, :rm_a, 3.);
p3_rma = plot_different_bp_lines(lam_range, params1, :lam, :rm_a, 3.);
p4_rma = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_a, 3.);
p5_rma = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_a, 3.);

l = @layout [a b; c d; e];
all_rma = plot(p1_rma, p2_rma, p3_rma, p4_rma, p5_rma, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rma, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rma.svg")

p1_rmr = plot_different_bp_lines(atp_range, params1, :atp, :rm_r, 3.);
p2_rmr = plot_different_bp_lines(kin_range, params1, :kin, :rm_r, 3.);
p3_rmr = plot_different_bp_lines(lam_range, params1, :lam, :rm_r, 3.);
p4_rmr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_r, 3.);
p5_rmr = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_r, 3.);

all_rmr = plot(p1_rmr, p2_rmr, p3_rmr, p4_rmr, p5_rmr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rmr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rmr.svg")

p1_rtca = plot_different_bp_lines(atp_range, params1, :atp, :rtca, 3.);
p2_rtca = plot_different_bp_lines(kin_range, params1, :kin, :rtca, 3.);
p3_rtca = plot_different_bp_lines(lam_range, params1, :lam, :rtca, 3.);
p4_rtca = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtca, 3.);
p5_rtca = plot_different_bp_lines(wr_range, params1, :ω_r, :rtca, 3.);

all_rtca = plot(p1_rtca, p2_rtca, p3_rtca, p4_rtca, p5_rtca, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rtca, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtca.svg")

p1_rtcr = plot_different_bp_lines(atp_range, params1, :atp, :rtcr, 3.);
p2_rtcr = plot_different_bp_lines(kin_range, params1, :kin, :rtcr, 3.);
p3_rtcr = plot_different_bp_lines(lam_range, params1, :lam, :rtcr, 3.);
p4_rtcr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtcr, 3.);
p5_rtcr = plot_different_bp_lines(wr_range, params1, :ω_r, :rtcr, 3.);

all_rtcr = plot(p1_rtcr, p2_rtcr, p3_rtcr, p4_rtcr, p5_rtcr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rtcr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtcr.svg")

p1_rh = plot_different_bp_lines(atp_range, params1, :atp, :rh, 3.);
p2_rh = plot_different_bp_lines(kin_range, params1, :kin, :rh, 3.);
p3_rh = plot_different_bp_lines(lam_range, params1, :lam, :rh, 3.);
p4_rh = plot_different_bp_lines(wab_range, params1, :ω_ab, :rh, 3.);
p5_rh = plot_different_bp_lines(wr_range, params1, :ω_r, :rh, 3.);

all_rh = plot(p1_rh, p2_rh, p3_rh, p4_rh, p5_rh, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rh.svg")

p1_rt = plot_different_bp_lines(atp_range, params1, :atp, :rt, 3.);
p2_rt = plot_different_bp_lines(kin_range, params1, :kin, :rt, 3.);
p3_rt = plot_different_bp_lines(lam_range, params1, :lam, :rt, 3.);
p4_rt = plot_different_bp_lines(wab_range, params1, :ω_ab, :rt, 3.);
p5_rt = plot_different_bp_lines(wr_range, params1, :ω_r, :rt, 3.);

all_rt = plot(p1_rt, p2_rt, p3_rt, p4_rt, p5_rt, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rt, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rt.svg")

p1_rd = plot_different_bp_lines(atp_range, params1, :atp, :rd, 3.);
p2_rd = plot_different_bp_lines(kin_range, params1, :kin, :rd, 3.);
p3_rd = plot_different_bp_lines(lam_range, params1, :lam, :rd, 3.);
p4_rd = plot_different_bp_lines(wab_range, params1, :ω_ab, :rd, 3.);
p5_rd = plot_different_bp_lines(wr_range, params1, :ω_r, :rd, 3.);

all_rd = plot(p1_rd, p2_rd, p3_rd, p4_rd, p5_rd, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
# savefig(all_rd, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rd.svg")


l = @layout [a b; c d; e f h]
atp_plots = plot(p1_rma, p1_rtca, p1_rmr, p1_rtcr, p1_rh, p1_rt, p1_rd, layout=l, size=(1200,800), plot_title="ATP = $([round.(i; digits=4) for i in atp_range])", titlefontsize=6, left_margin=4Plots.mm)
kin_plots = plot(p2_rma, p2_rtca, p2_rmr, p2_rtcr, p2_rh, p2_rt, p2_rd, layout=l, size=(1200,800), plot_title="kin = $([round.(i; digits=4) for i in kin_range])", titlefontsize=6, left_margin=4Plots.mm)
lam_plots = plot(p3_rma, p3_rtca, p3_rmr, p3_rtcr, p3_rh, p3_rt, p3_rd, layout=l, size=(1200,800), plot_title="λ = $([round.(i; digits=4) for i in lam_range])", titlefontsize=6, left_margin=4Plots.mm)
wab_plots = plot(p4_rma, p4_rtca, p4_rmr, p4_rtcr, p4_rh, p4_rt, p4_rd, layout=l, size=(1200,800), plot_title="ω_ab = $([round.(i; digits=4) for i in wab_range])", titlefontsize=6, left_margin=4Plots.mm)
wr_plots = plot(p5_rma, p5_rtca, p5_rmr, p5_rtcr, p5_rh, p5_rt, p5_rd, layout=l, size=(1200,800), plot_title="ω_r = $([round.(i; digits=4) for i in wr_range])", titlefontsize=6, left_margin=4Plots.mm)

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
