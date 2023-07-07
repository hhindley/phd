using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")

kin_range = range(0,stop=0.05,length=10)
atp_range = range(500,stop=5000,length=10)
lam_range = range(0.001,stop=0.04,length=10)
wr_range = (range(0.00001,stop=0.001,length=10))
wab_range = range(0.01, stop=4, length=10)

function plot_different_bp_lines(param_range, params1, param, specie)
    p = plot();
    copyparams = deepcopy(params1)
    for i in param_range
        params = merge(copyparams, (param=>i,))
        br = get_br(params, initial)
        plot!(p, br, vars = (:param, specie), legend=false)#, label="$param = $i", putspecialptlegend=false)
    end
    return p
end



p1_rma = plot_different_bp_lines(atp_range, params1, :atp, :rm_a);
p2_rma = plot_different_bp_lines(kin_range, params1, :kin, :rm_a);
p3_rma = plot_different_bp_lines(lam_range, params1, :lam, :rm_a);
p4_rma = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_a);
p5_rma = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_a);

l = @layout [a b; c d; e];
all_rma = plot(p1_rma, p2_rma, p3_rma, p4_rma, p5_rma, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rma, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rma.svg")

p1_rmr = plot_different_bp_lines(atp_range, params1, :atp, :rm_r);
p2_rmr = plot_different_bp_lines(kin_range, params1, :kin, :rm_r);
p3_rmr = plot_different_bp_lines(lam_range, params1, :lam, :rm_r);
p4_rmr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_r);
p5_rmr = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_r);

all_rmr = plot(p1_rmr, p2_rmr, p3_rmr, p4_rmr, p5_rmr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rmr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rmr.svg")

p1_rtca = plot_different_bp_lines(atp_range, params1, :atp, :rtca);
p2_rtca = plot_different_bp_lines(kin_range, params1, :kin, :rtca);
p3_rtca = plot_different_bp_lines(lam_range, params1, :lam, :rtca);
p4_rtca = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtca);
p5_rtca = plot_different_bp_lines(wr_range, params1, :ω_r, :rtca);

all_rtca = plot(p1_rtca, p2_rtca, p3_rtca, p4_rtca, p5_rtca, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rtca, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtca.svg")

p1_rtcr = plot_different_bp_lines(atp_range, params1, :atp, :rtcr);
p2_rtcr = plot_different_bp_lines(kin_range, params1, :kin, :rtcr);
p3_rtcr = plot_different_bp_lines(lam_range, params1, :lam, :rtcr);
p4_rtcr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtcr);
p5_rtcr = plot_different_bp_lines(wr_range, params1, :ω_r, :rtcr);

all_rtcr = plot(p1_rtcr, p2_rtcr, p3_rtcr, p4_rtcr, p5_rtcr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rtcr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtcr.svg")

p1_rh = plot_different_bp_lines(atp_range, params1, :atp, :rh);
p2_rh = plot_different_bp_lines(kin_range, params1, :kin, :rh);
p3_rh = plot_different_bp_lines(lam_range, params1, :lam, :rh);
p4_rh = plot_different_bp_lines(wab_range, params1, :ω_ab, :rh);
p5_rh = plot_different_bp_lines(wr_range, params1, :ω_r, :rh);

all_rh = plot(p1_rh, p2_rh, p3_rh, p4_rh, p5_rh, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rh.svg")

p1_rt = plot_different_bp_lines(atp_range, params1, :atp, :rt);
p2_rt = plot_different_bp_lines(kin_range, params1, :kin, :rt);
p3_rt = plot_different_bp_lines(lam_range, params1, :lam, :rt);
p4_rt = plot_different_bp_lines(wab_range, params1, :ω_ab, :rt);
p5_rt = plot_different_bp_lines(wr_range, params1, :ω_r, :rt);

all_rt = plot(p1_rt, p2_rt, p3_rt, p4_rt, p5_rt, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rt, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rt.svg")

p1_rd = plot_different_bp_lines(atp_range, params1, :atp, :rd);
p2_rd = plot_different_bp_lines(kin_range, params1, :kin, :rd);
p3_rd = plot_different_bp_lines(lam_range, params1, :lam, :rd);
p4_rd = plot_different_bp_lines(wab_range, params1, :ω_ab, :rd);
p5_rd = plot_different_bp_lines(wr_range, params1, :ω_r, :rd);

all_rd = plot(p1_rd, p2_rd, p3_rd, p4_rd, p5_rd, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rd, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rd.svg")


l = @layout [a b; c d; e f h]
atp_plots = plot(p1_rma, p1_rtca, p1_rmr, p1_rtcr, p1_rh, p1_rt, p1_rd, layout=l, size=(1200,800), plot_title="ATP = $([round.(i; digits=4) for i in atp_range])", titlefontsize=6, left_margin=4Plots.mm)
kin_plots = plot(p2_rma, p2_rtca, p2_rmr, p2_rtcr, p2_rh, p2_rt, p2_rd, layout=l, size=(1200,800), plot_title="kin = $([round.(i; digits=4) for i in kin_range])", titlefontsize=6, left_margin=4Plots.mm)
lam_plots = plot(p3_rma, p3_rtca, p3_rmr, p3_rtcr, p3_rh, p3_rt, p3_rd, layout=l, size=(1200,800), plot_title="λ = $([round.(i; digits=4) for i in lam_range])", titlefontsize=6, left_margin=4Plots.mm)
wab_plots = plot(p4_rma, p4_rtca, p4_rmr, p4_rtcr, p4_rh, p4_rt, p4_rd, layout=l, size=(1200,800), plot_title="ω_ab = $([round.(i; digits=4) for i in wab_range])", titlefontsize=6, left_margin=4Plots.mm)
wr_plots = plot(p5_rma, p5_rtca, p5_rmr, p5_rtcr, p5_rh, p5_rt, p5_rd, layout=l, size=(1200,800), plot_title="ω_r = $([round.(i; digits=4) for i in wr_range])", titlefontsize=6, left_margin=4Plots.mm)

savefig(atp_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/atp_plots.svg")
savefig(kin_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/kin_plots.svg")
savefig(lam_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/lam_plots.svg")
savefig(wab_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wab_plots.svg")
savefig(wr_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wr_plots.svg")
