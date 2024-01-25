using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using Plots, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");


atp_range = range(400,stop=5500,length=50)
lam_range = range(0.001,stop=0.04,length=50)

wab_range1 = 10 .^range(log10(0.056/5e4),log10(0.056/1e1),length=5)
results_wab = get_bs_region_results(atp_range, :atp, lam_range, :lam, wab_range1, :ω_ab, params_bf)
p_atplam_wab = plot_bs_region_same_plot(atp_range, lam_range, results_wab, "changing ω_ab when ω_r = $(round.(ω_r; sigdigits=2))", wab_range1, "ATP", "λ")

Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wab.svg")

wr_range1 = 10 .^ range(log10(0.01/5e7),log10(0.01/1e4),length=5)
results_wr = get_bs_region_results(atp_range, :atp, lam_range, :lam, wr_range1, :ω_r, params_bf)
p_atlam_wr = plot_bs_region_same_plot(atp_range, lam_range, results_wr, "changing ω_r when ω_ab = $(round.(ω_r; sigdigits=2))", wr_range1, "ATP", "λ")

Plots.savefig(p_atlam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wr.svg")

println(wab_range1)

ω_ab
ω_r
0.01/3e4










function numerical_param(param_range1, param1, param_range2, param2, params_trna1, kdam_range)
    # params_trna1 = deepcopy(params_trna)
    res1=[]
    param_combos=[]
    for i in ProgressBar(param_range1)
        # res1=[]
        params_trna1[param1] = i
        for j in param_range2
            params_trna1[param2] = j
            res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range)
            push!(res1, res_rtcb1)
            push!(param_combos, (i,j))
        end   
        # push!(res, res1) 
    end
    return res1, param_combos
end



# plot([scatter(x=kdam_range, y=res[i], name="$(param_combos[i])") for i in range(1,length(res))])

function bistable_q(res, param_combos)
    bistables=[]; 
    for j in range(1, length(res))
        for i in range(1,length(res[j])-1)
            if res[j][i+1] < 0.0001 * res[j][i]
                # push!(bistables, ((res[j][i+1], res[j][i])))
                push!(bistables, param_combos[j])
                break 
            # else
            #     print("monostable")
            end
        end
    end

    return bistables
end


# plot([scatter(x=[bistables[i][1]],y=[bistables[i][2]], mode="markers") for i in range(1, length(bistables))], Layout(xaxis_range=(400,5500), yaxis_range=(0.001,0.04)))


function get_max_min_yvals(param_range1, bistables)
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        vals=[]
        for j in range(1, length(bistables))
            if i == bistables[j][1]
                push!(vals, bistables[j][2])
            end
        end
        push!(max_, maximum(vals))
        push!(min_, minimum(vals))
    end
    return max_,min_
end


function numerical_triple_param_plotting(param_range1, param1, param_range2, param2, params_trna, kdam_range, w_range, w_change)
    p = Plots.plot()
    params_trna1 = deepcopy(params_trna)
    for i in w_range
        params_trna1[w_change] = i
        # @show params_trna1[w_change], params_trna1[param1], params_trna1[param2]
	    res, param_combos = numerical_param(param_range1, param1, param_range2, param2, params_trna1, kdam_range)
        bistables = bistable_q(res, param_combos)
        max_,min_ = get_max_min_yvals(param_range1, bistables)

        p = Plots.plot!(param_range1, max_; fillrange=(min_), fillalpha = 0.45, 
        label="$(@sprintf "%g" i)", xlims=(minimum(param_range1), maximum(param_range1)), 
        ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2", tickfontsize=15,
        guidefontsize=18,guidefont="sans-serif",)
    end

    
    return p 
end

using Plots 
kdam_range = range(0,400,length=50)
atp_range = range(400, stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)

wab_range1 = (10 .^ range(-4, stop=-1, length = 5)) #0.05 normally

p_atplam_wab = numerical_triple_param_plotting(atp_range, :atp, lam_range, :lam, params_trna, kdam_range, wab_range1, :ω_ab)
Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wab.svg")

wr_range1 = (10 .^ range(-5, stop=-2, length = 5))
p_atplam_wr = numerical_triple_param_plotting(atp_range, :atp, lam_range, :lam, params_trna, kdam_range, wr_range1, :ω_r)
Plots.savefig(p_atplam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wr.svg")
