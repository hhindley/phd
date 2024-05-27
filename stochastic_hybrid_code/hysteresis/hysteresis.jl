using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

PATH = "/home/hollie_hindley/Documents"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")

# include("/home/hollie_hindley/Documents/stochastic_hybrid/run_rtc_orig.jl")


n= 1000#200 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 0.1  # for sampling every x mins
)

X0 = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')

kdam_range1 = range(0,1.5,length=50)
kdam_range2 = reverse(kdam_range1)[2:end]

# for i in ProgressBar(kdam_range1)
#     time_taken = @elapsed run_stoch(X0, 50, i, "hysteresis/inc_kdam/kdam_$i.dat")
#     df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/inc_kdam/kdam_$i.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
#     init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
#     global X0 = collect(get_X0(indV, init1)')
# end

for i in ProgressBar(kdam_range2)
    time_taken = @elapsed run_stoch(X0, 50, i, "hysteresis/dec_kdam/kdam_$i.dat")
    df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/dec_kdam/kdam_$i.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
    init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    global X0 = collect(get_X0(indV, init1)')
end