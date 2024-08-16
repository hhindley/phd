using Dates, StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/threshold_analysis/histograms/make_hists.jl"))


n= 10000 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 10/60  # for sampling every x mins
)

X0 = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')

kdam_range1 = range(0,1.5,length=50)
kdam_range2 = reverse(kdam_range1)[2:end]

for i in ProgressBar(kdam_range1)
    time_taken = @elapsed run_stoch(X0, 150, i, "hysteresis/inc_kdam/kdam_$i.dat")
    df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/inc_kdam/kdam_$i.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
    # init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    global X0 = collect(get_X0(indV, init1)')
end

for i in ProgressBar(kdam_range2)
    time_taken = @elapsed run_stoch(X0, 150, i, "hysteresis/dec_kdam/kdam_$i.dat")
    df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/dec_kdam/kdam_$i.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
    # init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    global X0 = collect(get_X0(indV, init1)')
end

a = DataFrame(a=[1,2,3,4,5,6,7,8,9,10])

ss_region = Int(length(a.a)-length(a.a)*0.1)

a.a[ss_region:end]