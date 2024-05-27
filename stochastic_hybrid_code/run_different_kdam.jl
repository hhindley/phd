using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

PATH = "/home/hollie_hindley/Documents"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")
include("$PATH/stochastic_hybrid/indexing.jl")

n= 10000 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 0.1  # for sampling every x mins
)

X0 = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')


getssX0 = true
if getssX0
    fout=open("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    X0[vidx(:V)] = 1
    # run_stoch(X0, 0, 0, "X0")
    # df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
    # ss_df = df[1000:end,:]
    # ss = [mean(df[:,col]) for col in names(eachcol(df[:,3:end-2]))]
    # X0 = collect(get_X0(indV, ss)')
    CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end



# n=5
# threshold_vals =  (range(0.01^(1/n),80^(1/n),length=50)) .^ n
kdam_vals = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]

times=[]
for i in ProgressBar(kdam_vals)
    time_taken = @elapsed run_stoch(X0, 50, i, "kdam_test1/kdam_$i.dat")
    push!(times, time_taken)
end

df = DataFrame(kdam=kdam_vals, time=times)

CSV.write("$PATH/stochastic_hybrid/times.csv", df)

println("total time = $(sum(times)/60/60) hours")

