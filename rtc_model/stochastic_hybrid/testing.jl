using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

PATH = "/home/hollie_hindley/Documents"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")
include("$PATH/stochastic_hybrid/indexing.jl")

options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=>  [],       # Reactions to be treated determinisitically
    "tspan"     =>   1e6,     # Max time for cell cycle
    "samplingFreq"  => 10     # for sampling every x mins
)


X0 = collect(get_X0(indV)')
par = collect(get_par(indP)')


getssX0 = false

if getssX0
    fout=open("$PATH/stochastic_hybrid/all_X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    CSV.write("$PATH/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("$PATH/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end

 
function prop(X)
    nx = indV.nrOfItems - 1
    propen(X[1:nx]) 
end
function run_stoch(X0, thresh, kdam)
    par[pidx(:kdam)] = kdam
    threshold = thresh # set to zero for a deterministic result
    options["threshold"] = threshold
    fout=open("$PATH/stochastic_hybrid/threshold_testing/threshold_$thresh.dat","w")
    global propen, S, propList  = defineStochModel(par, indV)    
    solu = hybrid_algo(X0, options, prop, S, out=fout)
    close(fout)
end


time_taken = @belapsed run_stoch($X0, $i, 0.1)
