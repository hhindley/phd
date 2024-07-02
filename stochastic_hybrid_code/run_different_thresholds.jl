using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

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
    fout=open(joinpath(homedir(), "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat"),"w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    X0[vidx(:V)] = 1
    CSV.write(joinpath(homedir(), "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat"), DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read(joinpath(homedir(), "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat"), Tables.matrix, header=false)
end

Threads.nthreads()

threshold_vals = range(10,310,length=20)
threshold_vals_new = collect(range(threshold_vals[11], threshold_vals[13], length=3))
pushfirst!(threshold_vals_new, 160)
push!(threshold_vals_new, 210)

df = DataFrame(threshold=threshold_vals_new, time=zeros(length(threshold_vals_new)))

Threads.@threads for i in eachindex(threshold_vals_new)
    println("starting $i")
    time_taken = @elapsed run_stoch(X0, threshold_vals_new[i], 0.05, "/home/hollie_hindley/Documents/stochastic_hybrid/thresh_0107_repeat3/thresh_$(threshold_vals_new[i]).dat")
    df.time[i] = time_taken
    println("finished $i")
end

CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_times_0107.csv", df)

println("total time = $(sum(df.time)/60/60) hours")

arrow_conv(joinpath(homedir(), "Documents/stochastic_hybrid/thresh_0107_repeat3"), joinpath(homedir(), "Documents/stochastic_hybrid/thresh_0107_final_files"))

print("finished!")



i = Threads.Atomic{Int}(0)
ids = zeros(4)
old_is = zeros(4)
Threads.@threads for id in 1:4
    old_is[id] = Threads.atomic_add!(i, id)
    ids[id]= id
end
old_is
i[]
ids