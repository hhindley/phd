using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

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

println("starting X0 calc")
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

println("finished X0 calc")


# threshold_vals = range(10,310,length=20)
# threshold_vals_new = collect(range(threshold_vals[10], threshold_vals[14], length=5))
# pushfirst!(threshold_vals_new, 160)
# push!(threshold_vals_new, 210)

# threshold_vals = range(100,250,length=10) # 16/07 and 18/07 threshold values
threshold_vals = range(10, 150, length=15) # 18/07 smaller range threshold values 

mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/threshold_testing"
dir = "2308_vs_det_0dam" # change this! 
println(dir)
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"


df = DataFrame(threshold=threshold_vals, time=zeros(length(threshold_vals)))
# Threads.@threads :static for i in eachindex(threshold_vals_new)
for i in eachindex(threshold_vals)
    println("starting $i")
    time_taken = @elapsed run_stoch(X0, threshold_vals[i], 0, joinpath(folderpath,"thresh_$(threshold_vals[i]).dat"))
    df.time[i] = time_taken
    println("finished $i")
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

println("finished file conversion!")

println("making histograms")

create_histogram_files(mainpath, final_path)

println("finished making histograms!")



