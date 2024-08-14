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

println("starting X0 calc")
getssX0 = true
if getssX0
    fout=open("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    X0[vidx(:V)] = 1
    CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end

println("finished X0 calc")


mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing"
dir = "kdam_and_thresh_test_1408" # change this! 
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"


kdam_vals = [0.075, 0.1, 0.1, 0.3, 0.5]

# thresh_scaling = [14000, 13500, 13000, 12500, 12000, 11500, 11000, 3800, 2000, 440]

# thresh_vals = kdam_vals.*thresh_scaling

thresh_vals = [150, 150, 160, 160, 170]

df = DataFrame(kdam=kdam_vals, thresh=thresh_vals, time=zeros(length(kdam_vals)))
for i in eachindex(kdam_vals)
    println("starting $i, $(Dates.now())")
    time_taken = @elapsed run_stoch(X0, thresh_vals[i], kdam_vals[i], joinpath(folderpath,"kdam_$(kdam_vals[i])_thresh_$(thresh_vals[i]).dat"))
    df.time[i] = time_taken
    println("finished $i, $(Dates.now())")
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion for $dir")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

print("finished file conversion for $dir")

println("making histograms for $dir")

create_histogram_files(mainpath, final_path)

print("finished making histograms for $dir")



