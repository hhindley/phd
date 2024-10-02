using Dates, StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/histograms/make_hists.jl"))

high_kdam = false

n= 10000 # number of cell cycles
options = Dict(
"threshold"  =>  150.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 1  # for sampling every x mins
)

X0_init = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')

println("starting X0 calc")
X0 = run_stoch(X0_init, options["threshold"], 0, "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat")
X0[vidx(:V)] = 1
# getssX0 = true
# if getssX0
#     fout=open("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat","w")
#     propen, S, propList = defineStochModel(par, indV)
#     nx = indV.nrOfItems-1
#     prop(X) = propen(X[1:nx])
#     X0 = hybrid_algo(X0, options, prop, S, out=fout)
#     X0[vidx(:V)] = 1 #1e-15
#     CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
# else
#     X0 = CSV.read("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
# end
println("finished X0 calc")

if high_kdam
    X0_high = run_stoch(X0, options["threshold"], 1.5, "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat")
    X0_high[vidx(:V)] = 1
    X0 = X0_high
    mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/keyvals2_high_kdam"
else
    mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/keyvals2_low_kdam"
end


date = Dates.format(Dates.now(), "ddmm")
dir_num = 9 # change this! 
dir = "$(date)_$dir_num" 

folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

kdam_vals = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

df = DataFrame(kdam=kdam_vals, time=zeros(length(kdam_vals)))
for i in eachindex(kdam_vals)
    println("starting $i, $(Dates.now())")
    time_taken = @elapsed run_stoch(X0, options["threshold"], kdam_vals[i], joinpath(folderpath,"kdam_$(kdam_vals[i]).dat"))
    df.time[i] = time_taken
    println("finished $i, $(Dates.now())")
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion for $dir")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

println("finished file conversion for $dir")

println("making histograms for $dir")

create_histogram_files(mainpath, final_path)

print("finished making histograms for $dir")




dir = "$(date)_$(dir_num+1)" 

folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

df = DataFrame(kdam=kdam_vals, time=zeros(length(kdam_vals)))
for i in eachindex(kdam_vals)
    println("starting $i, $(Dates.now())")
    time_taken = @elapsed run_stoch(X0, 150, kdam_vals[i], joinpath(folderpath,"kdam_$(kdam_vals[i]).dat"))
    df.time[i] = time_taken
    println("finished $i, $(Dates.now())")
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion for $dir")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

println("finished file conversion for $dir")

println("making histograms for $dir")

create_histogram_files(mainpath, final_path)

print("finished making histograms for $dir")
