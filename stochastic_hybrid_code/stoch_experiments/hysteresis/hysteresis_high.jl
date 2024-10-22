using Dates, StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/histograms/make_hists.jl"))


n= 10000 # number of cell cycles
options = Dict(
"threshold"  =>  150.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 1 # for sampling every x mins
)

X0_init = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')

println("starting X0 calc")
X0 = run_stoch(X0_init, options["threshold"], 0, "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat")
X0[vidx(:V)] = 1
println("finished X0 calc, X0: $X0")

kdam_init_val_high = 1.5
kdam_middle = 0.8

date = Dates.format(Dates.now(), "ddmmyy")

mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/new/high_kdam"
dir = "$date" 
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

n = 10

df = DataFrame(kdam=fill(kdam_middle, n), time=zeros(n))

for i in 1:n
    println("starting X0 calc")
    X0_high = run_stoch(X0, options["threshold"], kdam_init_val_high, "/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat")
    X0_high[vidx(:V)] = 1
    println("X0_high: $X0_high")
    println("starting run $i at $(Dates.now())")
    time_taken = @elapsed run_stoch(X0_high, options["threshold"], kdam_middle, joinpath(folderpath,"run_$i.dat"))    
    println("finished run $i at $(Dates.now())")
    df.time[i] = time_taken
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion for $dir")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

print("finished file conversion for $dir")

println("making histograms for $dir")

create_histogram_files(mainpath, final_path)

print("finished making histograms for $dir")

move_time_file(mainpath)