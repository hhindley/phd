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

mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/"
dir = "dec_kdam_1908" # change this! 
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

kdam_range1 = range(0.01,1.5,length=15)
kdam_range2 = reverse(kdam_range1)[2:end]

df = DataFrame(kdam=kdam_range2, time=zeros(length(kdam_range2)))
for i in eachindex(kdam_range2)
    println("starting $i, $(Dates.now())")
    time_taken = @elapsed run_stoch(X0, 150, kdam_range2[i], joinpath(folderpath,"kdam_$(kdam_range2[i]).dat"))    
    df.time[i] = time_taken 
    # init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    df_ss = CSV.read(joinpath(folderpath,"kdam_$(kdam_range2[i]).dat"), DataFrame; header=false)
    ss_region_start = Int(round(length(df_ss[!,1])*0.9))
    # println("ss region start: ", ss_region_start)
    # println("length: ", length(df_ss[!,1]))
    ss_region = df_ss[ss_region_start:end,:]
    # println("ss region length: ", length(ss_region[!,1]))
    init1 = [mean(ss_region[!,col]) for col in names(ss_region[:,3:end-2])]
    global X0 = collect(get_X0(indV, init1)')
    println("calculated new X0 and finished $i, $(Dates.now())")
end

CSV.write(joinpath(mainpath, "$time_file"), df)

println("total time = $(sum(df.time)/60/60) hours")

println("starting file conversion for $dir")

arrow_conv(joinpath(mainpath, dir), joinpath(mainpath, final_path))

print("finished file conversion for $dir")

println("making histograms for $dir")

create_histogram_files(mainpath, final_path)

print("finished making histograms for $dir")
