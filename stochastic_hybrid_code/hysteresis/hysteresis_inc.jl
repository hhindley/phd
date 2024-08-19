using Dates, StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/threshold_analysis/histograms/make_hists.jl"))


n= 10 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 10/60  # for sampling every x mins
)

X0 = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')

mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/"
dir = "inc_kdam_1608" # change this! 
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

kdam_range1 = range(0,1.5,length=16)

df = DataFrame(kdam=kdam_range1, time=zeros(length(kdam_range1)))
for i in eachindex(kdam_range1)
    println("starting $i, $(Dates.now())")
    time_taken = @elapsed run_stoch(X0, 150, kdam_range1[i], joinpath(folderpath,"kdam_$(kdam_range1[i]).dat"))    
    df.time[i] = time_taken 
    # init1 = [mean(df[:,col]./df.volume) for col in names(eachcol(df[:,3:end-2]))]
    df_ss = DataFrame(CSV.File(IOBuffer(joinpath(folderpath,"kdam_$(kdam_range1[i]).dat")); header=false, skipto=2))
    ss_region = Int(round(length(df.time)-length(df.time)*0.1))
    println("ss_region ",df.time[ss_region+1:end])
    println("df ", df)
    init1 = [mean(df[ss_region+1:end,col]) for col in names(df[:,3:end-2])]
    println("init ",init1)
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



# df = DataFrame(a = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], b=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# names(df)
# [mean(df[18+1:end,col]) for col in names(eachcol(df[:,1:2]))]

# df[18+1:end,:a]