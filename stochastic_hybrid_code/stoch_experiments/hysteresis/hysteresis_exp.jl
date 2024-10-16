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
println("finished X0 calc, X0: $X0")

mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/"
date = Dates.format(Dates.now(), "ddmm")

high = false
num = 10

if high 
    dir_num = "high_$num"
    kdam_range1 = [1.5, 0.8]
else 
    dir_num = "low_$num"
    kdam_range1 = [0.01, 0.8]
end

dir = "$(date)_$dir_num" 
folderpath = joinpath(mainpath, dir)
if !isdir(folderpath)
    mkdir(folderpath)
end
time_file = dir * "_times.csv"
final_path = dir * "_final_files"

df = DataFrame(kdam=kdam_range1, time=zeros(length(kdam_range1)))
for i in eachindex(kdam_range1)
    println("starting $i, $(Dates.now())")
    println("X0: $X0")
    time_taken = @elapsed run_stoch(X0, 150, kdam_range1[i], joinpath(folderpath,"kdam_$(kdam_range1[i]).dat"))    
    df.time[i] = time_taken 
    df_ss = CSV.read(joinpath(folderpath,"kdam_$(kdam_range1[i]).dat"), DataFrame; header=false)
    init1 = [df_ss[end,col] for col in names(df_ss[:,3:end-2])]
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


