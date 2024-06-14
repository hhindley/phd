using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))

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
    fout=open(joinpath(homedir(), "Documents/stochastic_hybrid/X0.dat"),"w")
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
    CSV.write(joinpath(homedir(), "Documents/stochastic_hybrid/X0.dat"), DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read(joinpath(homedir(), "Documents/stochastic_hybrid/X0.dat"), Tables.matrix, header=false)
end



# n=5
# threshold_vals =  (range(0.01^(1/n),80^(1/n),length=50)) .^ n
threshold_vals = range(10,310,length=20)

times=[]
for i in ProgressBar(threshold_vals)
    time_taken = @elapsed run_stoch(X0, i, 0.05, "thresh_test/thresh_$i.dat")
    push!(times, time_taken)
end

df = DataFrame(threshold=threshold_vals, time=times)

CSV.write("$PATH/stochastic_hybrid/thresh_times.csv", df)

println("total time = $(sum(times)/60/60) hours")

# threshold = 50 takes roughly 12 hours i think 
# n=5
# threshold_vals =  ((range(0.01^(1/n),20^(1/n),length=50)) .^ n)

# times=[]
# for i in threshold_vals
#     t = ((((i/0.01)*10)/60)/60)
#     push!(times, t)
# end

# # # plot(scatter(x=threshold_vals, y=times, mode="markers"))

# sum(times)

# t = ((((1200/0.01)*10)/60)/60)/24