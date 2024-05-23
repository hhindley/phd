using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, CategoricalArrays

PATH = "/home/hollie_hindley/Documents"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")
include("$PATH/stochastic_hybrid/indexing.jl")


dfs = [DataFrame(CSV.File(joinpath("$PATH/stochastic_hybrid/threshold_testing", file), header=["time", "rm_a", "rh", "TotProp"])) for file in readdir("$PATH/stochastic_hybrid/threshold_testing")]
df_time = DataFrame(CSV.File("$PATH/stochastic_hybrid/times.csv", header=["threshold", "time"]))[2:end,:]

df_time.time = [parse.(Float64, subarray) for subarray in df_time.time];
df_time.threshold = [parse.(Float64, subarray) for subarray in df_time.threshold];

plot(scatter(x=df_time.threshold, y=df_time.time/60), Layout(xaxis_title="threshold", yaxis_title="time (minutes)"))
sum(df_time.time)/60/60

dfs[end].bins = ceil.(Int, (1:nrow(dfs[end]))/20)
df_grouped = combine(first, groupby(dfs[end], :bins))

traces_rma=[]
traces_rh=[]
for i in range(1,length(dfs))
    dfs[i].bins = ceil.(Int, (1:nrow(dfs[i]))/20)
    df_grouped = combine(first, groupby(dfs[i], :bins))
    push!(traces_rma,plot(histogram(x=df_grouped.rm_a, name="$(round.(df_time.threshold[i], digits=3))", nbinsx=50)))
    push!(traces_rh,plot(histogram(x=df_grouped.rh, name="$(round.(df_time.threshold[i], digits=3))", nbinsx=50)))
end

p = [traces_rma[1] traces_rma[2] traces_rma[3] traces_rma[4] traces_rma[5];
 traces_rma[6] traces_rma[7] traces_rma[8] traces_rma[9] traces_rma[10];
    traces_rma[11] traces_rma[12] traces_rma[13] traces_rma[14] traces_rma[15];
    traces_rma[16] traces_rma[17] traces_rma[18] traces_rma[19] traces_rma[20]];

relayout!(p, title_text="RtcA mRNA")
p

p1 = [traces_rh[1] traces_rh[2] traces_rh[3] traces_rh[4] traces_rh[5];
traces_rh[6] traces_rh[7] traces_rh[8] traces_rh[9] traces_rh[10];
    traces_rh[11] traces_rh[12] traces_rh[13] traces_rh[14] traces_rh[15];
    traces_rh[16] traces_rh[17] traces_rh[18] traces_rh[19] traces_rh[20]];
relayout!(p1, title_text="Rh")   
p1

plot([scatter(x=dfs[i].time, y=dfs[i].rm_a, name="$(round.(df_time.threshold[i], digits=3))") for i in range(1,length(dfs))], Layout(xaxis_title="time", yaxis_title="rm_a")) 
plot([scatter(x=dfs[i].time, y=dfs[i].rh, name="$(round.(df_time.threshold[i], digits=3))") for i in range(1,length(dfs))], Layout(xaxis_title="time", yaxis_title="rh")) 


means_rma = [mean(dfs[i].rm_a) for i in range(1,length(dfs))]
means_rh = [mean(dfs[i].rh) for i in range(1,length(dfs))]

cv_rma = [std(dfs[i].rm_a) / mean(dfs[i].rm_a) for i in range(1,length(dfs))]
cv_rh = [std(dfs[i].rh) / mean(dfs[i].rh) for i in range(1,length(dfs))]

plot(scatter(x=df_time.threshold, y=means_rma), Layout(xaxis_title="threshold", yaxis_title="mean rm_a"))
plot(scatter(x=df_time.threshold, y=means_rh), Layout(xaxis_title="threshold", yaxis_title="mean rh"))


plot(scatter(x=df_time.threshold, y=cv_rma), Layout(xaxis_title="threshold", yaxis_title="CV rm_a"))
plot(scatter(x=df_time.threshold, y=cv_rh), Layout(xaxis_title="threshold", yaxis_title="CV rh"))



p = plot(scatter(x=dfs[1].time, y=dfs[1].rm_a))#, name="$(round.(df_time.threshold[i], digits=3))") for i in range(1,length(dfs))], Layout(xaxis_title="time", yaxis_title="rm_a")) 
p1 = plot(scatter(x=dfs[1].time, y=dfs[1].rh))# name="$(round.(df_time.threshold[i], digits=3))") for i in range(1,length(dfs))], Layout(xaxis_title="time", yaxis_title="rh")) 

[p; p1]