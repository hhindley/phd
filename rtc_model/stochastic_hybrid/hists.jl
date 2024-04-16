using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, CategoricalArrays

PATH = "/home/holliehindley/phd"

include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/rtc_model/parameters/rtc_params_molecs.jl")
include("$PATH/rtc_model/stochastic_hybrid/indexing.jl")
include("$PATH/rtc_model/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/rtc_model/stochastic_hybrid/stoch_model.jl")
include("$PATH/rtc_model/stochastic_hybrid/indexing.jl")



dfs = [DataFrame(CSV.File(joinpath("$PATH/rtc_model/stochastic_hybrid/threshold_testing2", file), header=["time", "rm_a", "rh"])) for file in readdir("$PATH/rtc_model/stochastic_hybrid/threshold_testing2")]
df_time = DataFrame(CSV.File("$PATH/rtc_model/stochastic_hybrid/times.csv", header=["threshold", "time"]))[2:end,:]

df_time.time = [parse.(Float64, subarray) for subarray in df_time.time];
df_time.threshold = [parse.(Float64, subarray) for subarray in df_time.threshold];

plot(scatter(x=df_time.threshold, y=df_time.time/60), Layout(xaxis_title="threshold", yaxis_title="time (minutes)"))
sum(df_time.time)/60/60

dfs[end].bins = ceil.(Int, (1:nrow(dfs[end]))/20)
df_grouped = combine(first, groupby(dfs[end], :bins))


dfs[end].bins = ceil.(Int, (1:nrow(dfs[end]))/20)
df_grouped = combine(first, groupby(dfs[end], :bins))

traces_rma=[]
traces_rh=[]
for i in range(1,length(dfs))
    dfs[i].bins = ceil.(Int, (1:nrow(dfs[i]))/20)
    df_grouped = combine(first, groupby(dfs[i], :bins))
    push!(traces_rma,(histogram(x=df_grouped.rm_a, name="$(round.(df_time.threshold[i], digits=3))", nbinsx=50)))
    push!(traces_rh,(histogram(x=df_grouped.rh, name="$(round.(df_time.threshold[i], digits=3))", nbinsx=50)))
end

p = make_subplots(rows=10, cols=5, shared_xaxes=false)
positions = [(row, col) for row in 1:10 for col in 1:5]
for (trace, (row, col)) in zip(traces_rma, positions)
    add_trace!(p, trace, row=row, col=col)
end
relayout!(p, title_text="RtcA mRNA")
p

p1 = make_subplots(rows=10, cols=5, shared_xaxes=true)
for (trace, (row, col)) in zip(traces_rh, positions)
    add_trace!(p1, trace, row=row, col=col)
end
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


