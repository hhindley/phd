using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics

include("/home/hollie_hindley/Documents/stochastic_hybrid/analysis_funcs.jl")

kdam_vals = range(0,2,length=21) # vals for kdam_test

kdam_vals = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4] # vals for kdam_test1

df_times = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/times.csv")) # times for kdam_test1

plot(scatter(x=df_times.kdam, y=df_times.time/60))

# load all dataframes
dfs = [DataFrame(CSV.File(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_test1", file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])) for file in readdir("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_test1")]


# sort dataframes
props, df_rs, df_ps = df_sort_all(dfs)


# stochastic reactions
df_reacts, tot_stoch = all_react_dfs(dfs, df_rs)

p_stochreacts = plot([bar(x=df_reacts[i].reaction, y=df_reacts[i].count, name="$(kdam_vals[i])") for i in eachindex(kdam_vals)], Layout(barmode="group", yaxis_type="log"))
p_totstoch = plot(bar(x=kdam_vals, y=tot_stoch))

# histograms 
groups, bins = all_hists(df_ps, kdam_vals, 1000)
for specie in species_rtc
    display(plot([histogram(x=groups[i][:,specie], nbins=bins[i], opacity=0.6, name="$(kdam_vals[i])") for i in eachindex(kdam_vals)], Layout(barmode="overlay", title="$specie", yaxis_type="log")))
end

# plotting individual species
x = 8
plot(scatter(x=df_ps[x].time, y=df_ps[x].rm_a))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rm_b))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rm_r))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rtca))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rtcb))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rtcr))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rh))#./df_ps[x].volume))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rd))#./df_ps[x].volume))
plot(scatter(x=df_ps[x].time, y=df_ps[x].rt))#./df_ps[x].volume))

plot([scatter(x=df_ps[x].time, y=df_ps[x][:,i]./df_ps[x].volume, name="$i") for i in [:rh, :rd, :rt]])
plot([scatter(x=df_ps[x].time, y=df_ps[x][:,i]./df_ps[x].volume, name="$i") for i in [:rm_a, :rtca]])
plot([scatter(x=df_ps[x].time, y=df_ps[x][:,i]./df_ps[x].volume, name="$i") for i in [:rm_b, :rtcb]])
plot([scatter(x=df_ps[x].time, y=df_ps[x][:,i]./df_ps[x].volume, name="$i") for i in [:rm_r, :rtcr]])

# % of expression
exp_df = all_exp(df_ps, kdam_vals)

p_exp = plot([scatter(x=exp_df.kdam, y=exp_df[:,i], name="$(names(exp_df)[i])") for i in 2:7], Layout(xaxis_title="kdam", yaxis_title="% of expression (when species > 1)"))

# mean values
df_av = mean_vals(df_ps, kdam_vals, 5000)

p_av = plot([scatter(x=df_av.kdam, y=df_av[:,i], name="$(names(df_av)[i])") for i in 2:10], Layout(xaxis_title="kdam", yaxis_title="mean molecule number"))

# expression + stoch reactions 
plot([scatter(x=df_ps[x].time, y=df_ps[x].rm_r, name="rm_r"), 
scatter(x=filter(row -> row[:event] == 2, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 2, df_rs[x]).rm_r)*1.1], length(df_rs[x].time)), name="$(react_names[2])", mode="markers"),
scatter(x=filter(row -> row[:event] == 13, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 13, df_rs[x]).rm_r)*10.1], length(df_rs[x].time)), name="$(react_names[13])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rm_a, name="rm_a"), 
scatter(x=filter(row -> row[:event] == 1, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 1, df_rs[x]).rm_a)*1.1], length(df_rs[x].time)), name="$(react_names[1])", mode="markers"),
scatter(x=filter(row -> row[:event] == 11, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 11, df_rs[x]).rm_a)*1.2], length(df_rs[x].time)), name="$(react_names[11])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rm_b, name="rm_b"), 
scatter(x=filter(row -> row[:event] == 1, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 1, df_rs[x]).rm_b)*1.1], length(df_rs[x].time)), name="$(react_names[1])", mode="markers"),
scatter(x=filter(row -> row[:event] == 12, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 12, df_rs[x]).rm_b)*1.25], length(df_rs[x].time)), name="$(react_names[12])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rtcr, name="rtcr"), 
scatter(x=filter(row -> row[:event] == 5, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 5, df_rs[x]).rtcr)*1.1], length(df_rs[x].time)), name="$(react_names[5])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rtca, name="rtca"), 
scatter(x=filter(row -> row[:event] == 3, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 3, df_rs[x]).rtca)*50.1], length(df_rs[x].time)), name="$(react_names[3])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rtcb, name="rtcb"), 
scatter(x=filter(row -> row[:event] == 4, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 4, df_rs[x]).rtcb)*38.1], length(df_rs[x].time)), name="$(react_names[4])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rh, name="rh"), 
# scatter(x=filter(row -> row[:event] == 6, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 6, df_rs[x]).rh)*1.1], length(df_rs[x].time)), name="$(react_names[6])", mode="markers"),
# scatter(x=filter(row -> row[:event] == 7, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 7, df_rs[x]).rh)*1.1], length(df_rs[x].time)), name="$(react_names[7])", mode="markers"),
scatter(x=filter(row -> row[:event] == 8, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 8, df_rs[x]).rh)*1.2], length(df_rs[x].time)), name="$(react_names[8])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rt, name="rt"), 
# scatter(x=filter(row -> row[:event] == 6, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 6, df_rs[x]).rh)*1.1], length(df_rs[x].time)), name="$(react_names[6])", mode="markers"),
scatter(x=filter(row -> row[:event] == 9, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 9, df_rs[x]).rt)*7.7], length(df_rs[x].time)), name="$(react_names[9])", mode="markers"),
scatter(x=filter(row -> row[:event] == 8, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 8, df_rs[x]).rt)*180.2], length(df_rs[x].time)), name="$(react_names[8])", mode="markers")])

plot([scatter(x=df_ps[x].time, y=df_ps[x].rd, name="rd"), 
# scatter(x=filter(row -> row[:event] == 6, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 6, df_rs[x]).rh)*1.1], length(df_rs[x].time)), name="$(react_names[6])", mode="markers"),
scatter(x=filter(row -> row[:event] == 10, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 10, df_rs[x]).rd)*1.7], length(df_rs[x].time)), name="$(react_names[10])", mode="markers"),
scatter(x=filter(row -> row[:event] == 8, df_rs[x]).time, y=repeat([maximum(filter(row -> row[:event] == 8, df_rs[x]).rd)*1.2], length(df_rs[x].time)), name="$(react_names[8])", mode="markers")])





df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_test1/kdam_0.05.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

df_props, df_r, df_p = df_sort(df)

plot(scatter(x=df_p.time, y=df_p.rtca))

plot(scatter(y=df_p.rm_a))


p_props = plot([scatter(x=df_p.time[1:1000:end], y=df_props[1:1000:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])

display(p_props)

plot(scatter(x=[1,2,3,4,],y=[2,2,2,2]))


