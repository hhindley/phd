using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics

include("/home/hollie_hindley/Documents/stochastic_hybrid/analysis_funcs.jl")

kdam_range1 = range(0,1.5,length=50)
kdam_range2 = reverse(kdam_range1)[2:end]

dfs1 = [DataFrame(CSV.File(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/inc_kdam", file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])) for file in readdir("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/inc_kdam")]
dfs2 = [DataFrame(CSV.File(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/dec_kdam", file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])) for file in readdir("/home/hollie_hindley/Documents/stochastic_hybrid/hysteresis/dec_kdam")]

# sort dataframes
props1, df_rs1, df_ps1 = df_sort_all(dfs1)
props2, df_rs2, df_ps2 = df_sort_all(dfs2)

# stochastic reactions
df_reacts1, tot_stoch1 = all_react_dfs(dfs1, df_rs1)
df_reacts2, tot_stoch2 = all_react_dfs(dfs2, df_rs2)

p_stochreacts1 = plot([bar(x=df_reacts1[i].reaction, y=df_reacts1[i].count, name="$(kdam_range1[i])") for i in eachindex(kdam_range1)], Layout(barmode="group", yaxis_type="log"))
p_stochreacts2 = plot([bar(x=df_reacts2[i].reaction, y=df_reacts2[i].count, name="$(kdam_range2[i])") for i in eachindex(kdam_range2)], Layout(barmode="group", yaxis_type="log"))

p_totstoch1 = plot(bar(x=kdam_range1, y=tot_stoch1))
p_totstoch2 = plot(bar(x=kdam_range2, y=tot_stoch2))

# histograms 
groups1, bins1 = all_hists(df_ps1, kdam_range1, 1000)
for specie in species_rtc
    display(plot([histogram(x=groups1[i][:,specie], nbins=bins1[i], opacity=0.6, name="$(kdam_range1[i])") for i in eachindex(kdam_range1)], Layout(barmode="overlay", title="$specie", yaxis_type="log")))
end

groups2, bins2 = all_hists(df_ps2, kdam_range2, 1000)
for specie in species_rtc
    display(plot([histogram(x=groups2[i][:,specie], nbins=bins2[i], opacity=0.6, name="$(kdam_range2[i])") for i in eachindex(kdam_range2)], Layout(barmode="overlay", title="$specie", yaxis_type="log")))
end

# % of expression
exp_df1 = all_exp(df_ps1, kdam_range1)
exp_df2 = all_exp(df_ps2, kdam_range2)

p_exp1 = plot([scatter(x=exp_df1.kdam, y=exp_df1[:,i], name="$(names(exp_df1)[i])") for i in 2:7], Layout(xaxis_title="kdam", yaxis_title="% of expression (when species > 1)"))
p_exp2 = plot([scatter(x=reverse(exp_df2.kdam), y=exp_df2[:,i], name="$(names(exp_df2)[i])") for i in 2:7], Layout(xaxis_title="kdam", yaxis_title="% of expression (when species > 1)"))

# mean values
df_av1 = mean_vals(df_ps1, kdam_range1, 5000)
df_av2 = mean_vals(df_ps2, kdam_range2, 5000)

p_av1 = plot([scatter(x=df_av1.kdam, y=df_av1[:,i], name="$(names(df_av1)[i])") for i in 2:10], Layout(xaxis_title="kdam", yaxis_title="mean molecule number"))
p_av2 = plot([scatter(x=reverse(df_av2.kdam), y=df_av2[:,i], name="$(names(df_av2)[i])") for i in 2:10], Layout(xaxis_title="kdam", yaxis_title="mean molecule number"))

# plotting individual species
x = 8
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rm_a))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rm_b))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rm_r))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rtca))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rtcb))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rtcr))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rh))#./df_ps[x].volume))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rd))#./df_ps[x].volume))
plot(scatter(x=df_ps1[x].time, y=df_ps1[x].rt))#./df_ps[x].volume))

plot([scatter(x=df_ps1[x].time, y=df_ps1[x][:,i]./df_ps1[x].volume, name="$i") for i in [:rh, :rd, :rt]])
plot([scatter(x=df_ps1[x].time, y=df_ps1[x][:,i]./df_ps1[x].volume, name="$i") for i in [:rm_a, :rtca]])
plot([scatter(x=df_ps1[x].time, y=df_ps1[x][:,i]./df_ps1[x].volume, name="$i") for i in [:rm_b, :rtcb]])
plot([scatter(x=df_ps1[x].time, y=df_ps1[x][:,i]./df_ps1[x].volume, name="$i") for i in [:rm_r, :rtcr]])




function compare_sims(i, j, specie)
    return plot([scatter(x=df_ps[i].time, y=df_ps[i][:,specie], name="10000 cc"), scatter(x=df_ps1[j].time, y=df_ps1[j][:,specie], name="1000 cc")])
end


kdam_vals[1]
kdam_range1[3]
compare_sims(1, 3, :rtca)

kdam_vals[2]
kdam_range1[4]
compare_sims(2, 4, :rtca)

kdam_vals[3]
kdam_range1[6]
compare_sims(3, 6, :rtca)

kdam_vals[4]
kdam_range1[8]
compare_sims(4, 8, :rtca)

kdam_vals[5]
kdam_range1[9]
compare_sims(5, 9, :rtca)

kdam_vals[6]    
kdam_range1[11]
compare_sims(6, 11, :rtca)

kdam_vals[7]
kdam_range1[13]
compare_sims(7, 13, :rtca)

kdam_vals[8]
kdam_range1[14]
compare_sims(8, 14, :rtca)

plot(scatter(x=df_ps1[1].time, y=df_ps1[1].rtca))
using Plots 
Plots.plot(df_ps1[1].time, df_ps1[1].rtca)
Plots.plot(df_ps[1].time, df_ps[1].rtca)

using InteractiveViz, GLMakie
ilines(sin, 0, 100)
ilines(df_ps[1].time, df_ps[1].rtca)


p_av = plot([scatter(x=df_av.kdam, y=df_av[:,i], name="$(names(df_av)[i])") for i in 2:10], Layout(xaxis_title="kdam", yaxis_title="mean molecule number", title="10,000 cell cycles"));
p_av1 = plot([scatter(x=df_av_sp.kdam, y=df_av_sp[:,i], name="$(names(df_av_sp)[i])") for i in 2:10], Layout(xaxis_title="kdam", yaxis_title="mean molecule number", title="1000 cell cycles"));

p = [p_av p_av1]

sp = [3,4,6,8,9,11,13,14]

df_av_sp = df_av1[sp,:]

open("/home/hollie_hindley/Documents/stochastic_hybrid/cc_comp.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end