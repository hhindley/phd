using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing/keyvals2")
folders_dict = Dict(filter(pair -> pair.first in [6,7,8,9], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false)

# plot one result
folder = 6; index = 1; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false, conc=true)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)

display(GLMakie.Screen(), f)

# display(GLMakie.Screen(), f_all)
display(GLMakie.Screen(), f_rib)
display(GLMakie.Screen(), f_mrna)
display(GLMakie.Screen(), f_prot)

DataInspector(f_rib)
DataInspector(f_mrna)
DataInspector(f_prot)


zoom = dict_results[folder][2][1:100000,:]

# two versions of determining the switching rate
# if rtca and rtcb are both above a certain threshold
# or if rtca and rtcb are both increasing then we are in the on state, if they are not changing or are decreasing then we are in the off state 
threshold=5
df, start_times, stop_times = switch_times(dict_results[7][2], threshold=threshold)
switch_rate_on, switch_rate_off = calc_switch_rate(df, start_times, stop_times)


# plot switch times
f = Figure()
ax = Axis(f[1,1], title="Threshold for on state = $threshold, switch rate = $(round(switch_rate_on, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax, df.time, df.rtca, label="rtca")
lines!(ax, df.time, df.rtcb, label="rtcb")
scatter!(ax, start_times, ones(length(start_times)), markersize=10, label="start on")
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10, label="stop on")
axislegend()
DataInspector(f)
display(GLMakie.Screen(), f)    

# plot on states
fdf = filter(row -> row.on == 1, df)
f_rtca = Figure()
ax2 = Axis(f_rtca[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax2, df.time, df.rtca, label="rtca")
scatter!(ax2, fdf.time, fdf.rtca, markersize=5, color=:pink, label="rtca on")
axislegend()
DataInspector(f_rtca)
display(GLMakie.Screen(), f_rtca) 

f_rtcb = Figure()
ax3 = Axis(f_rtcb[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax3, df.time, df.rtcb, label="rtcb")
scatter!(ax3, fdf.time, fdf.rtcb, markersize=5, color=:pink, label="rtcb on")
axislegend()
DataInspector(f_rtcb)
display(GLMakie.Screen(), f_rtcb) 


switch_vals_on = []
switch_vals_off = []
for i in eachindex(dict_results)
    println(i)
    switch_val_on, switch_val_off = get_all_switch_rates(i, threshold=threshold)
    push!(switch_vals_on, switch_val_on)
    push!(switch_vals_off, switch_val_off)
end


f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", title="threshold method", yscale=log10)
[lines!(ax, dict_kdamvals[i][:kdam], switch_vals_on[i], label="on→off $i") for i in eachindex(dict_kdamvals)]
[lines!(ax, dict_kdamvals[i][:kdam], switch_vals_off[i], label="off→on $i") for i in eachindex(dict_kdamvals)]

axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate", title="difference method", yscale=log10)
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on1, label="on→off")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off1, label="off→on")
axislegend()


fracs_on = []
fracs_off = []
for i in eachindex(switch_vals_on)
    println(i)
    frac_on, frac_off = calc_frac_times(switch_vals_on[i], switch_vals_off[i])
    push!(fracs_on, frac_on)
    push!(fracs_off, frac_off)
end

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, dict_kdamvals[i][:kdam], fracs_on[i], label="on state $i") for i in eachindex(dict_kdamvals)]
[lines!(ax, dict_kdamvals[i][:kdam], fracs_off[i], label="off state $i") for i in eachindex(dict_kdamvals)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)




# average conc for on and off states 
rtca_on, rtcb_on = get_av_conc_state(dict_results[6][1], on=true)
rtca_off, rtcb_off = get_av_conc_state(dict_results[6][1], on=false)

rtca_ons = []
rtcb_ons = []
rtca_offs = []
rtcb_offs = []
for i in eachindex(dict_results)
    println(i)
    rtca_on, rtcb_on, rtca_off, rtcb_off = get_all_av_conc(i)
    push!(rtca_ons, rtca_on)
    push!(rtcb_ons, rtcb_on)
    push!(rtca_offs, rtca_off)
    push!(rtcb_offs, rtcb_off)
end


rtc_on6 = [get_av_conc_state(dict_results[6][i], on=true) for i in eachindex(dict_results[6])]
rtc_off6 = [get_av_conc_state(dict_results[6][i], on=false) for i in eachindex(dict_results[6])]

rtc_on7 = [get_av_conc_state(dict_results[7][i], on=true) for i in eachindex(dict_results[7])]
rtc_off7 = [get_av_conc_state(dict_results[7][i], on=false) for i in eachindex(dict_results[7])]

rtc_on8 = [get_av_conc_state(dict_results[8][i], on=true) for i in eachindex(dict_results[8])]
rtc_off8 = [get_av_conc_state(dict_results[8][i], on=false) for i in eachindex(dict_results[8])]

rtc_on9 = [get_av_conc_state(dict_results[9][i], on=true) for i in eachindex(dict_results[9])]
rtc_off9 = [get_av_conc_state(dict_results[9][i], on=false) for i in eachindex(dict_results[9])]


rtca_on6 = [rtc_on6[i][1] for i in eachindex(rtc_on6)]
rtca_on7 = [rtc_on7[i][1] for i in eachindex(rtc_on7)]
rtca_on8 = [rtc_on8[i][1] for i in eachindex(rtc_on8)]
rtca_on9 = [rtc_on9[i][1] for i in eachindex(rtc_on9)]

rtcb_on6 = [rtc_on6[i][2] for i in eachindex(rtc_on6)]
rtcb_on7 = [rtc_on7[i][2] for i in eachindex(rtc_on7)]
rtcb_on8 = [rtc_on8[i][2] for i in eachindex(rtc_on8)]
rtcb_on9 = [rtc_on9[i][2] for i in eachindex(rtc_on9)]

rtca_off6 = [rtc_off6[i][1] for i in eachindex(rtc_off6)]
rtca_off7 = [rtc_off7[i][1] for i in eachindex(rtc_off7)]
rtca_off8 = [rtc_off8[i][1] for i in eachindex(rtc_off8)]
rtca_off9 = [rtc_off9[i][1] for i in eachindex(rtc_off9)]

rtcb_off6 = [rtc_off6[i][2] for i in eachindex(rtc_off6)]
rtcb_off7 = [rtc_off7[i][2] for i in eachindex(rtc_off7)]
rtcb_off8 = [rtc_off8[i][2] for i in eachindex(rtc_off8)]
rtcb_off9 = [rtc_off9[i][2] for i in eachindex(rtc_off9)]

normalized_fracs_on = [clamp.(fracs_on[i], 0, 1) for i in eachindex(fracs_on)]
normalized_fracs_off = [clamp.(fracs_off[i], 0, 1) for i in eachindex(fracs_off)]

f = Figure()
ax = Axis(f[1,1])#, yscale=log10)
[scatter!(ax, dict_kdamvals[i][:kdam], rtca_ons[i], color=normalized_fracs_on[i]) for i in eachindex(rtca_ons)]
[scatter!(ax, dict_kdamvals[i][:kdam], rtca_offs[i], color=normalized_fracs_off[i]) for i in eachindex(rtca_offs)]
Colorbar(f[1, 2], limits = (0,1), colormap = :viridis)

min_on = minimum([minimum(fracs_on[i]) for i in eachindex(fracs_on)])
max_on = maximum([maximum(fracs_on[i]) for i in eachindex(fracs_on)])
min_off = minimum([minimum(fracs_off[i]) for i in eachindex(fracs_on)])
max_off = maximum([maximum(fracs_off[i]) for i in eachindex(fracs_on)])



df = DataFrame(:rtca_ons=>rtca_ons, :rtcb_ons=>rtcb_ons, :rtca_offs=>rtca_offs, :rtcb_offs=>rtcb_offs, :fracs_on=>fracs_on, :fracs_off=>fracs_off)
df_kdam_vals = DataFrame(:kdam=>dict_kdamvals[6][:kdam])

CSV.write("/Users/s2257179/Desktop/res.csv", df)
CSV.write("/Users/s2257179/Desktop/kdam.csv", df_kdam_vals)