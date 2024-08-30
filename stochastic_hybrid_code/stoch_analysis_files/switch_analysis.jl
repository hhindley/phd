using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [6,7,8,9], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)

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



switch_vals_on6, switch_vals_off6 = get_all_switch_rates(6, threshold=threshold)
switch_vals_on7, switch_vals_off7 = get_all_switch_rates(7, threshold=threshold)
switch_vals_on8, switch_vals_off8 = get_all_switch_rates(8, threshold=threshold)
switch_vals_on9, switch_vals_off9 = get_all_switch_rates(9, threshold=threshold)


f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", title="threshold method", yscale=log10)
lines!(ax, dict_kdamvals[6][:kdam], switch_vals_on6, label="on→off 6")
lines!(ax, dict_kdamvals[6][:kdam], switch_vals_off6, label="off→on 6")
lines!(ax, dict_kdamvals[7][:kdam], switch_vals_on7, label="on→off 7")
lines!(ax, dict_kdamvals[7][:kdam], switch_vals_off7, label="off→on 7")
lines!(ax, dict_kdamvals[8][:kdam], switch_vals_on8, label="on→off 8")
lines!(ax, dict_kdamvals[8][:kdam], switch_vals_off8, label="off→on 8")
lines!(ax, dict_kdamvals[9][:kdam], switch_vals_on9, label="on→off 9")
lines!(ax, dict_kdamvals[9][:kdam], switch_vals_off9, label="off→on 9")
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate", title="difference method", yscale=log10)
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on1, label="on→off")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off1, label="off→on")
axislegend()


frac_on6, frac_off6 = calc_frac_times(switch_vals_on6, switch_vals_off6)
frac_on7, frac_off7 = calc_frac_times(switch_vals_on7, switch_vals_off7)
frac_on8, frac_off8 = calc_frac_times(switch_vals_on8, switch_vals_off8)
frac_on9, frac_off9 = calc_frac_times(switch_vals_on9, switch_vals_off9)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
lines!(ax, dict_kdamvals[6][:kdam], frac_on6, label="on state6")
lines!(ax, dict_kdamvals[6][:kdam], frac_off6, label="off state6")
lines!(ax, dict_kdamvals[7][:kdam], frac_on7, label="on state7")
lines!(ax, dict_kdamvals[7][:kdam], frac_off7, label="off state7")
lines!(ax, dict_kdamvals[8][:kdam], frac_on8, label="on state8")
lines!(ax, dict_kdamvals[8][:kdam], frac_off8, label="off state8")
lines!(ax, dict_kdamvals[9][:kdam], frac_on9, label="on state9")
lines!(ax, dict_kdamvals[9][:kdam], frac_off9, label="off state9")
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)
frac_on6+frac_off6
dict_reacts[6][1]


# average conc for on and off states 
rtca_on, rtcb_on = get_av_conc_state(dict_results[6][1], on=true)
rtca_off, rtcb_off = get_av_conc_state(dict_results[6][1], on=false)

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

f = Figure()
ax = Axis(f[1,1])#, yscale=log10)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_on6, color=frac_on6)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_on7, color=frac_on7)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_on8, color=frac_on8)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_on9, color=frac_on9)

scatter!(ax, dict_kdamvals[6][:kdam], rtca_off6, color=frac_off6)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_off7, color=frac_off7)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_off8, color=frac_off8)
scatter!(ax, dict_kdamvals[6][:kdam], rtca_off9, color=frac_off9)