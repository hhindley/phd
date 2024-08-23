using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"
# mount_path = "/Users/s2257179/stoch_files/hysteresis/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [9], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

# plot one result
folder = 9; index = 1; species = "rh"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)


zoom = dict_results[folder][1][1000:5000,:]
f_all = plot_results("plot_results", zoom, 1, folders_dict[folder], species=all_species, titles=[dict_titles[folder][index]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true, yscale=log10)
f_rib = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_mrna = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_prot = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)

# display(GLMakie.Screen(), f_all)
display(GLMakie.Screen(), f_rib)
display(GLMakie.Screen(), f_mrna)
display(GLMakie.Screen(), f_prot)

DataInspector(f_all)
DataInspector(f_mrna)
DataInspector(f_prot)

v1 = lines(dict_results[folder][1].time, dict_results[folder][1].volume)
dict_results[folder][1].volume
v7 = lines(dict_results[folder][7].time, dict_results[folder][7].volume)

display(GLMakie.Screen(), v1)
display(GLMakie.Screen(), v7)

# if you want to plot with log scale then need to remove the zeros and replace with a very small value
epsilon = 1e-5
dict_results[folder][index] = replace_zeros_in_dataframe(dict_results[folder][index])
nodam = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=all_species, titles=[dict_titles[folder][index]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, yscale=identity, conc=false)
nodam_conc = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=all_species, titles=[dict_titles[folder][index]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, yscale=identity, conc=true)

display(GLMakie.Screen(), nodam)
display(GLMakie.Screen(), nodam_conc)

DataInspector(nodam)
DataInspector(nodam_conc)


dict_kdamvals[folder][:kdam]

f_rib1 = plot_results("plot_results", dict_results[folder][1], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][1]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna1 = plot_results("plot_results", dict_results[folder][1], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][1]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot1 = plot_results("plot_results", dict_results[folder][1], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][1]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib2 = plot_results("plot_results", dict_results[folder][2], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][2]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna2 = plot_results("plot_results", dict_results[folder][2], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][2]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot2 = plot_results("plot_results", dict_results[folder][2], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][2]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib3 = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna3 = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot3 = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib4 = plot_results("plot_results", dict_results[folder][4], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][4]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna4 = plot_results("plot_results", dict_results[folder][4], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][4]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot4 = plot_results("plot_results", dict_results[folder][4], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][4]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib5 = plot_results("plot_results", dict_results[folder][5], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][5]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna5 = plot_results("plot_results", dict_results[folder][5], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][5]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot5 = plot_results("plot_results", dict_results[folder][5], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][5]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib6 = plot_results("plot_results", dict_results[folder][6], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][6]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna6 = plot_results("plot_results", dict_results[folder][6], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][6]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot6 = plot_results("plot_results", dict_results[folder][6], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][6]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib7 = plot_results("plot_results", dict_results[folder][7], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][7]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna7 = plot_results("plot_results", dict_results[folder][7], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][7]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot7 = plot_results("plot_results", dict_results[folder][7], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][7]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

display(GLMakie.Screen(), f_rib5)
display(GLMakie.Screen(), f_mrna1)
display(GLMakie.Screen(), f_prot1)

DataInspector(f_rib1)
DataInspector(f_mrna1)
DataInspector(f_prot1)


display(GLMakie.Screen(), f_rib_on)
display(GLMakie.Screen(), f_rib_off)
display(GLMakie.Screen(), f_rib_off2)

# working out which regions of dataframe are rtc on and calculating percentage of total dataframe 
dict_results[1][1].rm_a
zoom = dict_results[folder][1][1510000:1800000,:]
f_rib_on_zoom = plot_results("plot_results", dict_results[1][1], 1, folders_dict[1], species=[:rh, :rd, :rt], titles=[dict_titles[1][1]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
display(GLMakie.Screen(), f_rib_on_zoom)

DataInspector(f_rib_on_zoom)

function calc_perc_on(res, title)
    df = res[!,[:time, :rh, :rd, :rt]]
    fdf = filter(row -> row.rh > row.rd && row.rt > 10, df)

    f = Figure()
    ax = Axis(f[1,1], title=title)
    [lines!(ax, df.time, get_column(df, s)) for s in [:rh, :rd, :rt]]
    
    f2 = Figure()
    ax2 = Axis(f2[1,1], title=title)
    [lines!(ax2, df.time, get_column(df, s)) for s in [:rh, :rd, :rt]]
    [lines!(ax2, fdf.time, get_column(fdf, s)) for s in [:rh, :rd, :rt]]

    return (length(fdf.time)/length(df.time))*100, f, f2
end

perc, f, f2 = calc_perc_on(zoom, dict_titles[1][end])
display(GLMakie.Screen(), f)
display(GLMakie.Screen(), f2)

zoom1 = dict_results[2][4][1510000:1800000,:]
perc1, f1, f21 = calc_perc_on(zoom1, dict_titles[2][4])
display(GLMakie.Screen(), f1)
display(GLMakie.Screen(), f21)


zoom2 = dict_results[2][16][1510000:1800000,:]
perc2, f2, f22 = calc_perc_on(zoom2, dict_titles[2][16])
display(GLMakie.Screen(), f2)
display(GLMakie.Screen(), f22)