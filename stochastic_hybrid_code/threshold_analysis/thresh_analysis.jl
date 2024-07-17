using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files"

folders = ["new_thresh_vals_0507_nofloor_final_files", #1
           "run_individually_0407_final_files", #2
           "LOW_thresh_0807_final_files", #3
           "test_0307_5_values_from_original_range_final_files", #4
           "new_thresh_vals_0507_final_files", #5
           "thresh_0107_noround_final_files", #6
           "thresh_0107_final_files", #7
           "thresh_1906_final_files", #8
           "thresh_test_arrow_files_14_06", #9
           "thresh_analysis_fixed_xi_final_files", #10
           "stoch_division_thresh_analysis_final_files" #11
        ]

df_times1, threshold_vals1, titles1, df_results1, df_reacts1, df_props1 = LoadDataVars(folders[1], props=false);
df_times2, threshold_vals2, titles2, df_results2, df_reacts2, df_props2 = LoadDataVars(folders[2], props=false);
df_times3, threshold_vals3, titles3, df_results3, df_reacts3, df_props3 = LoadDataVars(folders[3], props=true);
df_times4, threshold_vals4, titles4, df_results4, df_reacts4, df_props4 = LoadDataVars(folders[4], props=false);
df_times5, threshold_vals5, titles5, df_results5, df_reacts5, df_props5 = LoadDataVars(folders[5], props=false);
df_times6, threshold_vals6, titles6, df_results6, df_reacts6, df_props6 = LoadDataVars(folders[6], props=false);
df_times7, threshold_vals7, titles7, df_results7, df_reacts7, df_props7 = LoadDataVars(folders[7], props=false);
df_times8, threshold_vals8, titles8, df_results8, df_reacts8, df_props8 = LoadDataVars(folders[8], props=false); 
df_times9, threshold_vals9, titles9, df_results9, df_reacts9, df_props9 = LoadDataVars(folders[9], props=false, timefilepath="/Users/s2257179/stoch_files/thresh_1906_final_files/thresh_1906_times.csv");
df_times10, threshold_vals10, titles10, df_results10, df_reacts10, df_props10 = LoadDataVars(folders[10], props=false);
df_times11, threshold_vals11, titles11, df_results11, df_reacts11, df_props11 = LoadDataVars(folders[11], props=true);

f_times1 = plot_times(df_times1, "05/07 no flooring", folder=folders[1])
f_times2 = plot_times(df_times2, "04/07", folder=folders[2])
f_times3 = plot_times(df_times3, "08/07 LOW vals no flooring", folder=folders[3])
f_times4 = plot_times(df_times4, "03/07 just 5 values no flooring", folder=folders[4])
f_times5 = plot_times(df_times5, "05/07", folder=folders[5])
f_times6 = plot_times(df_times6, "05/07 new vals", folder=folders[6])
f_times7 = plot_times(df_times7, "01/07 no flooring", folder=folders[7])
f_times8 = plot_times(df_times8, "01/07 with flooring", folder=folders[8])
f_times9 = plot_times(df_times9, "19/06 no time file", folder=folders[9])
f_times10 = plot_times(df_times10, "12/07 xi fixed", folder=folders[10])
f_times11 = plot_times(df_times11, "stoch_division_thresh_analysis", folder=folders[11])

tot_counts1 = prod_tot_count(df_reacts1)
tot_counts2 = prod_tot_count(df_reacts2)
tot_counts3 = prod_tot_count(df_reacts3)
tot_counts4 = prod_tot_count(df_reacts4)
tot_counts5 = prod_tot_count(df_reacts5)
tot_counts6 = prod_tot_count(df_reacts6)
tot_counts7 = prod_tot_count(df_reacts7)
tot_counts8 = prod_tot_count(df_reacts8)
tot_counts9 = prod_tot_count(df_reacts9)
tot_counts10 = prod_tot_count(df_reacts10)
tot_counts11 = prod_tot_count(df_reacts11)

f_tsc1 = plot_totstochcount(threshold_vals1, tot_counts1, "05/07 no flooring", folder=folders[1])
f_tsc2 = plot_totstochcount(threshold_vals2, tot_counts2, "04/07", folder=folders[2])
f_tsc3 = plot_totstochcount(threshold_vals3, tot_counts3, "08/07 LOW vals no flooring", folder=folders[3])
f_tsc4 = plot_totstochcount(threshold_vals4, tot_counts4, "03/07 just 5 values no flooring", folder=folders[4])
f_tsc5 = plot_totstochcount(threshold_vals5, tot_counts5, "05/07", folder=folders[5])
f_tsc6 = plot_totstochcount(threshold_vals6, tot_counts6, "05/07 new vals", folder=folders[6])
f_tsc7 = plot_totstochcount(threshold_vals7, tot_counts7, "01/07 no flooring", folder=folders[7])
f_tsc8 = plot_totstochcount(threshold_vals8, tot_counts8, "01/07 with flooring", folder=folders[8])
f_tsc9 = plot_totstochcount(threshold_vals9, tot_counts9, "19/06", folder=folders[9])
f_tsc10 = plot_totstochcount(threshold_vals10, tot_counts10, "12/07 xi fixed", folder=folders[10])
f_tsc11 = plot_totstochcount(threshold_vals11, tot_counts11, "stoch division thresh analysis", folder=folders[11])

hists1 = load_hist_files(joinpath(mount_path, folders[1], "hists"))
hists2 = load_hist_files(joinpath(mount_path, folders[2], "hists"))
hists3 = load_hist_files(joinpath(mount_path, folders[3], "hists"))
hists4 = load_hist_files(joinpath(mount_path, folders[4], "hists"))
hists5 = load_hist_files(joinpath(mount_path, folders[5], "hists"))
hists6 = load_hist_files(joinpath(mount_path, folders[6], "hists"))
hists7 = load_hist_files(joinpath(mount_path, folders[7], "hists"))
hists8 = load_hist_files(joinpath(mount_path, folders[8], "hists"))
hists9 = load_hist_files(joinpath(mount_path, folders[9], "hists"))
hists10 = load_hist_files(joinpath(mount_path, folders[10], "hists"))
hists11 = load_hist_files(joinpath(mount_path, folders[11], "hists"))


# here put the folders you want to plot/save
for i in all_species
    f_rtca1 = plot_results("plot_results", df_results11, length(threshold_vals11), species=i, xlabel="time", ylabel="$i", titles=titles11, size=(1000,650), folder=folders[11]);
    h_rh1 = plot_results("plot_hists", hists11, length(threshold_vals11), species=i, xlabel="$i", ylabel="frequency", titles=titles11, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[11]);

end

f_rtca1 = plot_results("plot_results", df_results11, length(threshold_vals11), species=:rtca, xlabel="time", ylabel="rtca", titles=titles11, size=(1000,650), folder=folders[11], hidelabels=[false, false], linkaxes=false);


f_sr1 = plot_results("plot_stoch_reacts", df_reacts1, length(threshold_vals1), xlabel="reaction", ylabel="count", titles=titles1, folder=folders[1], size=(1000,650))
f_sr2 = plot_results("plot_stoch_reacts", df_reacts2, length(threshold_vals2), xlabel="reaction", ylabel="count", titles=titles2, folder=folders[2], size=(1000,650))
f_sr3 = plot_results("plot_stoch_reacts", df_reacts3, length(threshold_vals3), xlabel="reaction", ylabel="count", titles=titles3, folder=folders[3], size=(1000,650))
f_sr4 = plot_results("plot_stoch_reacts", df_reacts4, length(threshold_vals4), xlabel="reaction", ylabel="count", titles=titles4, folder=folders[4], size=(1000,650))
f_sr5 = plot_results("plot_stoch_reacts", df_reacts5, length(threshold_vals5), xlabel="reaction", ylabel="count", titles=titles5, folder=folders[5], size=(1000,650))
f_sr6 = plot_results("plot_stoch_reacts", df_reacts6, length(threshold_vals6), xlabel="reaction", ylabel="count", titles=titles6, folder=folders[6], size=(1000,650))
f_sr7 = plot_results("plot_stoch_reacts", df_reacts7, length(threshold_vals7), xlabel="reaction", ylabel="count", titles=titles7, folder=folders[7], size=(1000,650))
f_sr8 = plot_results("plot_stoch_reacts", df_reacts8, length(threshold_vals8), xlabel="reaction", ylabel="count", titles=titles8, folder=folders[8], size=(1000,650))
f_sr9 = plot_results("plot_stoch_reacts", df_reacts9, length(threshold_vals9), xlabel="reaction", ylabel="count", titles=titles9, folder=folders[9], size=(1000,650))
f_sr10 = plot_results("plot_stoch_reacts", df_reacts10, length(threshold_vals10), xlabel="reaction", ylabel="count", titles=titles10, folder=folders[10], size=(1000,650))
f_sr11 = plot_results("plot_stoch_reacts", df_reacts11, length(threshold_vals11), xlabel="reaction", ylabel="count", titles=titles11, folder=folders[11], size=(1000,650))



f1 = plot_prop(df_results11, df_props11, 1, "threshold_$(threshold_vals11[1])", threshold_vals11, maxval=3000, folder=folders[11])
f2 = plot_prop(df_results11, df_props11, 2, "threshold_$(threshold_vals11[2])", threshold_vals11, folder=folders[11])
f3 = plot_prop(df_results11, df_props11, 3, "threshold_$(threshold_vals11[3])", threshold_vals11, folder=folders[11])
f4 = plot_prop(df_results11, df_props11, 4, "threshold_$(threshold_vals11[4])", threshold_vals11, maxval=3000, folder=folders[11])
f5 = plot_prop(df_results11, df_props11, 5, "threshold_$(threshold_vals11[5])", threshold_vals11, maxval=3000, folder=folders[11])


display(f1)
f3_3 = plot_prop(df_results3, df_props3, 2, "threshold_$(threshold_vals3[2])", threshold_vals3)


df_results11[1].rm_a

h2 = hist(df_results11[2].rtca[100000:200000])
h3 = hist!(df_results11[2].rtca[200000:300000])
h4 = hist!(df_results11[2].rtca[300000:400000])
h5 = hist!(df_results11[2].rtca[400000:500000])
h6 = hist!(df_results11[2].rtca[500000:600000])
h7 = hist!(df_results11[2].rtca[700000:800000])
h8 = hist!(df_results11[2].rtca[800000:end])



10000*log(2)/lam_val

df_results11[2][1:end,:rm_a]