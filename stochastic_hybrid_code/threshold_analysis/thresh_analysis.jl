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
           "thresh_test_last5" #10
        ]

df_times1, threshold_vals1, titles1, df_results1, df_reacts1 = LoadDataVars(folders[1], props=false);
df_times2, threshold_vals2, titles2, df_results2, df_reacts2 = LoadDataVars(folders[2], props=false);
df_times3, threshold_vals3, titles3, df_results3, df_reacts3 = LoadDataVars(folders[3], props=false);
df_times4, threshold_vals4, titles4, df_results4, df_reacts4 = LoadDataVars(folders[4], props=false);
df_times5, threshold_vals5, titles5, df_results5, df_reacts5 = LoadDataVars(folders[5], props=false);
df_times6, threshold_vals6, titles6, df_results6, df_reacts6 = LoadDataVars(folders[6], props=false);
df_times7, threshold_vals7, titles7, df_results7, df_reacts7 = LoadDataVars(folders[7], props=false);
df_times8, threshold_vals8, titles8, df_results8, df_reacts8 = LoadDataVars(folders[8], props=false); 
df_times9, threshold_vals9, titles9, df_results9, df_reacts9 = LoadDataVars(folders[9], props=false, timefilepath="/Users/s2257179/stoch_files/thresh_1906_final_files/thresh_1906_times.csv");

f_times1 = plot_times(df_times1, "05/07 no flooring", folder=folders[1])
f_times2 = plot_times(df_times2, "04/07", folder=folders[2])
f_times3 = plot_times(df_times3, "08/07 LOW vals no flooring", folder=folders[3])
f_times4 = plot_times(df_times4, "03/07 just 5 values no flooring", folder=folders[4])
f_times5 = plot_times(df_times5, "05/07", folder=folders[5])
f_times6 = plot_times(df_times6, "05/07 new vals", folder=folders[6])
f_times7 = plot_times(df_times7, "01/07 no flooring", folder=folders[7])
f_times8 = plot_times(df_times8, "01/07 with flooring", folder=folders[8])
f_times9 = plot_times(df_times9, "19/06 no time file", folder=folders[9])

tot_counts1 = prod_tot_count(df_reacts1)
tot_counts2 = prod_tot_count(df_reacts2)
tot_counts3 = prod_tot_count(df_reacts3)
tot_counts4 = prod_tot_count(df_reacts4)
tot_counts5 = prod_tot_count(df_reacts5)
tot_counts6 = prod_tot_count(df_reacts6)
tot_counts7 = prod_tot_count(df_reacts7)
tot_counts8 = prod_tot_count(df_reacts8)
tot_counts9 = prod_tot_count(df_reacts9)

f_tsc1 = plot_totstochcount(threshold_vals1, tot_counts1, "05/07 no flooring", folder=folders[1])
f_tsc2 = plot_totstochcount(threshold_vals2, tot_counts2, "04/07", folder=folders[2])
f_tsc3 = plot_totstochcount(threshold_vals3, tot_counts3, "08/07 LOW vals no flooring", folder=folders[3])
f_tsc4 = plot_totstochcount(threshold_vals4, tot_counts4, "03/07 just 5 values no flooring", folder=folders[4])
f_tsc5 = plot_totstochcount(threshold_vals5, tot_counts5, "05/07", folder=folders[5])
f_tsc6 = plot_totstochcount(threshold_vals6, tot_counts6, "05/07 new vals", folder=folders[6])
f_tsc7 = plot_totstochcount(threshold_vals7, tot_counts7, "01/07 no flooring", folder=folders[7])
f_tsc8 = plot_totstochcount(threshold_vals8, tot_counts8, "01/07 with flooring", folder=folders[8])
f_tsc9 = plot_totstochcount(threshold_vals9, tot_counts9, "19/06", folder=folders[9])


for i in all_species
    f_rtca1 = plot_results("plot_results", df_results1, length(threshold_vals1), species=i, xlabel="time", ylabel="$i", titles=titles1, size=(1000,650), folder=folders[1]);
    f_rtca2 = plot_results("plot_results", df_results2, length(threshold_vals2), species=i, xlabel="time", ylabel="$i", titles=titles2, size=(1000,650), folder=folders[2]);
    f_rtca3 = plot_results("plot_results", df_results3, length(threshold_vals3), species=i, xlabel="time", ylabel="$i", titles=titles3, size=(1000,650), folder=folders[3]);
    f_rtca4 = plot_results("plot_results", df_results4, length(threshold_vals4), species=i, xlabel="time", ylabel="$i", titles=titles4, size=(1000,650), folder=folders[4]);
    f_rtca5 = plot_results("plot_results", df_results5, length(threshold_vals5), species=i, xlabel="time", ylabel="$i", titles=titles5, size=(1000,650), folder=folders[5]);
    f_rtca6 = plot_results("plot_results", df_results6, length(threshold_vals6), species=i, xlabel="time", ylabel="$i", titles=titles6, size=(1000,650), folder=folders[6]);
    f_rtca7 = plot_results("plot_results", df_results7, length(threshold_vals7), species=i, xlabel="time", ylabel="$i", titles=titles7, size=(1000,650), folder=folders[7]);
    f_rtca8 = plot_results("plot_results", df_results8, length(threshold_vals8), species=i, xlabel="time", ylabel="$i", titles=titles8, size=(1000,650), folder=folders[8]);
    f_rtca9 = plot_results("plot_results", df_results9, length(threshold_vals9), species=i, xlabel="time", ylabel="$i", titles=titles9, size=(1000,650), folder=folders[9]);

    h_rh1 = plot_results("plot_hists", hists1, length(threshold_vals1), species=i, xlabel="$i", ylabel="frequency", titles=titles1, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[1]);
    h_rh2 = plot_results("plot_hists", hists2, length(threshold_vals2), species=i, xlabel="$i", ylabel="frequency", titles=titles2, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[2]);
    h_rh3 = plot_results("plot_hists", hists3, length(threshold_vals3), species=i, xlabel="$i", ylabel="frequency", titles=titles3, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[3]);
    h_rh4 = plot_results("plot_hists", hists4, length(threshold_vals4), species=i, xlabel="$i", ylabel="frequency", titles=titles4, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[4]);
    h_rh5 = plot_results("plot_hists", hists5, length(threshold_vals5), species=i, xlabel="$i", ylabel="frequency", titles=titles5, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[5]);
    h_rh6 = plot_results("plot_hists", hists6, length(threshold_vals6), species=i, xlabel="$i", ylabel="frequency", titles=titles6, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[6]);
    h_rh7 = plot_results("plot_hists", hists7, length(threshold_vals7), species=i, xlabel="$i", ylabel="frequency", titles=titles7, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[7]);
    h_rh8 = plot_results("plot_hists", hists8, length(threshold_vals8), species=i, xlabel="$i", ylabel="frequency", titles=titles8, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[8]);
    h_rh9 = plot_results("plot_hists", hists9, length(threshold_vals9), species=i, xlabel="$i", ylabel="frequency", titles=titles9, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[9]);

end


hists1 = load_hist_files(joinpath(mount_path, folders[1], "hists"))
hists2 = load_hist_files(joinpath(mount_path, folders[2], "hists"))
hists3 = load_hist_files(joinpath(mount_path, folders[3], "hists"))
hists4 = load_hist_files(joinpath(mount_path, folders[4], "hists"))
hists5 = load_hist_files(joinpath(mount_path, folders[5], "hists"))
hists6 = load_hist_files(joinpath(mount_path, folders[6], "hists"))
hists7 = load_hist_files(joinpath(mount_path, folders[7], "hists"))
hists8 = load_hist_files(joinpath(mount_path, folders[8], "hists"))
hists9 = load_hist_files(joinpath(mount_path, folders[9], "hists"))


h_rh1 = plot_results("plot_hists", hists1, length(threshold_vals1), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles1, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[1])
h_rh2 = plot_results("plot_hists", hists2, length(threshold_vals2), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles2, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[2])
h_rh3 = plot_results("plot_hists", hists3, length(threshold_vals3), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles3, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[3])
h_rh4 = plot_results("plot_hists", hists4, length(threshold_vals4), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles4, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[4])
h_rh5 = plot_results("plot_hists", hists5, length(threshold_vals5), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles5, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[5])
h_rh6 = plot_results("plot_hists", hists6, length(threshold_vals6), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles6, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[6])
h_rh7 = plot_results("plot_hists", hists7, length(threshold_vals7), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles7, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[7])
h_rh8 = plot_results("plot_hists", hists8, length(threshold_vals8), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles8, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[8])
h_rh9 = plot_results("plot_hists", hists9, length(threshold_vals9), species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles9, hidelabels=[false, false], linkaxes=false, size=(1000,650), folder=folders[9])





p_totprop1 = plot_results("plot_results", df_results1, length(threshold_vals1), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles1, size=(1000,650), folder=folders[1]);
p_totprop2 = plot_results("plot_results", df_results2, length(threshold_vals2), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles2, size=(1000,650), folder=folders[2]);
p_totprop3 = plot_results("plot_results", df_results3, length(threshold_vals3), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles3, size=(1000,650), folder=folders[3]);
p_totprop4 = plot_results("plot_results", df_results4, length(threshold_vals4), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles4, size=(1000,650), folder=folders[4]);
p_totprop5 = plot_results("plot_results", df_results5, length(threshold_vals5), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles5, size=(1000,650), folder=folders[5]);
p_totprop6 = plot_results("plot_results", df_results6, length(threshold_vals6), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles6, size=(1000,650), folder=folders[6]);
p_totprop7 = plot_results("plot_results", df_results7, length(threshold_vals7), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles7, size=(1000,650), folder=folders[7]);
p_totprop8 = plot_results("plot_results", df_results8, length(threshold_vals8), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles8, size=(1000,650), folder=folders[8]);
p_totprop9 = plot_results("plot_results", df_results9, length(threshold_vals9), species=:totprop, xlabel="time", ylabel="tot prop", titles=titles9, size=(1000,650), folder=folders[9]);




df_reacts
f = Figure()
ax = Axis(f[1,1],xlabel="reactions", ylabel="reaction count", title="thresh_val", xticks=(1:13, react_names_str), xticklabelrotation=45)
barplot!(df_reacts[1].event, df_reacts[1].count)


f = plot_results("plot_stoch_reacts", df_reacts, 5, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, true], linkaxes=true, species="stoch_reacts")#, folder="thresh_plots")
display(f)


dfs = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts, rn, threshold_vals)
    if !haskey(dfs, rn)
        dfs[rn] = df
    end
end

titles_reacts = ["$i" for i in react_names[1:end-1]]

p1 = plot_results("plot_individual_reacts", dfs, 12, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots", species="individual_reacts")





f_rtca = plot_results("plot_results", df_results, 5, species=:rtca, xlabel="time", ylabel="rtca", titles=titles, folder="testing_figs")



f1 = plot_prop(df_results, df_props, 1, "thresh_plots2/props", "threshold_$(threshold_vals[1])", threshold_vals, 31856.296174439733)
f2 = plot_prop(df_results, df_props, 2, "thresh_plots2/props", "threshold_$(threshold_vals[2])", threshold_vals, 31856.296174439733)
f3 = plot_prop(df_results, df_props, 3, "thresh_plots2/props", "threshold_$(threshold_vals[3])", threshold_vals, 31856.296174439733)
f4 = plot_prop(df_results, df_props, 4, "thresh_plots2/props", "threshold_$(threshold_vals[4])", threshold_vals, 31856.296174439733)
f5 = plot_prop(df_results, df_props, 5, "thresh_plots2/props", "threshold_$(threshold_vals[5])", threshold_vals, 31856.296174439733)


f = plot_props(df_results, df_props, 5, threshold_vals)



