using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files"

# need to fix the cases where time files arent properly labelled
df_times0, threshold_vals0, titles0, df_results0, df_reacts0 = LoadDataVars("new_thresh_vals_0507_nofloor_final_files", props=false,)
df_times1, threshold_vals1, titles1, df_results1, df_reacts1 = LoadDataVars("run_individually_0407_final_files", props=false)
df_times2, threshold_vals2, titles2, df_results2, df_reacts2 = LoadDataVars("LOW_thresh_0807_final_files", props=false)
df_times3, threshold_vals3, titles3, df_results3, df_reacts3 = LoadDataVars("test_0307_5_values_from_original_range_final_files", props=false)
df_times4, threshold_vals4, titles4, df_results4, df_reacts4 = LoadDataVars("new_thresh_vals_0507_final_files", props=false)
df_times5, threshold_vals5, titles5, df_results5, df_reacts5 = LoadDataVars("thresh_0107_noround_final_files", props=false)
df_times6, threshold_vals6, titles6, df_results6, df_reacts6 = LoadDataVars("thresh_0107_final_files", props=false)
df_times7, threshold_vals7, titles7, df_results7, df_reacts7 = LoadDataVars("thresh_1906_final_files", props=false)
df_times8, threshold_vals8, titles8, df_results8, df_reacts8 = LoadDataVars("thresh_test_arrow_files_14_06", props=false)
df_times9, threshold_vals9, titles9, df_results9, df_reacts9 = LoadDataVars("thresh_test_last5", props=false)
df_times10, threshold_vals10, titles10, df_results10, df_reacts10 = LoadDataVars("thresh_test_arrow", props=false)

f_times0 = plot_times(df_times0, "05/07 no flooring")
f_times1 = plot_times(df_times1, "04/07")
f_times2 = plot_times(df_times1, "08/07 LOW vals no flooring")
f_times3 = plot_times(df_times3, "03/07 just 5 values no flooring")
f_times4 = plot_times(df_times4, "05/07")
f_times5 = plot_times(df_times5, "05/07 new vals")
f_times6 = plot_times(df_times6, "01/07 no flooring")
f_times7 = plot_times(df_times7, "01/07 with flooring")
f_times8 = plot_times(df_times8, "19/06")
f_times9 = plot_times(df_times9, "no date last 5") # do these results fit on the end of results 10 below?
f_times10 = plot_times(df_times10, "no date")


tot_counts0 = prod_tot_count(df_reacts0)
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

f_tsc0 = plot_totstochcount(threshold_vals0, tot_counts0, "05/07 no flooring")
f_tsc1 = plot_totstochcount(threshold_vals1, tot_counts1, "04/07")
f_tsc2 = plot_totstochcount(threshold_vals2, tot_counts2, "08/07 LOW vals no flooring")
f_tsc3 = plot_totstochcount(threshold_vals3, tot_counts3, "03/07 just 5 values no flooring")
f_tsc4 = plot_totstochcount(threshold_vals4, tot_counts4, "05/07")
f_tsc5 = plot_totstochcount(threshold_vals5, tot_counts5, "05/07 new vals")
f_tsc6 = plot_totstochcount(threshold_vals6, tot_counts6, "01/07 no flooring")
f_tsc7 = plot_totstochcount(threshold_vals7, tot_counts7, "01/07 with flooring")
f_tsc8 = plot_totstochcount(threshold_vals8, tot_counts8, "19/06")
f_tsc9 = plot_totstochcount(threshold_vals9, tot_counts9, "no date last 5")
f_tsc10 = plot_totstochcount(threshold_vals10, tot_counts10, "no date")

# need to work out how to put titles into main function to make this complete
f_rtca0 = plot_results("plot_results", df_results0, length(threshold_vals0), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca1 = plot_results("plot_results", df_results1, length(threshold_vals1), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca2 = plot_results("plot_results", df_results2, length(threshold_vals2), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca3 = plot_results("plot_results", df_results3, length(threshold_vals3), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca4 = plot_results("plot_results", df_results4, length(threshold_vals4), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca5 = plot_results("plot_results", df_results5, length(threshold_vals5), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca6 = plot_results("plot_results", df_results6, length(threshold_vals6), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca7 = plot_results("plot_results", df_results7, length(threshold_vals7), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca8 = plot_results("plot_results", df_results8, length(threshold_vals8), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca9 = plot_results("plot_results", df_results9, length(threshold_vals9), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")
f_rtca10 = plot_results("plot_results", df_results10, length(threshold_vals10), species=:rtca, xlabel="time", ylabel="rtca")#, folder="testing_figs")




# need to finish this file so that for all files we also plot different reactions that happen stochastically and also propensities somehow

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



