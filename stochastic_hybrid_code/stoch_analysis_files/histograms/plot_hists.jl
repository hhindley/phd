using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mainfolder = "/Users/s2257179/stoch_files/hists2"

load_hist_files(mainfolder)

mainfolder0 = "/Users/s2257179/stoch_files/hists"
foldernames0 = readdir(mainfolder0)
all_filenames0 = [readdir(joinpath(mainfolder0, foldernames0[i])) for i in eachindex(foldernames0)]

dfs0 = Dict{String, Vector{DataFrame}}()
for folder in foldernames0
    folderpath = joinpath(mainfolder0, folder)
    filenames = readdir(folderpath)
    for file in filenames
        prefix = split(file, ".")[1]
        filepath = joinpath(folderpath, file)
        df = DataFrame(CSV.File(filepath))
        df.bin = Array{Float64}.(JSON.parse.(df.bin))
        if haskey(dfs0, prefix)
            push!(dfs0[prefix], df)
        else
            dfs0[prefix] = [df]
        end
    end
end

threshold_vals1 = range(10,310,length=20)
threshold_vals2 = range(246,310,length=5)
threshold_vals0 = [threshold_vals1[1:15]; threshold_vals2]
titles = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]
titles0 = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]

# to plot one at a time
f = Figure()
ax = Axis(f[1,1], yscale=log10)
plot_hist(dfs["rh"][1], f[1,1])





rma_hist = plot_results("plot_hists", dfs, 2, 3, yscale=log10, species=:rm_a, xlabel="rm_a conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rmb_hist = plot_results("plot_hists", dfs, 2, 3, yscale=log10, species=:rm_b, xlabel="rm_b conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rmr_hist = plot_results("plot_hists", dfs, 2, 2, yscale=log10, species=:rm_r, xlabel="rm_r conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")

rtca_hist = plot_results("plot_hists", dfs, 2, 2, yscale=log10, species=:rtca, xlabel="rtca conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rtcb_hist = plot_results("plot_hists", dfs, 2, 2, species=:rtcb, yscale=log10, xlabel="rtcb conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rtcr_hist = plot_results("plot_hists", dfs, 2, 2, species=:rtcr, yscale=log10, xlabel="rtcr conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")

rh_hist = plot_results("plot_hists", dfs, 2, 2, yscale=identity, species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rd_hist = plot_results("plot_hists", dfs, 2, 2, yscale=log10, species=:rd, xlabel="rd conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")
rt_hist = plot_results("plot_hists", dfs, 2, 2, yscale=log10, species=:rt, xlabel="rt conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")


rma_hist0 = plot_results("plot_hists", dfs0, :rm_a, 5, 4, yscale=log10, xlabel="rm_a conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rmb_hist0 = plot_results("plot_hists", dfs0, :rm_b, 5, 4, yscale=log10, xlabel="rm_b conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rmr_hist0 = plot_results("plot_hists", dfs0, :rm_r, 5, 4, yscale=log10, xlabel="rm_r conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")

rtca_hist0 = plot_results("plot_hists", dfs0, :rtca, 5, 4, yscale=log10, xlabel="rtca conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rtcb_hist0 = plot_results("plot_hists", dfs0, :rtcb, 5, 4, yscale=log10, xlabel="rtcb conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rtcr_hist0 = plot_results("plot_hists", dfs0, :rtcr, 5, 4, yscale=log10, xlabel="rtcr conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")

rh_hist0 = plot_results("plot_hists", dfs0, :rh, 5, 4, yscale=log10, xlabel="rh conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rd_hist0 = plot_results("plot_hists", dfs0, :rd, 5, 4, yscale=log10, xlabel="rd conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")
rt_hist0 = plot_results("plot_hists", dfs0, :rt, 5, 4, yscale=log10, xlabel="rt conc", ylabel="frequency", titles=titles0, folder="thresh_plots/hists")

