using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_props0 = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/props")
df_props = load_files("/Users/s2257179/stoch_files/thresh_1906_final_files/props")
df_results0 = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/results")
df_results = load_files("/Users/s2257179/stoch_files/thresh_1906_final_files/results")

color_list = [
    RGB(0.121, 0.466, 0.705),  # Blue
    RGB(1.0, 0.498, 0.054),    # Orange
    RGB(0.172, 0.627, 0.172),  # Green
    RGB(0.839, 0.152, 0.156),  # Red
    RGB(0.580, 0.403, 0.741),  # Purple
    RGB(0.549, 0.337, 0.294),  # Brown
    RGB(0.890, 0.466, 0.760),  # Pink
    RGB(0.498, 0.498, 0.498),  # Gray
    RGB(22/255, 219/255, 180/255),  # Olive
    RGB(13/255, 48/255, 222/255),  # Cyan
    RGB(57/255, 235/255, 21/255),        # Black
    RGB(255/255, 219/255, 38/255),        # Yellow
    RGB(1.0, 0.0, 1.0)         # Magenta
]

threshold_vals1 = range(10,310,length=20)
threshold_vals2 = range(246,310,length=5)
threshold_vals0 = [threshold_vals1[1:15]; threshold_vals2]


# f = create_subplots("plot_props", 5, 4)
# [iscatter!(f[1,1], df_results[1].time, get_column(df_props[1], react_names[i]), label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
# display(f)

# plot_props(df_results)

# f = Figure(size=(1450, 900))
# ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", limits=(nothing,(0, 31856.296174439733)))
# [iscatter!(f[1,1],df_results[1].time, get_column(df_props[1], react_names[i]), label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]

# axislegend(ax, framevisible=true)
# lines!(range(minimum(df_results[1].time), maximum(df_results[1].time), length=2), [threshold_vals1[1], threshold_vals1[1]], linewidth=4, color=RGB(0.0, 0.0, 0.0))

# f = plot_props_subplots()
# [iscatter!(f[1,1], df_results[1].time, get_column(df_props[1], react_names[i]), label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]


# add_subplots(f, "plot_props", df_results, 5, 4)

# display(f)

# for df in df_props
#     println(maximum(maximum(column) for column in eachcol(df)))
# end


# f = Figure()
# ax = Axis(f[1,1])
# for r in names(df_props[1])
#     iscatter!(ax, df_results[1].time, df_props[1][:,r])
# end
# ylims!(ax, 0,maximum([maximum(i) for i in eachcol(df_props[1])]))



# f = Figure(size=(1450, 900))

# for j in 1:4
#     for i in 1:5
#         data_ind = i + 5 * (j - 1)
#         ax = Axis(f[i,j], limits=(nothing,(0, 31856.296174439733)))
#         [iscatter!(ax,df_results[data_ind].time, get_column(df_props[data_ind], react_names[r]), label="$(react_names[r])", color=color_list[r]) for r in eachindex(react_names[1:end-1])]
#     end
# end


function plot_prop(df_results, df_props, res_ind, savein, title, threshold_vals, max_val)
    f = Figure()
    ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", title=title, limits=(nothing,(0, max_val)))
    [iscatter!(ax,df_results[res_ind].time, get_column(df_props[res_ind], react_names[i]), label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
    axislegend(ax, framevisible=true)
    lines!(range(minimum(df_results[res_ind].time), maximum(df_results[res_ind].time), length=2), [threshold_vals[res_ind], threshold_vals[res_ind]], linewidth=4, color=RGB(0.0, 0.0, 0.0))
    return f
    # save("/Users/s2257179/phd/stochastic_hybrid_code/$savein/$title.png", f)
end

df_results[1][:time]


WGLMakie.activate!()

f = plot_prop(df_results, df_props, 1, "thresh_plots2/props", "threshold_$(threshold_vals1[1])", threshold_vals1, 31856.296174439733)

f = create_subplots("plot_results", 5, 4);
[iscatter!(f[1,1],df_results[1].time, df_props[1][i], label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
axislegend(f[1,1], framevisible=true)
lines!(range(minimum(df_results[res_ind].time), maximum(df_results[res_ind].time), length=2), [threshold_vals[res_ind], threshold_vals[res_ind]], linewidth=4, color=RGB(0.0, 0.0, 0.0))

save("/Users/s2257179/phd/stochastic_hybrid_code/threshold_analysis/test.html", f)

for i in eachindex(df_props)
    plot_prop(df_results, df_props, i, "thresh_plots2/props", "threshold_$(threshold_vals1[i])", threshold_vals1, 31856.296174439733)
end

for i in eachindex(df_props0)
    plot_prop(df_results0, df_props0, i, "thresh_plots/props", "threshold_$(threshold_vals0[i])", threshold_vals0, 31856.296174439733)
end

for i in eachindex(df_props)
    plot_prop(df_results, df_props, i, "thresh_plots2/props_zoom", "threshold_$(threshold_vals1[i])", threshold_vals1, threshold_vals1[i]*2.5)
end

for i in eachindex(df_props0)
    plot_prop(df_results0, df_props0, i, "thresh_plots/props_zoom", "threshold_$(threshold_vals0[i])", threshold_vals0, threshold_vals0[i]*2.5)
end




plot_props(df_results, df_props, threshold_vals1)
