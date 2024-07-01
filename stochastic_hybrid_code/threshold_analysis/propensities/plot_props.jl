using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, GLMakie

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
    f = Figure(size=(1450, 800))
    ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", title=title, limits=(nothing,(0, max_val)))
    [iscatter!(ax,df_results[res_ind].time, df_props[res_ind][i], label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
    axislegend(ax, framevisible=true)
    lines!(range(minimum(df_results[res_ind].time), maximum(df_results[res_ind].time), length=2), [threshold_vals[res_ind], threshold_vals[res_ind]], linewidth=4, color=RGB(0.0, 0.0, 0.0))
    # return f
    save("/Users/s2257179/phd/stochastic_hybrid_code/$savein/$title.html", f)
end

df_results[1][:time]



f = plot_prop(df_results, df_props, 11, "thresh_plots2/props", "threshold_$(threshold_vals1[11])", threshold_vals1, 31856.296174439733)


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



WGLMakie.activate!()

f = plot_props(df_results, df_props, threshold_vals1)
save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots2/props/props.html", f)


f0 = plot_props(df_results0, df_props0, threshold_vals0)
save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/props/props.html", f0)



f = Figure(size=(1450, 800))
ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", limits=(nothing,(0, 31856.296174439733)))
[iscatter!(ax,df_results[11].time, df_props[11][i], label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
axislegend(ax, framevisible=true)
lines!(range(minimum(df_results[11].time), maximum(df_results[11].time), length=2), [threshold_vals[11], threshold_vals[11]], linewidth=4, color=RGB(0.0, 0.0, 0.0))

save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots2/props/props_test.html", f)





WGLMakie.activate!()


f = Figure()
ax = Axis(f[1,1])
lines!(ax, 1:4)

save("/Users/s2257179/Desktop/test.html", f)


using GLMakie
GLMakie.activate!()

f = Figure(size=(1450, 800))
ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", limits=(nothing,(0, 31856.296174439733)))
[iscatter!(ax,df_results[11].time, df_props[11][i], label="$(react_names[i])", color=color_list[i]) for i in eachindex(react_names[1:end-1])]
axislegend(ax, framevisible=true)
lines!(range(minimum(df_results[11].time), maximum(df_results[11].time), length=2), [threshold_vals[11], threshold_vals[11]], linewidth=4, color=RGB(0.0, 0.0, 0.0))



min_time = minimum(df_results[11].time)
max_time = maximum(df_results[11].time)
threshold_val = [threshold_vals1[11], threshold_vals1[11]]

f = Figure(size=(1450, 800))
ax = Axis(f[1,1], xlabel="time", ylabel="propensity", limits=(nothing, (0, 31856.296174439733)))

for i in eachindex(react_names[1:end-1])
    iscatter!(ax, df_results[11].time, df_props[11][i], label="$(react_names[i])", color=color_list[i])
end

ax = Axis(f[1,2], xlabel="time", ylabel="propensity", limits=(nothing, (0, 31856.296174439733)))

for i in eachindex(react_names[1:end-1])
    iscatter!(ax, df_results[1].time, df_props[1][i], label="$(react_names[i])", color=color_list[i])
end

axislegend(ax, framevisible=true)

# Use precomputed values
lines!(range(min_time, max_time, length=2), threshold_val, linewidth=4, color=RGB(0.0, 0.0, 0.0))

cols=2
rows=2
f = Figure(size=(1450, 800))
for j in 1:cols
    for i in 1:rows
        data_ind = i + rows * (j - 1)
        ax = Axis(f[i,j], xlabel = "time", ylabel = "propensity", limits=(nothing,(0, 31856)))
        for i in eachindex(react_names[1:end-1])
            iscatter!(ax, df_results[data_ind].time, df_props[1][i], label="$(react_names[i])", color=color_list[i])
        end
        lines!(range(minimum(df_results[data_ind].time), maximum(df_results[data_ind].time), length=2), [threshold_vals1[data_ind], threshold_vals1[data_ind]], linewidth=4, color=RGB(0.0, 0.0, 0.0))

        # Hide x-axis decorations for axes not in the bottom row
        if i != rows
            hidexdecorations!(ax, grid=false)
        end

        # Hide y-axis decorations for axes not in the first column
        if j > 1
            hideydecorations!(ax, grid=false)
        end
    end
end
axislegend(ax, framevisible=true)

linkaxes!(filter(x -> x isa Axis, f.content)...)










open("index.html", "w") do io
    println(io, """
    <html>
        <head>
        </head>
        <body>
    """)
    Page(exportable=true, offline=true)
    # Then, you can just inline plots or whatever you want :)
    # Of course it would make more sense to put this into a single app
    app = App() do
        f = Figure()
        ax = Axis(f[1,1])
        lines!(ax, 1:4)
    end
    show(io, MIME"text/html"(), app)
    # or anything else from Bonito, or that can be displayed as html:
    println(io, """
        </body>
    </html>
    """)
end