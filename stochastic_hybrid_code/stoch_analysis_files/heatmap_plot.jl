using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes, KernelDensity, Arrow, StatsBase

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/hist_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/data_processing_funcs.jl"))

kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

function produce_heatmap_data(simulation; threshold=2)
    states = []
    for i in simulation
        if i >= threshold
            push!(states, 1)
        else
            push!(states, 0)
        end
    end
    heatmap_data = reshape(states, 1, length(states))
    return heatmap_data
end

type_kdam = "low_kdam"
@load "/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_stops.jld2" df_lengths df_stops
df_rtca = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_times.arrow"))

res, times_res = remove_missing(df_rtca, df_times)

sims, times = split_into_simulations(res, times_res)

lines(times[1][0.08], sims[1][0.08])

sims[1][0.08][1:10:end]

heatmap_data = produce_heatmap_data(sims[1][0.02][1:50:end], threshold=2)
heatmap_data2 = produce_heatmap_data(sims[1][1.5][1:50:end], threshold=2)

# test = round.(rand(10000))
# heatmap_data = reshape(test, 1, length(test))
f = Figure()
ax = Axis(f[1,1])
img = image!(ax, heatmap_data';
    colormap = [:red, :green], 
    colorrange = (0, 1), alpha=0.5     # Set the range of values between 0 and 1
)
ax2 = Axis(f[1,2])
img = image!(ax2, heatmap_data2';
    colormap = [:red, :green], 
    colorrange = (0, 1), alpha=0.5      # Set the range of values between 0 and 1
)
hideydecorations!(ax)
hideydecorations!(ax2)
hidexdecorations!(ax, ticks=false, ticklabels=false)
hidexdecorations!(ax2, ticks=false, ticklabels=false)