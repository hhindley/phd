using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes, KernelDensity, Arrow, StatsBase

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

df_rtca = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_times.arrow"))

res, times_res = remove_missing(df_rtca, df_times)

sims, times = split_into_simulations(res, times_res)

lines(times[1][0.08], sims[1][0.08])



heatmap_data = produce_heatmap_data(sims[1][0.08], threshold=2)

test = round.(rand(50000))
heatmap_data = reshape(test, 1, length(test))
f = Figure()
img = image(f[1,1], heatmap_data';
    colormap = [:red, :green],  # Blue for off, Red for on
    colorrange = (0, 1),       # Set the range of values between 0 and 1
)