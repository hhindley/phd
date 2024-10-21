using JLD2, InteractiveViz, Statistics, DataFrames, KernelDensity, Arrow, StatsBase

using GLMakie
fontsize_theme = Theme(fontsize = 10)
set_theme!(fontsize_theme)

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/hist_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/data_processing_funcs.jl"))


kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

server = false
if server
    mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/"
else
    mainpath = "/Users/s2257179/Desktop/saved_variables/"
end

type_kdam = "low_kdam"
@load "$mainpath/$type_kdam/hist_data.jld2" x_all_low y_all_low barlines_all_low x_on_low y_on_low barlines_on_low x_off_low y_off_low barlines_off_low
df_res_low = DataFrame(Arrow.Table("$mainpath/$type_kdam/all_results_grouped.arrow"))

type_kdam = "high_kdam"
@load "$mainpath/$type_kdam/hist_data.jld2" x_all_high y_all_high barlines_all_high x_on_high y_on_high barlines_on_high x_off_high y_off_high barlines_off_high
df_res_high = DataFrame(Arrow.Table("$mainpath/$type_kdam/all_results_grouped.arrow"))

# violin plot 
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min⁻¹)", ylabel="RtcA in on state (μM)", title="Hysteresis experiement")
violin!(ax, df_res_high.group, df_res_high.rtca, side=:left)


f = Figure()
ax = Axis(f[1,1], xlabel="Log Molecules", ylabel = "Log frequency")
barplot!(ax, x_all_high[0.08], y_all_high[0.08], gap=0)
# lines!(ax, x_all[0.08], y_all[0.08], color = :red)
lines!(ax, barlines_all_high[0.08][1], barlines_all_high[0.08][2], color=:red)
display(GLMakie.Screen(), f)

# work out how to do loglog kde calculation 
# kde_all = Dict{Float64, KernelDensity.UnivariateKDE}()
# for kdam in kdams
#     kde_all[kdam] = kde(res_log[kdam])
# end

# lines(kde_all[0.08].x, kde_all[0.08].density)


plot3d(barlines_all_high, xlab="Log molecules", zlab="Log frequency", title="whole dataset - high", tosave=true)
plot3d(barlines_on_high, xlab="Log molecules", title="on state frequencies - high", tosave=true)
plot3d(barlines_off_high, title="off state frequencies - high", tosave=true)

plot3d(barlines_all_low, xlab="Log molecules", zlab="Log frequency", title="whole dataset - low", tosave=true)
plot3d(barlines_on_low, xlab="Log molecules", title="on state frequencies - low", tosave=true)
plot3d(barlines_off_low, title="off state frequencies - low", tosave=true)

plot3d(barlines_all_high, second_dataset=barlines_all_low, xlab="Log molecules", zlab="Log frequency", title="whole dataset", tosave=true)
plot3d(barlines_on_high, second_dataset=barlines_on_low, xlab="Log molecules", title="on state frequencies", tosave=true)
plot3d(barlines_off_high, second_dataset=barlines_off_low, title="off state frequencies", tosave=true)





# work out time in each bin (not sure it was worth it to do this)

include("/Users/s2257179/phd/stochastic_hybrid_code/stoch_analysis_files/histograms/make_hists.jl")

bin_edges = determine_common_bin_edges(res_log, 20)

function segment_simulations(df)
    # Identify the start of new simulations
    df = dropmissing(df, :t)
    new_sim_starts = findall(diff([0; df.t]) .< 0)
    new_sim_starts = [1; new_sim_starts; nrow(df) + 1]  # Add the first and last indices

    # Segment the DataFrame by each simulation
    segments = [df[new_sim_starts[i]:(new_sim_starts[i+1]-1), :] for i in 1:(length(new_sim_starts)-1)]
    return segments
end

function total_time_in_state_multiple_simulations(df)
    segments = segment_simulations(df)
    total_time = 0.0
    for segment in segments
        total_time += vec_time_in_state(segment)
    end
    return total_time
end


function group_dataframe(res_log, kdam)
    df = DataFrame(s=res_log[kdam],t=times_res["$kdam"])

    df.bin = cut(df.s, bin_edges, extend=true)

    df.bin_index = levelcode.(df.bin)

    df.actual_index = 1:length(df.s)

    all_bins_df = DataFrame(bin_index = 1:20)

    full_df = leftjoin(all_bins_df, df, on = :bin_index)

    grouped_df = groupby(full_df, :bin_index)
    return grouped_df
end

function calc_time_in_bin(res_log)
    state_times = Dict(kdam=>Dict() for kdam in kdams)
    for kdam in kdams
        grouped_df = group_dataframe(res_log, kdam)
        for i in 1:length(grouped_df)
            # println(i)
            df = grouped_df[i]
            tot_time_in_state = total_time_in_state_multiple_simulations(df)
            state_times[kdam]["bin$i"] = tot_time_in_state
        end
    end
    return state_times
end

state_times_log = calc_time_in_bin(res_log)

state_times_log[0.0]    
for kdam in kdams
    println(sum([state_times_log[kdam]["bin$i"] for i in 1:length(state_times)]))
end

norm_states = Dict(kdam=>[] for kdam in kdams)
for kdam in kdams
    max_t = maximum(collect(values(state_times_log[kdam])))
    min_t = minimum(collect(values(state_times_log[kdam])))
    # norm_states[kdam] = ((collect(values(state_times_log[kdam]))) .- min_t) ./(max_t-min_t)
    norm_states[kdam] = ([state_times_log[0.0]["bin$i"] for i in 1:length(state_times_log[0.0])] .- min_t) ./(max_t-min_t)

end
norm_states
[state_times_log[0.0]["bin$i"] for i in 1:length(state_times_log[0.0])]
maximum(collect(values(state_times_log[0.0])))
minimum(collect(values(state_times_log[0.0])))

norm_states[0.0]
lines(kdams,norm_states[0.0])
t=lines(x,y)
display(GLMakie.Screen(), t)

color_vals = Dict()
for kdam in kdams
    color_vals[kdam] = [get(ColorSchemes.rainbow_bgyr_35_85_c72_n256, log10(frac+1)) for frac in collect(values(norm_states[kdam]))]
end
color_vals[0.0]
color_vals[0.08]

lines(barlines_on[0.08][1], barlines_on[0.08][2], color = color_vals[0.08])


x = barlines_all[0.08][1]
y = barlines_all[0.08][2]    
color_vals[0.08]
# Create points for segments, each pair of points defines a segment
points = Point2f[(x[i], y[i]) for i in 1:length(x)-1 for j in [i, i+1]]

# Plot with LineSegments and assign colors
f = Figure()
ax = Axis(f[1, 1])
linesegments!(ax, points, color=color_vals[0.08], linewidth=2)
