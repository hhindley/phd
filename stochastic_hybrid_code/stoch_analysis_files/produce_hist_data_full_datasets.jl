using JLD2, InteractiveViz, Statistics, DataFrames, KernelDensity, Arrow, StatsBase

using GLMakie
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/hist_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/data_processing_funcs.jl"))


kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

server = true
if server
    mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/"
else
    mainpath = "/Users/s2257179/Desktop/saved_variables/"
end

type_kdam = "low_kdam"
@load "$mainpath/$type_kdam/$(type_kdam)_stops.jld2" df_lengths df_stops 
df_lengths_low = df_lengths; df_stops_low = df_stops
df_rtca_low = DataFrame(Arrow.Table("$mainpath/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times_low = DataFrame(Arrow.Table("$mainpath/$type_kdam/$(type_kdam)_times.arrow"))

type_kdam = "high_kdam"
@load "$mainpath/$type_kdam/$(type_kdam)_stops.jld2" df_lengths df_stops 
df_lengths_high = df_lengths; df_stops_high = df_stops
df_rtca_high = DataFrame(Arrow.Table("$mainpath/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times_high = DataFrame(Arrow.Table("$mainpath/$type_kdam/$(type_kdam)_times.arrow"))

res_high, times_res_high = remove_missing(df_rtca_high, df_times_high)
res_log_high = log_results(res_high)
res_on_high, res_off_high = determine_state(res_high, threshold=2)
res_on_log_high = log_results(res_on_high)
res_off_log_high = log_results(res_off_high)

res_low, times_res_low = remove_missing(df_rtca_low, df_times_low)
res_log_low = log_results(res_low)
res_on_low, res_off_low = determine_state(res_low, threshold=2)
res_on_log_low = log_results(res_on_low)
res_off_log_low = log_results(res_off_low)

# all data 
x_all_high, y_all_high, barlines_all_high = produce_hist_data(res_log_high, true)
x_all_low, y_all_low, barlines_all_low = produce_hist_data(res_log_low, true)

# on data 
x_on_high, y_on_high, barlines_on_high = produce_hist_data(res_on_log_high, false)
x_on_low, y_on_low, barlines_on_low = produce_hist_data(res_on_log_low, false)

# off data - not really relevant using a threshold of 2 as most of these are just 0?
x_off_high, y_off_high, barlines_off_high = produce_hist_data(res_off_high, false)
x_off_low, y_off_low, barlines_off_low = produce_hist_data(res_off_low, false)

@save "$mainpath/high_kdam/hist_data.jld2" x_all_high y_all_high barlines_all_high x_on_high y_on_high barlines_on_high x_off_high y_off_high barlines_off_high
@save "$mainpath/low_kdam/hist_data.jld2" x_all_low y_all_low barlines_all_low x_on_low y_on_low barlines_on_low x_off_low y_off_low barlines_off_low

"scp hollie_hindley@resistance.bio.ed.ac.uk:/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/high_kdam/hist_data.jld2 /Users/s2257179/Desktop/saved_variables/high_kdam/hist_data.jld2"
"scp hollie_hindley@resistance.bio.ed.ac.uk:/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/low_kdam/hist_data.jld2 /Users/s2257179/Desktop/saved_variables/low_kdam/hist_data.jld2"




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
