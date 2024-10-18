using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes, KernelDensity, Arrow, StatsBase

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

type_kdam = "high_kdam"
@load "/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_stops.jld2" df_lengths df_stops 

df_rtca = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_times.arrow"))

kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

res = Dict()
times_res = Dict()
for i in eachindex(names(df_rtca))
    res[names(df_rtca)[i]] = collect(skipmissing(df_rtca[!,names(df_rtca)[i]]))
    times_res[names(df_times)[i]] = collect(skipmissing(df_times[!,names(df_times)[i]]))
end
res
res_log = Dict()
for kdam in kdams
    res_log[kdam] = log.(res["$kdam"].+1)
end

res_on = Dict([kdam=>Float64[] for kdam in kdams])
res_on_log = Dict([kdam=>Float64[] for kdam in kdams])
for kdam in kdams
    res_on[kdam] = [i for i in res["$kdam"] if i >= 2]
    res_on_log[kdam] = [log.(i.+1) for i in res["$kdam"] if i>=2]
end

res_off = Dict([kdam=>Float64[] for kdam in kdams])
res_off_log = Dict([kdam=>Float64[] for kdam in kdams])
for kdam in kdams
    res_off[kdam] = [i for i in res["$kdam"] if i < 2]
    res_off_log[kdam] = [log.(i.+1) for i in res["$kdam"] if i < 2]
end

all_res=[]
for i in kdams
    push!(all_res, res["$i"])
end
tot_lengths = [sum(df_lengths[!,"$kdam"]) for kdam in kdams]

groups = [fill(i, tot_lengths[i]) for i in eachindex(kdams)]

df_res = DataFrame(rtca=vcat(all_res...), group=vcat(groups...))

f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min⁻¹)", ylabel="RtcA in on state (μM)", title="Hysteresis experiement")
violin!(ax, df_res.group, df_res.rtca, side=:left)


# split into original simulations
sims = Dict(i=>Dict(kdam=>Float64[] for kdam in kdams) for i in 1:20)
times = Dict(i=>Dict(kdam=>Float64[] for kdam in kdams) for i in 1:20)
for kdam in eachindex(kdams)
    sims[1][kdams[kdam]] = res["$(kdams[kdam])"][1:df_stops[1,"$(kdams[kdam])"]]
    times[1][kdams[kdam]] = times_res["$(kdams[kdam])"][1:df_stops[1,"$(kdams[kdam])"]]
    for j in eachindex(kdams)
        if j < 20
            ind = j+1
        else 
            continue
        end
        sims[ind][kdams[kdam]] = res["$(kdams[kdam])"][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
        times[ind][kdams[kdam]] = times_res["$(kdams[kdam])"][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
    end
end

lines(times[1][0.08], sims[1][0.08])

on_times = [getindex(times[1][0.0][i]) for i in eachindex(sims[1][0.0]) if sims[1][0.0][i] >= 2]
off_times = [getindex(times[1][0.0][i]) for i in eachindex(sims[1][0.0]) if sims[1][0.0][i] < 2]

states = []
for i in sims[1][0.0]
    if i >= 2
        push!(states, 1)
    else
        push!(states, 0)
    end
end
heatmap_data = reshape(states, 1, length(states))

test = round.(rand(50000))
heatmap_data = reshape(test, 1, length(test))
f = Figure()
img = image(f[1,1], heatmap_data';
    colormap = [:red, :green],  # Blue for off, Red for on
    colorrange = (0, 1),       # Set the range of values between 0 and 1
)


function determine_common_bin_edges(res_log, nbins)
    all_data = vcat(values(res_log)...)
    min_val = minimum(all_data)
    max_val = maximum(all_data)
    return range(min_val, max_val, length=nbins+1)
end

function prod_hist_data(res_log, logdensity::Bool)
    bin_edges = determine_common_bin_edges(res_log, 20)
    hist_data = fit(Histogram, res_log, bin_edges)
    x = hist_data.edges[1][1:end-1] .+ diff(hist_data.edges[1]) ./ 2 #bar edges 
    if logdensity 
        y = log.(hist_data.weights.+1)
    else
        y = hist_data.weights
    end
    return x, y, hist_data
end

function produce_hist_data(res_log, logdensity::Bool)
    bin_edges = determine_common_bin_edges(res_log, 20)
    x_all = Dict(); y_all = Dict(); barlines = Dict();
    for kdam in kdams
        # println(kdam)
        x, y, hist_data = prod_hist_data(res_log[kdam], logdensity)
        barline = create_barline(hist_data, y, bin_edges)
        x_all[kdam] = x
        y_all[kdam] = y
        barlines[kdam] = barline
    end
    return x_all, y_all, barlines
end

function create_barline(hist_data, y, bin_edges)
    barline=Float64[]
    ys = Float64[]
    for i in 1:length(hist_data.edges[1][1:end-1])
        push!(barline, hist_data.edges[1][1:end-1][i])
        push!(barline, hist_data.edges[1][1:end-1][i] .+ diff(hist_data.edges[1])[i])
        push!(ys, y[i])
        push!(ys, y[i])
    end
    if maximum(barline) < maximum(bin_edges)
        max_barline = maximum(barline)
        push!(barline, max_barline)
        push!(ys, 0)
        push!(barline, maximum(bin_edges))
        push!(ys, 0)
    end
    return (barline, ys)
end


# all data 
x_all, y_all, barlines_all = produce_hist_data(res_log, true)

f = Figure()
ax = Axis(f[1,1], xlabel="Log Molecules", ylabel = "Log frequency")
barplot!(ax, x_all[0.08], y_all[0.08], gap=0)
# lines!(ax, x_all[0.08], y_all[0.08], color = :red)
lines!(ax, barlines_all[0.08][1], barlines_all[0.08][2], color=:red)
display(GLMakie.Screen(), f)

# work out how to do loglog kde calculation 
kde_all = Dict{Float64, KernelDensity.UnivariateKDE}()
for kdam in kdams
    kde_all[kdam] = kde(res_log[kdam])
end

lines(kde_all[0.08].x, kde_all[0.08].density)


# on data 
x_on, y_on, barlines_on = produce_hist_data(res_on_log, false)

f = Figure()
ax = Axis(f[1,1], xlabel="Log Molecules", ylabel = "Frequency", title="on state frequencies")
barplot!(ax, x_on[0.08], y_on[0.08], gap=0)
lines!(ax, x_on[0.08], y_on[0.08], color = :red)
display(GLMakie.Screen(), f)

kde_on = Dict{Float64, KernelDensity.UnivariateKDE}()
for kdam in kdams
    kde_on[kdam] = kde(res_on_log[kdam])
end

lines(kde_on[0.08].x, kde_on[0.08].density)

# off data - not really relevant using a threshold of 2 as most of these are just 0?
x_off, y_off, barlines_off = produce_hist_data(res_off, false)

f = Figure()
ax = Axis(f[1,1], xlabel="Log Molecules", ylabel = "Frequency", title="off state frequencies")
barplot!(ax, x_off[0.08], y_off, gap=0)
lines!(ax, x_off[0.08], y_off[0.08], color = :red)
display(GLMakie.Screen(), f)

kde_off = Dict{Float64, KernelDensity.UnivariateKDE}()
for kdam in kdams
    kde_off[kdam] = kde(res_off_log[kdam])
end

lines(kde_on[0.08].x, kde_on[0.08].density)




# function plot3d(x, y, z; zmin = minimum(z), lw = 1., colmap = :linear_bgy_10_95_c74_n256, colorband = (:lightblue, 1.), xlab = "x", ylab = "y", zlab = "z")
#     # Initialisation
#     fig = Figure()
#     ax = Axis3(fig[1,1], xlabel = xlab, ylabel = ylab, zlabel = zlab)
#     for (j, yv) in enumerate(y)
        
#         zj = z[j, :]
#         println(length(zj))
#         println(length(yv))
#         lower = Point3f.(x, yv, zmin)
#         upper = Point3f.(x, yv, zj)
#         # edge_start = [Point3f(x[1], yv, zmin), Point3f(x[1], yv, zj[1])]
#         # edge_end = [Point3f(x[end], yv, zmin), Point3f(x[end], yv, zj[end])]

#         # # Surface
#         band!(ax, lower, upper, color = colorband)

#         # Line
#         lines!(ax, upper, color = zj)#, colormap = colmap, linewidth = lw)
        
#         # Edges
#         # lines!(ax, edge_start, color = zj[1]*ones(2), colormap = colmap, linewidth = lw)
#         # lines!(ax, edge_end, color = zj[end]*ones(2), colormap = colmap, linewidth = lw)
#     end

#     # Set axes limits
#     # xlims!(ax, minimum(x), maximum(x))
#     # ylims!(ax, minimum(y), maximum(y))
#     # zlims!(ax, zmin, maximum(z))

#     fig
# end


# kdes_res = [kde(res[i]) for i in eachindex(res)]
# kdes_res = [kde(log.(res[i] .+ 1)) for i in eachindex(res)]  # Apply log transformation

# x = kdes_res[1].x
# y=1:length(res)
# z = hcat([kde.density for kde in kdes_res]...)
# z = permutedims(z, (2, 1))
# nx = length(x)
# ny = length(y)


# x = kde_all[0.0].x
# length(x)
# y = 1:length(res)
# length(y)
# z = hcat([kde_all[kdam].density for kdam in kdams]...)
# z = permutedims(z, (2, 1))

# fig = plot3d(x, y, z, xlab="conc", ylab="dam", zlab="freq")


# td = lines(barlines_all[0.0][1], barlines_all[0.0][2])
# display(GLMakie.Screen(), td)   

# x = range(minimum(barlines_all[0.08][1]),maximum(barlines_all[0.08][1]), length=40)
# length(x)
# y = 1:length(res)
# barlines_all
# z = hcat([barlines_all[kdam][2] for kdam in kdams]...)
# z = permutedims(z, (2, 1))

# fig = plot3d(x, y, z, xlab="conc", ylab="dam", zlab="freq")
# display(GLMakie.Screen(), fig)


Makie.wong_colors()[1]


function hist3d_data(barlines_all, kdam)
    x = barlines_all[kdam][1]
    y = barlines_all[kdam][2]
    lower = Point3f.(x, kdam, 0)
    upper = Point3f.(x, kdam, y)
    return x, y, lower, upper 
end
function plot3d(barlines_all; title="", xlab="Molecules", ylab="Damage rate", zlab="Frequency", tosave=false)
    fig = Figure()
    ax = Axis3(fig[1, 1], title="$title", xlabel = "$xlab", ylabel = "$ylab", zlabel = "$zlab")
    c = Makie.wong_colors()[1]
    for (i, kdam) in enumerate(kdams)
        x, y, lower, upper = hist3d_data(barlines_all, kdam)
        lines!(ax, x, fill(kdam, length(x)), y, linewidth = 2, color=c)
        band!(ax, lower, upper, transparency = true, color=c, alpha=0.2)
    end
    if tosave
        mainpath = "/Users/s2257179/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Documents/rtc/stochastic/plots/analysis/histograms/"
        if type_kdam == "high_kdam"
            save(joinpath(mainpath, "high_kdam_histogram.png"), fig)
        else
            folder = "init_vals_low_kdam/thresh$thresh"
        end
    end
    return fig
end

plot3d(barlines_all, xlab="Log molecules", zlab="Log frequency", title="whole dataset")

plot3d(barlines_on, xlab="Log molecules", title="on state frequencies")



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
