
function determine_common_bin_edges(res_log, nbins)
    all_data = vcat(values(res_log)...)
    min_val = minimum(all_data)
    max_val = maximum(all_data)
    return range(min_val, max_val, length=nbins+1)
end

function prod_hist_data(res_log, time_data, logdensity::Bool)
    bin_edges = determine_common_bin_edges(res_log, 20)
    time_weights = diff(time_data)
    push!(time_weights, time_weights[end])
    hist_data = fit(Histogram, res_log, weights(time_weights), bin_edges)
    x = hist_data.edges[1][1:end-1] .+ diff(hist_data.edges[1]) ./ 2 #bar edges 
    if logdensity 
        y = log.(hist_data.weights.+1)
    else
        y = hist_data.weights
    end
    return x, y, hist_data
end

function produce_hist_data(res_log, time_data, logdensity::Bool)
    bin_edges = determine_common_bin_edges(res_log, 20)
    x_all = Dict(); y_all = Dict(); barlines = Dict();
    for kdam in kdams
        # println(kdam)
        x, y, hist_data = prod_hist_data(res_log[kdam], time_data[kdam], logdensity)
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



function hist3d_data(barlines_all, kdam)
    x = barlines_all[kdam][1]
    y = barlines_all[kdam][2]
    lower = Point3f.(x, kdam, 0)
    upper = Point3f.(x, kdam, y)
    return x, y, lower, upper 
end
function plot3d(barlines_all; second_dataset=[], title="", xlab="Molecules", ylab="Damage rate", zlab="Frequency", tosave=false)
    fig = Figure()
    ax = Axis3(fig[1, 1], title="$title", xlabel = "$xlab", ylabel = "$ylab", zlabel = "$zlab")
    c = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    for (i, kdam) in enumerate(kdams)
        x, y, lower, upper = hist3d_data(barlines_all, kdam)
        lines!(ax, x, fill(kdam, length(x)), y, linewidth = 2, color=c)
        band!(ax, lower, upper, transparency = true, color=c, alpha=0.2)
        if second_dataset != []
            x2, y2, lower2, upper2 = hist3d_data(second_dataset, kdam)
            lines!(ax, x2, fill(kdam, length(x2)), y2, linewidth = 2, color=c2)
            band!(ax, lower2, upper2, transparency = true, color=c2, alpha=0.5)
        end
    end
    if tosave
        mainpath = "/Users/s2257179/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Documents/rtc/stochastic/plots/analysis/histograms/"
        filename = if type_kdam == "high_kdam"
            "high_kdam_histogram.png"
        else
            "low_kdam_histogram.png"
        end
    
        # Check if the file already exists and modify the filename if it does
        filepath = joinpath(mainpath, filename)
        counter = 1
        while isfile(filepath)
            base, ext = splitext(filename)
            filepath = joinpath(mainpath, "$(base)_$(counter)$(ext)")
            counter += 1
        end
        save(filepath, fig)
    end
    return fig
end
