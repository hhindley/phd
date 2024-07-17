all_species = [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rd, :rt]

function plotBIG(x, y; xtitle="", ytitle="", title="")
    xvals = range(minimum(x), maximum(x), length=length(x))
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "$xtitle", ylabel = "$ytitle", title="$title")
    ilines!(f[1,1], xvals, y)

    return f
end
function plot_times(df_times, title; folder="")
    threshold_vals = df_times.threshold
    title = "$title, thresh vals: $(minimum(threshold_vals)), $(maximum(threshold_vals)), step=$(threshold_vals[2]-threshold_vals[1]), length=$(length(threshold_vals))"
    f=Figure()
    ax=Axis(f[1,1],xlabel="threshold", ylabel="time (hours)", title=title)
    lines!(ax,df_times.threshold, df_times.time/60/60)
    if folder != ""
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        if !isdir(plots_folder)
            mkdir(plots_folder)
        end
        save(joinpath(plots_folder, "times.png"), f)
    end
    return f 
end
function plot_totstochcount(threshold_vals, tot_counts, title; folder="")
    title = "$title, thresh vals: $(minimum(threshold_vals)), $(maximum(threshold_vals)), step=$(threshold_vals[2]-threshold_vals[1]), length=$(length(threshold_vals))"
    f = Figure()
    ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count", title=title)
    barplot!(threshold_vals, tot_counts)
    if folder != ""
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        if !isdir(plots_folder)
            mkdir(plots_folder)
        end
        save(joinpath(plots_folder, "tot_stoch_count.png"), f)
    end
    return f 
end




struct PlotWrapper
    plot::Any
end
function get_column(df::DataFrame, column_name::Union{Symbol, String})
    return getproperty(df, column_name)
end

function create_subplots(plotting_func, num_plots; size=(600, 450), xlabel="", ylabel="", titles=[], hidelabels=[true, true], yscale=identity)
    f = Figure(size=size)
    base = ceil(sqrt(num_plots))
    columns = Int(base)
    rows = Int(ceil(num_plots / columns))
    
    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            # if length(titles) % 2 == 1
            #     titles = push!(titles, "empty")
            # end
            if data_ind <= num_plots
                title = titles != [] ? titles[data_ind] : "Plot $data_ind"
                if plotting_func == "plot_results" || plotting_func == "plot_hists" || plotting_func == "plot_individual_reacts"
                    ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale, xticklabelrotation=45)
                elseif plotting_func == "plot_stoch_reacts"
                    ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale, xticks=(1:13, react_names_str), xticklabelrotation=45)
    
                end
                # Hide x-axis decorations for axes not in the bottom row
                if i != rows && (hidelabels == [true, true] || hidelabels == [true, false])
                    hidexdecorations!(ax, grid=false)
                end

                # Hide y-axis decorations for axes not in the first column
                if j > 1 && (hidelabels == [true, true] || hidelabels == [false, true])
                    hideydecorations!(ax, grid=false)
                    
                end
            end
        end
    end
    return f
end
function add_subplots(f, plotting_func, df_results, num_plots; species=:rm_a, linkaxes=true)
    if plotting_func == "plot_results"
        preprocessed_times = [get_column(df, :time) for df in df_results]
        preprocessed_results = [get_column(df, species) for df in df_results]
    elseif plotting_func == "plot_hists"
        dfs = df_results[string(species)]
    end
    base = ceil(sqrt(num_plots))
    columns = Int(base)
    rows = Int(ceil(num_plots / columns))
    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            if data_ind <= num_plots
                if plotting_func == "plot_results"
                    plot_timeres(preprocessed_times[data_ind], preprocessed_results[data_ind], f[i,j])
                elseif plotting_func == "plot_hists"
                    plot_hist(dfs[data_ind], f[i,j])
                elseif plotting_func == "plot_stoch_reacts" 
                    barplot!(f[i,j], df_results[data_ind].event, df_results[data_ind].count)
                elseif plotting_func == "plot_individual_reacts"
                    if data_ind <= length(react_names[1:end-1])
                        barplot!(f[i,j], df_results[react_names[1:end-1][data_ind]].threshold, df_results[react_names[1:end-1][data_ind]].count)
                    end
                elseif plotting_func == "plot_props"
                    time_data = df_results[data_ind].time
                    prop_data = [df_props[data_ind][react] for react in eachindex(react_names[1:end-1])]

                    for r in eachindex(react_names[1:end-1])
                        iscatter!(f[i,j], time_data, prop_data[r], label="$(react_names[r])", color=color_list[r])
                    end    
                    axislegend(f[i,j], framevisible=true)
                    lines!(f[i,j], range(minimum(time_data), maximum(time_data), length = 2), [threshold_vals[data_ind], threshold_vals[data_ind]], linewidth = 4, color = :black)
                end
            end
        end
    end
    if linkaxes == true
        linkaxes!(filter(x -> x isa Axis, f.content)...)
    end
    return f
end
function plot_results(plotting_func, df_results, num_plots; species=:rm_a, size=(600, 450), xlabel="", ylabel="", titles=[], yscale=identity, hidelabels=[true,true], linkaxes=true, folder="")
    f = create_subplots(plotting_func, num_plots; size=size, xlabel=xlabel, ylabel=ylabel, titles=titles, hidelabels=hidelabels, yscale=yscale)
    f = add_subplots(f, plotting_func, df_results, num_plots; species=species, linkaxes=linkaxes)
    
    if !isempty(folder)
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        if plotting_func == "plot_stoch_reacts"
            save(joinpath(plots_folder, "stoch_reacts.png"), f)
        else
            results_folder = joinpath(plots_folder, plotting_func)
            if !isdir(results_folder)
                mkdir(results_folder)
            end
            save(joinpath(results_folder, "$species.png"), f)
        end
            
        # return PlotWrapper(f)
    else
        return f
    end
end


function plot_timeres(time, res, loc)
    x = makiex(time)
    ilines!(loc, x, res)
end
function plot_hist(df, loc)
    bins = [i[1] for i in df.bin]
    push!(bins, df.bin[end][2])
    bin_c = (bins[1:end-1] .+ bins[2:end]) ./ 2
    barplot!(loc, bin_c, df.freq, width=diff(bins), gap=0)
    # return f
end
function makiex(x)
    return range(minimum(x), maximum(x), length=length(x))
end

function build_reaction_count_df(df_reacts, react, threshold_vals)
    df_tag=[]
    for i in eachindex(df_reacts)
        df = filter(row -> row.reaction == react, df_reacts[i])
        df.threshold = [threshold_vals[i]]
        push!(df_tag, df)
    end
    combined_df = reduce(vcat, df_tag)
    return combined_df
end


function plot_prop(df_results, df_props, res_ind, title, threshold_vals; maxval=31856, folder="")
    time_data = df_results[res_ind].time
    prop_data = [df_props[res_ind][i] for i in eachindex(react_names[1:end-1])]
    
    f = Figure()
    ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", title=title, limits=((2.2e5, 2.6e5),(0, maxval)))
    for i in eachindex(react_names[1:end-1])
        iscatter!(ax, time_data, prop_data[i], label="$(react_names[i])", color=color_list[i])
    end    
    axislegend(ax, framevisible=true)
    lines!(ax, range(minimum(time_data), maximum(time_data), length = 2), [threshold_vals[res_ind], threshold_vals[res_ind]], linewidth = 4, color = :black)
    
    if !isempty(folder)
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        results_folder = joinpath(plots_folder, "propensities")
        if !isdir(results_folder)
            mkdir(results_folder)
        end
        
        save(joinpath(results_folder, "$title.png"), f)
        # return PlotWrapper(f)
    else
        return f
    end
end

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
];

react_names_str = [string(i) for i in react_names]
