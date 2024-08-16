all_species = [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rd, :rt]

function plotBIG(x, y; xtitle="", ytitle="", title="")
    xvals = range(minimum(x), maximum(x), length=length(x))
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "$xtitle", ylabel = "$ytitle", title="$title")
    ilines!(f[1,1], xvals, y)

    return f
end
function plot_times(df_times, title; folder="")
    if "thresh" in names(df_times)
        threshold_vals = df_times.threshold
        title = "$title, \n thresh vals: $(minimum(threshold_vals)), $(maximum(threshold_vals)), step=$(threshold_vals[2]-threshold_vals[1]), length=$(length(threshold_vals))"
        f=Figure()
        ax=Axis(f[1,1],xlabel="threshold", ylabel="time (hours)", title=title)
        lines!(ax,df_times.threshold, df_times.time/60/60)
    elseif "kdam" in names(df_times) && "thresh" in names(df_times)
        threshold_vals = df_times.thresh
        kdam_vals = df_times.kdam
        title = "$title, \n kdam vals: $kdam_vals \n thresh vals: $threshold_vals"
        f=Figure()
        ax=Axis(f[1,1],xlabel="run", ylabel="time (hours)", title=title)
        lines!(ax,1:length(df_times.kdam), df_times.time/60/60)
    elseif "kdam" in names(df_times)
        threshold_vals = df_times.kdam
        title = "$title, \n kdam vals: $(minimum(threshold_vals)), $(maximum(threshold_vals)), step=$(threshold_vals[2]-threshold_vals[1]), length=$(length(threshold_vals))"
        f=Figure()
        ax=Axis(f[1,1],xlabel="kdam", ylabel="time (hours)", title=title)
        lines!(ax,df_times.kdam, df_times.time/60/60)
    end
    
    
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
    title = "$title, \n thresh vals: $(minimum(threshold_vals)), $(maximum(threshold_vals)), step=$(threshold_vals[2]-threshold_vals[1]), length=$(length(threshold_vals))"
    f = Figure()
    ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count", title=title, xticks=(1:length(threshold_vals), [string(i) for i in threshold_vals]), xticklabelrotation=45)
    barplot!(f[1,1], eachindex(threshold_vals), tot_counts)
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

function create_subplots(plotting_func, num_plots, folder; size=(600, 450), xlabel="", ylabel="", titles=[], hidelabels=[true, true], yscale=identity)
    f = Figure(size=size)
    base = ceil(sqrt(num_plots))
    columns = Int(base)
    rows = Int(ceil(num_plots / columns))

    f[0, :] = Label(f, "$folder", fontsize=14)
    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            # if length(titles) % 2 == 1
            #     titles = push!(titles, "empty")
            # end
            if data_ind <= num_plots
                if num_plots > 1
                    title = titles != [] ? titles[data_ind] : "Plot $data_ind"
                else
                    title = titles != [] ? titles[1] : "Plot $data_ind"
                end
                if plotting_func == "plot_results" || plotting_func == "plot_hists" || plotting_func == "plot_individual_reacts"
                    ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale, xticklabelrotation=45, xlabelsize=10, ylabelsize=10, xticklabelsize=10, yticklabelsize=10, titlesize=12)
                elseif plotting_func == "plot_stoch_reacts"
                    ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale, xticks=(1:13, react_names_str), xticklabelrotation=45, xlabelsize=10, ylabelsize=10, xticklabelsize=10, yticklabelsize=10, titlesize=12)
    
                end
                colsize!(f.layout, j, Relative(1/columns))
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
        if num_plots > 1
            preprocessed_times = [get_column(df, :time) for df in df_results]
            preprocessed_results = [get_column(df, species) for df in df_results]
        else
            preprocessed_times = [get_column(df_results, :time)]
            preprocessed_results = [get_column(df_results, species)]
        end
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
                    plot_hist(dfs[data_ind], loc=f[i,j])
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
    if linkaxes == true && num_plots > 1
        linkaxes!(filter(x -> x isa Axis, f.content)...)
    end
    return f
end
function plot_results(plotting_func, df_results, num_plots, folder; species=:rm_a, size=(600, 450), xlabel="", ylabel="", titles=[], yscale=identity, hidelabels=[true,true], linkaxes=true, tosave=false)
    f = create_subplots(plotting_func, num_plots, folder; size=size, xlabel=xlabel, ylabel=ylabel, titles=titles, hidelabels=hidelabels, yscale=yscale)
    f = add_subplots(f, plotting_func, df_results, num_plots; species=species, linkaxes=linkaxes)
    
    if tosave
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
    end

    return f
end


function plot_timeres(time, res, loc)
    x = makiex(time)
    ilines!(loc, x, res)
end
function plot_hist(df; loc=nothing, label="", specie="", maxval=8e4)
    bins = [i[1] for i in df.bin]
    push!(bins, df.bin[end][2])
    bin_c = (bins[1:end-1] .+ bins[2:end]) ./ 2
    
    if isnothing(loc)
        f = Figure()
        ax = Axis(f[1,1], xlabel="bin", ylabel="frequency", title="$specie", limits=((nothing, nothing), (0, maxval)))
        barplot!(ax, bin_c, df.freq, width=diff(bins), gap=0, label=label)
        return f, ax
    else
        barplot!(loc, bin_c, df.freq, width=diff(bins), gap=0, label=label)
        ax = loc
        return ax
    end
end
function hists_with_less_bins(df, num_groups)
    df.group_id = ceil.(Int, (1:length(df.bin)) / (length(df.bin) / num_groups))
    agg_df = combine(groupby(df, :group_id), :freq => sum, :bin => (x -> (first(x)[1], last(x)[end])))
    rename!(agg_df, :freq_sum => :freq, :bin_function => :bin)
    agg_df.bin = [[x[1], x[2]] for x in agg_df.bin]
    return agg_df
end

function plot_hists_overlay(folder_num, specie, start; last=[], folder="", maxval=8e4)
    data = dict_hists[folder_num]["$specie"]
    data_start = data[start]
    f, ax = plot_hist(data_start, label="$(dict_kdamvals[folder_num][start])", specie="$specie", maxval=maxval)
    if isempty(last)
        theend = length(data)
    else
        theend = last[1]
    end
    for i in start+1:theend
        ax = plot_hist(data[i], loc=ax, label="$(dict_kdamvals[folder_num][i])")
    end
    axislegend(ax, framevisible=true)
    display(f)
    if folder != ""
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        if !isdir(plots_folder)
            mkdir(plots_folder)
        end
        hists_folder = joinpath(plots_folder, "plot_hists")
        if !isdir(hists_folder)
            mkdir(hists_folder)
        end
        save(joinpath(hists_folder, "hists_$specie.png"), f)
    end
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


function plot_prop(df_results, df_props, res_ind, title, threshold_vals, folder; maxval=31856, size=(600, 450), set_thresh=150, tosave=false)
    time_data = df_results[res_ind].time
    prop_data = [df_props[res_ind][i] for i in eachindex(react_names[1:end-1])]
    
    f = Figure(size=size)
    ax = Axis(f[1,1], xlabel = "time", ylabel = "propensity", title="$folder, \n $title", limits=((nothing, nothing), (0, maxval)))#, limits=((2.2e5, 2.6e5),(0, maxval)))
    for i in eachindex(react_names[1:end-1])
        iscatter!(ax, time_data, prop_data[i], label="$(react_names[i])", color=color_list[i])
    end    
    axislegend(ax, framevisible=true)
    if mount_path == "/Users/s2257179/stoch_files/threshold_testing/"
        lines!(ax, range(minimum(time_data), maximum(time_data), length = 2), [threshold_vals[res_ind], threshold_vals[res_ind]], linewidth = 4, color = :black)
    elseif mount_path == "/Users/s2257179/stoch_files/kdam_testing/"
        lines!(ax, range(minimum(time_data), maximum(time_data), length = 2), [set_thresh, set_thresh], linewidth = 4, color = :black)
    end
    
    if tosave
        plots_folder = joinpath(joinpath(mount_path, folder), "plots")
        results_folder = joinpath(plots_folder, "propensities")
        if !isdir(results_folder)
            mkdir(results_folder)
        end
        
        save(joinpath(results_folder, "$title.png"), f)
        # return PlotWrapper(f)
    end

    return f
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


function setup_plot_dicts()
    dict_plot_times = Dict{Int, Any}()
    dict_plot_counts = Dict{Int, Any}()
    dict_plot_results = Dict{Tuple{Int64, Symbol}, Any}()
    dict_plot_hists = Dict{Tuple{Int64, Symbol}, Any}()
    dict_stoch_reacts = Dict{Int, Any}()
    dict_plot_props = Dict{Any, Any}()

    return dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props
end