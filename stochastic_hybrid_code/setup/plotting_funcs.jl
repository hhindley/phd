function plotBIG(x, y; xtitle="", ytitle="", title="")
    xvals = range(minimum(x), maximum(x), length=length(x))

    f = Figure()
    ax = Axis(f[1, 1], xlabel = "$xtitle", ylabel = "$ytitle", title="$title")
    ilines!(f[1,1], xvals, y)

    return f
end

struct PlotWrapper
    plot::Any
end
function get_column(df::DataFrame, column_name::Union{Symbol, String})
    return getproperty(df, column_name)
end

function create_subplots(plotting_func, rows, columns; size=(1450, 800), xlabel="", ylabel="", titles=[], hidelabels=[true, true], yscale=identity)
    f = Figure(size=size)

    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            # println(data_ind)
            title = titles != [] ? titles[data_ind] : "Plot $data_ind"
            if plotting_func == "plot_results" || plotting_func == "plot_hists" || plotting_func == "plot_individual_reacts"
                ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale)
            elseif plotting_func == "plot_stoch_reacts"
                ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title, yscale=yscale, xticks=(1:13, react_names_str), xticklabelrotation=45)

            end
            # ilines!(f[i,j], makiex(df_results[data_ind].time), df_results[data_ind][:,species])

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
    return f
end
function add_subplots(f, plotting_func, df_results, rows, columns; species=:rm_a, linkaxes=true)
    if plotting_func == "plot_results"
        preprocessed_times = [get_column(df, :time) for df in df_results]
        preprocessed_results = [get_column(df, species) for df in df_results]
    elseif plotting_func == "plot_hists"
        dfs = df_results[string(species)]
    end
    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            if plotting_func == "plot_results"
                # x = makiex(preprocessed_times[data_ind])
                # ilines!(f[i,j], x, preprocessed_results[data_ind])
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
                # for r in eachindex(react_names[1:end-1])
                #     iscatter!(f[i,j], df_results[data_ind].time, get_column(df_props[data_ind], react_names[r]), label="$(react_names[r])", color=color_list[r])
                # end
                [iscatter!(f[i,j], df_results[data_ind].time, get_column(df_props[data_ind], react_names[r]), label="$(react_names[r])", color=color_list[r]) for r in eachindex(react_names[1:end-1])]
            end
        end
    end
    if linkaxes == true
        linkaxes!(filter(x -> x isa Axis, f.content)...)
    end
    return f
end
function plot_results(plotting_func, df_results, rows, columns; species=:rm_a, size=(1450, 900), xlabel="", ylabel="", titles=[], yscale=identity, hidelabels=[true,true], linkaxes=true, folder="")
    f = create_subplots(plotting_func, rows, columns; size=size, xlabel=xlabel, ylabel=ylabel, titles=titles, hidelabels=hidelabels, yscale=yscale)
    f = add_subplots(f, plotting_func, df_results, rows, columns; species=species, linkaxes=linkaxes)
    if !isempty(folder)
        folder_path = joinpath("/Users/s2257179/phd/stochastic_hybrid_code", folder)
        save(joinpath(folder_path, "$species.png"), f)
        return PlotWrapper(f)
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
    return f
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



function plot_props(df_results, df_props, threshold_vals, max_value=31856.296174439733)
    f = Figure(size=(1450, 800))
    rows = 5
    total_columns = 4

    for j in 1:total_columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            ax = Axis(f[i,j], xlabel = "time", ylabel = "propensity", limits=(nothing,(0, max_value)))
            [iscatter!(ax,df_results[data_ind].time, df_props[data_ind][r], label="$(react_names[r])", color=color_list[r]) for r in eachindex(react_names[1:end-1])]
            axislegend(ax, framevisible=true)
            lines!(range(minimum(df_results[data_ind].time), maximum(df_results[data_ind].time), length=2), [threshold_vals[data_ind], threshold_vals[data_ind]], linewidth=4, color=RGB(0.0, 0.0, 0.0))

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
    linkaxes!(filter(x -> x isa Axis, f.content)...)
    return f
    
end