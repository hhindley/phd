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

function create_subplots(rows, columns; size=(1450, 900), xlabel="", ylabel="", titles=[], hidelabels=true)
    f = Figure(size=size)

    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            # println(data_ind)
            title = titles != [] ? titles[data_ind] : "Plot $data_ind"

            ax = Axis(f[i, j], xlabel = xlabel, ylabel = ylabel, title=title)
            # ilines!(f[i,j], makiex(df_results[data_ind].time), df_results[data_ind][:,species])

            # Hide x-axis decorations for axes not in the bottom row
            if i != rows && hidelabels == true
                hidexdecorations!(ax, grid=false)
            end

            # Hide y-axis decorations for axes not in the first column
            if j > 1 && hidelabels == true
                hideydecorations!(ax, grid=false)
            end
        end
    end
    return f
end
function add_subplots(f, df_results, species, rows, columns; linkaxes=true)
    preprocessed_times = [df.time for df in df_results]
    preprocessed_results = [df[:,species] for df in df_results]
    for j in 1:columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            x = makiex(preprocessed_times[data_ind])
            ilines!(f[i,j], x, preprocessed_results[data_ind])
        end
    end
    if linkaxes == true
        linkaxes!(filter(x -> x isa Axis, f.content)...)
    end
    return f
end
function plot_results(df_results, species, rows, columns; size=(1450, 900), xlabel="", ylabel="", titles=[], hidelabels=true, linkaxes=true, folder="")
    f = create_subplots(rows, columns; size=size, xlabel=xlabel, ylabel=ylabel, titles=titles, hidelabels=hidelabels)
    f = add_subplots(f, df_results, species, rows, columns; linkaxes=linkaxes)
    if !isempty(folder)
        folder_path = joinpath("/Users/s2257179/phd/stochastic_hybrid_code", folder)
        save(joinpath(folder_path, "$species.png"), f)
        return PlotWrapper(f)
    else
        return f
    end
    
end



function makiex(x)
    return range(minimum(x), maximum(x), length=length(x))
end


function plot_props(df_results)
    max_value = mapreduce(df -> maximum(maximum(df.event)), max, df_results)

    f = Figure(size=(1450, 900))
    rows = 5
    total_columns = 4

    for j in 1:total_columns
        for i in 1:rows
            data_ind = i + rows * (j - 1)
            println(data_ind)
            ax = Axis(f[i, j], xlabel = "time", ylabel = "propensity", title="threshold: $(round(threshold_vals[data_ind], digits=2))")
            for r in names(df_props[data_ind])
                iscatter!(ax, df_results[data_ind].time, df_props[data_ind][:,r])
            end
            ylims!(ax, 0,max_value)

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
    save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/props.png", f)
end