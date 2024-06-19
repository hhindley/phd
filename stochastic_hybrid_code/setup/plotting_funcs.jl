function plotBIG(x, y; xtitle="", ytitle="", title="")
    xvals = range(minimum(x), maximum(x), length=length(x))

    f = Figure()
    ax = Axis(f[1, 1], xlabel = "$xtitle", ylabel = "$ytitle", title="$title")
    ilines!(f[1,1], xvals, y)

    return f
end


function plot_subplots(df_results, species, threshold_vals)

    f = Figure(size=(1450, 900))
    total_rows = 5
    total_columns = 4
    
    for j in 1:total_columns
        for i in 1:total_rows
            data_ind = i + total_rows * (j - 1)
            println(data_ind)
            ax = Axis(f[i, j], xlabel = "time", ylabel = "$species", title="threshold: $(round(threshold_vals[data_ind], digits=2))")
            ilines!(f[i,j], makiex(df_results[data_ind].time), df_results[data_ind][:,species])
    
            # Hide x-axis decorations for axes not in the bottom row
            if i != total_rows
                hidexdecorations!(ax, grid=false)
            end
    
            # Hide y-axis decorations for axes not in the first column
            if j > 1
                hideydecorations!(ax, grid=false)
            end
        end
    end
    save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/$species.png", f)

    return f
end

function plot_subplots_hists(df_reacts, threshold_vals)

    f = Figure(size=(1450, 900))
    total_rows = 5
    total_columns = 4
    for j in 1:total_columns
        for i in 1:total_rows
            data_ind = i + total_rows * (j - 1)
            println(data_ind)
            ax = Axis(f[i,j], ylabel="reaction count", xticks=(1:length(df_reacts[1].event), string.(collect(df_reacts[1].reaction))),xticklabelrotation=45, title= title="threshold: $(round(threshold_vals[data_ind], digits=2))")#, yscale=log10)
            barplot!(ax, 1:length(collect(df_reacts[data_ind].event)), df_reacts[data_ind].count)
            # ylims!(ax, 0, 1000000)
            # Hide x-axis decorations for axes not in the bottom row
            if i != total_rows
                hidexdecorations!(ax, grid=false)
            end
    
            # Hide y-axis decorations for axes not in the first column
            # if j > 1
            #     hideydecorations!(ax, grid=false)
            # end
            
        end
    end
    linkaxes!(filter(x -> x isa Axis, f.content)...)
    save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/stoch_reacts.png", f)

    return f
end

function makiex(x)
    return range(minimum(x), maximum(x), length=length(x))
end