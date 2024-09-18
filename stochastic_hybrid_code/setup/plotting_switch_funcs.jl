function plot_mean_std(thresh, switch::String)

    f = Figure()

    if typeof(thresh) == String
        if switch == "switch"
            ax = Axis(f[1, 1], xlabel="kdam", ylabel="switch rate", title="threshold $thresh", yscale=log10)
            labels = ["on → off", "off → on"]
        else
            ax = Axis(f[1, 1], xlabel="kdam", ylabel="fraction of time in state", title="threshold $thresh", yscale=identity)
            labels = ["on state", "off state"]
        end

        sorted_keys = sort(collect(keys(mean_switch_frac["off"][thresh]["switch"])))

        errorbars!(ax, kdams, [mean_switch_frac["off"][thresh][switch][key] for key in sorted_keys], [std_switch_frac["off"][thresh][switch][key] for key in sorted_keys])
        lines!(ax, kdams, [mean_switch_frac["off"][thresh][switch][key] for key in sorted_keys], label=labels[1])
        errorbars!(ax, kdams, [mean_switch_frac["on"][thresh][switch][key] for key in sorted_keys], [std_switch_frac["on"][thresh][switch][key] for key in sorted_keys])
        lines!(ax, kdams, [mean_switch_frac["on"][thresh][switch][key] for key in sorted_keys], label=labels[2])
    else
        if switch == "switch"
            ax = Axis(f[1, 1], xlabel="kdam", ylabel="switch rate", title="threshold $thresh", yscale=log10)
        else
            ax = Axis(f[1, 1], xlabel="kdam", ylabel="fraction of time in state", title="threshold $thresh", yscale=identity)
        end

        sorted_keys = sort(collect(keys(mean_switch_frac["off"][thresh[1]]["switch"])))

        for i in thresh
            if switch == "switch"
                labels = ["$i: on→off", "$i: off→on"]
            else
                labels = ["$i: on", "$i: off"]
            end

            errorbars!(ax, kdams, [mean_switch_frac["off"][i][switch][key] for key in sorted_keys], [std_switch_frac["off"][i][switch][key] for key in sorted_keys])
            lines!(ax, kdams, [mean_switch_frac["off"][i][switch][key] for key in sorted_keys], label=labels[1])
            errorbars!(ax, kdams, [mean_switch_frac["on"][i][switch][key] for key in sorted_keys], [std_switch_frac["on"][i][switch][key] for key in sorted_keys])
            lines!(ax, kdams, [mean_switch_frac["on"][i][switch][key] for key in sorted_keys], label=labels[2])
        end
    end
    axislegend(position=:rc)

    return f
end

function plot_conc_frac(specie_plot, thresh; logz=false)
    specie_plot = specie_plot; species = [:rtca, :rm_a, :rh];
    cmap = :rainbow_bgyr_35_85_c72_n256
    f = Figure()
    ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[$(species[specie_plot])] (μM)", title="threshold $(thresh)")
    
    if !logz
        max_val = maximum([fracs["off"][thresh][i][key] for i in eachindex(fracs["off"][thresh]) for key in kdams])
        min_val = minimum([fracs["on"][thresh][i][key] for i in eachindex(fracs["on"][thresh]) for key in kdams])

        [scatter!(ax, kdams, [species_mean["on"][thresh][i][key][specie_plot] for key in kdams], color=[fracs["on"][thresh][i][key] for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(species_mean["on"][thresh])]
        [scatter!(ax, kdams, [species_mean["off"][thresh][i][key][specie_plot] for key in kdams], color=[fracs["off"][thresh][i][key] for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(species_mean["off"][thresh])]
        Colorbar(f[1, 2], limits = (min_val,max_val), colormap = cmap, label="fraction of time in state")
    else 
        max_val = log10(maximum([fracs["off"][thresh][i][key] for i in eachindex(fracs["off"][thresh]) for key in kdams]))
        min_val = log10(minimum([fracs["on"][thresh][i][key] for i in eachindex(fracs["on"][thresh]) for key in kdams]))
        [scatter!(ax, kdams, [species_mean["on"][thresh][i][key][specie_plot] for key in kdams], color=[log10(fracs["on"][thresh][i][key]) for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(species_mean["on"][thresh])]
        [scatter!(ax, kdams, [species_mean["off"][thresh][i][key][specie_plot] for key in kdams], color=[log10(fracs["off"][thresh][i][key]) for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(species_mean["off"][thresh])]
        # Colorbar(f[1, 2], limits=(min_val,max_val), colormap=cmap, label="log10(fraction of time in state)")
        colorbar_ticks = LinRange(min_val, max_val, 5)
        Colorbar(f[1, 2], limits=(min_val,max_val), colormap=cmap, label="fraction of time in state (logscale)", ticks=(colorbar_ticks, string.(round.(10 .^colorbar_ticks, digits=3))))
    end
    return f
end