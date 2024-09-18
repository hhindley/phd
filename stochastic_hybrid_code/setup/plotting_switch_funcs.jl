function plot_mean_std(mean_switch_frac_on, std_switch_frac_on, mean_switch_frac_off, std_switch_frac_off, switch::String)
    sorted_keys = sort(collect(keys(mean_switch_frac_off10["switch"])))

    f = Figure()
    if switch == "switch"
        ax = Axis(f[1, 1], xlabel="kdam", ylabel="switch rate", yscale=log10)
        labels = ["on → off", "off → on"]
    else
        ax = Axis(f[1, 1], xlabel="kdam", ylabel="fraction of time in state", yscale=identity)
        labels = ["on state", "off state"]
    end
    errorbars!(ax, kdams, [mean_switch_frac_off[switch][key] for key in sorted_keys], [std_switch_frac_off[switch][key] for key in sorted_keys])
    lines!(ax, kdams, [mean_switch_frac_off[switch][key] for key in sorted_keys], label=labels[1])
    errorbars!(ax, kdams, [mean_switch_frac_on[switch][key] for key in sorted_keys], [std_switch_frac_on[switch][key] for key in sorted_keys])
    lines!(ax, kdams, [mean_switch_frac_on[switch][key] for key in sorted_keys], label=labels[2])
    axislegend(position=:rc)
    return f
end

function plot_conc_frac(specie_plot, all_species_mean_on, all_species_mean_off, all_fracs_on, all_fracs_off, var_name; logz=false)
    specie_plot = specie_plot; species = [:rtca, :rm_a, :rh];
    cmap = :rainbow_bgyr_35_85_c72_n256
    f = Figure()
    ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[$(species[specie_plot])] (μM)", title="threshold $(var_name)")
    
    if !logz
        max_val = maximum([all_fracs_off[i][key] for i in eachindex(all_fracs_off) for key in kdams])
        min_val = minimum([all_fracs_on[i][key] for i in eachindex(all_fracs_on) for key in kdams])

        [scatter!(ax, kdams, [all_species_mean_on[i][key][specie_plot] for key in kdams], color=[all_fracs_on[i][key] for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(all_species_mean_on)]
        [scatter!(ax, kdams, [all_species_mean_off[i][key][specie_plot] for key in kdams], color=[all_fracs_off[i][key] for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(all_species_mean_off)]
        Colorbar(f[1, 2], limits = (min_val,max_val), colormap = cmap, label="fraction of time in state")
    else 
        max_val = log10(maximum([all_fracs_off[i][key] for i in eachindex(all_fracs_off) for key in kdams]))
        min_val = log10(minimum([all_fracs_on[i][key] for i in eachindex(all_fracs_on) for key in kdams]))
        [scatter!(ax, kdams, [all_species_mean_on[i][key][specie_plot] for key in kdams], color=[log10(all_fracs_on[i][key]) for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(all_species_mean_on)]
        [scatter!(ax, kdams, [all_species_mean_off[i][key][specie_plot] for key in kdams], color=[log10(all_fracs_off[i][key]) for key in kdams], colorrange=(min_val,max_val), colormap=cmap) for i in eachindex(all_species_mean_off)]
        # Colorbar(f[1, 2], limits=(min_val,max_val), colormap=cmap, label="log10(fraction of time in state)")
        colorbar_ticks = LinRange(min_val, max_val, 5)
        Colorbar(f[1, 2], limits=(min_val,max_val), colormap=cmap, label="fraction of time in state (logscale)", ticks=(colorbar_ticks, string.(round.(10 .^colorbar_ticks, digits=3))))
    end
    return f
end