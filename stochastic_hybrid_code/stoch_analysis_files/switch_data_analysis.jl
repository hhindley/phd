using JLD2, InteractiveViz, GLMakie

@load "/Users/s2257179/phd/stochastic_hybrid_code/saved_data/raw_data.jld2" folders_dict dict_times dict_kdamvals dict_titles dict_results dict_reacts dict_props dict_counts dict_hists
@load "/Users/s2257179/phd/stochastic_hybrid_code/saved_data/data.jld2" all_start_indices all_stop_indices all_switch_rates_on all_switch_rates_off all_fracs_on all_fracs_off all_species_mean_on all_species_mean_off all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs

# sorting results for plotting
kdams = sort(collect(keys(all_switch_rates_on[collect(keys(folders_dict))[1]])))

# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams, [all_switch_rates_on[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on)]
[lines!(ax, kdams, [all_switch_rates_off[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, kdams, [all_fracs_on[i][key] for key in kdams], label="on state $i") for i in eachindex(all_fracs_on)]
[lines!(ax, kdams, [all_fracs_off[i][key] for key in kdams], label="off state $i") for i in eachindex(all_fracs_off)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)

specie_plot = 1
cmap = :rainbow_bgyr_35_85_c72_n256
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[$([species][specie_plot])] (μM)")
[scatter!(ax, kdams, [all_species_mean_on[i][key][specie_plot] for key in kdams], color=[all_fracs_on[i][key] for key in kdams], colorrange=(0,1), colormap=cmap) for i in eachindex(all_species_mean_on)]
[scatter!(ax, kdams, [all_species_mean_off[i][key][specie_plot] for key in kdams], color=[all_fracs_off[i][key] for key in kdams], colorrange=(0,1), colormap=cmap) for i in eachindex(all_species_mean_off)]
Colorbar(f[1, 2], limits = (0,1), colormap = cmap, label="fraction of time in state")