using JLD2, InteractiveViz, GLMakie, Statistics

@load "/Users/s2257179/Desktop/saved_variables/data_thresh_2.jld2" all_start_indices2 all_stop_indices2 all_switch_rates_on2 all_switch_rates_off2 all_fracs_on2 all_fracs_off2 all_species_mean_on2 all_species_mean_off2
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_5.jld2" all_start_indices5 all_stop_indices5 all_switch_rates_on5 all_switch_rates_off5 all_fracs_on5 all_fracs_off5 all_species_mean_on5 all_species_mean_off5
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_10.jld2" all_start_indices10 all_stop_indices10 all_switch_rates_on10 all_switch_rates_off10 all_fracs_on10 all_fracs_off10 all_species_mean_on10 all_species_mean_off10
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_bs.jld2" all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs thresholds_rtca thresholds_rtcb

# sorting results for plotting
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

# calculate the average and std switch rate but first need to rearrange data so they are in a kdam dict 20 kdam values and in each value there should be 19 numbers 
all_switch_rates_on2[1][1.3]
mean_switch_rates_on2 = Dict()
std_switch_rates_on2 = Dict()
mean_switch_rates_off2 = Dict()
std_switch_rates_off2 = Dict()
mean_switch_rates_on5 = Dict()
std_switch_rates_on5 = Dict()
mean_switch_rates_off5 = Dict()
std_switch_rates_off5 = Dict()
mean_switch_rates_on10 = Dict()
std_switch_rates_on10 = Dict()
mean_switch_rates_off10 = Dict()
std_switch_rates_off10 = Dict()
mean_switch_rates_on_bs = Dict()
std_switch_rates_on_bs = Dict()
mean_switch_rates_off_bs = Dict()
std_switch_rates_off_bs = Dict()

for kdam in kdams
    mean_switch_rates_on2[kdam] = mean([all_switch_rates_on2[i][kdam] for i in eachindex(all_switch_rates_on2)])
    std_switch_rates_on2[kdam] = std([all_switch_rates_on2[i][kdam] for i in eachindex(all_switch_rates_on2)])
    mean_switch_rates_off2[kdam] = mean([all_switch_rates_off2[i][kdam] for i in eachindex(all_switch_rates_off2)])
    std_switch_rates_off2[kdam] = std([all_switch_rates_off2[i][kdam] for i in eachindex(all_switch_rates_off2)])
    mean_switch_rates_on5[kdam] = mean([all_switch_rates_on5[i][kdam] for i in eachindex(all_switch_rates_on5)])
    std_switch_rates_on5[kdam] = std([all_switch_rates_on5[i][kdam] for i in eachindex(all_switch_rates_on5)])
    mean_switch_rates_off5[kdam] = mean([all_switch_rates_off5[i][kdam] for i in eachindex(all_switch_rates_off5)])
    std_switch_rates_off5[kdam] = std([all_switch_rates_off5[i][kdam] for i in eachindex(all_switch_rates_off5)])
    mean_switch_rates_on10[kdam] = mean([all_switch_rates_on10[i][kdam] for i in eachindex(all_switch_rates_on10)])
    std_switch_rates_on10[kdam] = std([all_switch_rates_on10[i][kdam] for i in eachindex(all_switch_rates_on10)])
    mean_switch_rates_off10[kdam] = mean([all_switch_rates_off10[i][kdam] for i in eachindex(all_switch_rates_off10)])
    std_switch_rates_off10[kdam] = std([all_switch_rates_off10[i][kdam] for i in eachindex(all_switch_rates_off10)])
    mean_switch_rates_on_bs[kdam] = mean([all_switch_rates_on_bs[i][kdam] for i in eachindex(all_switch_rates_on_bs)])
    std_switch_rates_on_bs[kdam] = std([all_switch_rates_on_bs[i][kdam] for i in eachindex(all_switch_rates_on_bs)])
end


# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams, [all_switch_rates_on2[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on2)]
[lines!(ax, kdams, [all_switch_rates_off2[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on2)]
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