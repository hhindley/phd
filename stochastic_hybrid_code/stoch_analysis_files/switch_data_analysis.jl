using JLD2, InteractiveViz, GLMakie, Statistics

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))

@load "/Users/s2257179/Desktop/saved_variables/data_thresh_2.jld2" all_start_indices2 all_stop_indices2 all_switch_rates_on2 all_switch_rates_off2 all_fracs_on2 all_fracs_off2 all_species_mean_on2 all_species_mean_off2
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_5.jld2" all_start_indices5 all_stop_indices5 all_switch_rates_on5 all_switch_rates_off5 all_fracs_on5 all_fracs_off5 all_species_mean_on5 all_species_mean_off5
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_10.jld2" all_start_indices10 all_stop_indices10 all_switch_rates_on10 all_switch_rates_off10 all_fracs_on10 all_fracs_off10 all_species_mean_on10 all_species_mean_off10
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_bs.jld2" all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs thresholds_rtca thresholds_rtcb

@load "/Users/s2257179/Desktop/saved_variables/data_thresh_all.jld2" all_start_indices all_stop_indices all_switch_rates_on all_switch_rates_off all_fracs_on all_fracs_off all_species_mean_on all_species_mean_off thresholds_rtca thresholds_rtcb

# sorting results for plotting
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]


# calculate the average and std switch rate but first need to rearrange data so they are in a kdam dict 20 kdam values and in each value there should be 19 numbers 
mean_switch_frac_on = Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict()))
std_switch_frac_on = Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict()))

for thresh in eachindex(mean_switch_frac_on)
    for kdam in kdams
        mean_switch_frac_on[thresh]["switch"][kdam] = mean([all_switch_rates_on[thresh][i][kdam] for i in eachindex(all_switch_rates_on[thresh])])
        std_switch_frac_on[thresh]["switch"][kdam] = std([all_switch_rates_on[thresh][i][kdam] for i in eachindex(all_switch_rates_on[thresh])])
        mean_switch_frac_on[thresh]["frac"][kdam] = mean([all_fracs_on[thresh][i][kdam] for i in eachindex(all_fracs_on[thresh])])
        std_switch_frac_on[thresh]["frac"][kdam] = std([all_fracs_on[thresh][i][kdam] for i in eachindex(all_fracs_on[thresh])])
    end
end

mean_switch_frac_on10["frac"]
sorted_keys = sort(collect(keys(mean_switch_frac_off10["switch"])))
sorted_mean_values = [mean_switch_frac_off10["switch"][key] for key in sorted_keys]
sorted_std_values = [std_switch_frac_off10["switch"][key] for key in sorted_keys]



f_switch2 = plot_mean_std(mean_switch_frac_on2, std_switch_frac_on2, mean_switch_frac_off2, std_switch_frac_off2, "switch")
f_frac2 = plot_mean_std(mean_switch_frac_on2, std_switch_frac_on2, mean_switch_frac_off2, std_switch_frac_off2, "frac")

f_switch5 = plot_mean_std(mean_switch_frac_on5, std_switch_frac_on5, mean_switch_frac_off5, std_switch_frac_off5, "switch")
f_frac5 = plot_mean_std(mean_switch_frac_on5, std_switch_frac_on5, mean_switch_frac_off5, std_switch_frac_off5, "frac")

f_switch10 = plot_mean_std(mean_switch_frac_on10, std_switch_frac_on10, mean_switch_frac_off10, std_switch_frac_off10, "switch")
f_frac10 = plot_mean_std(mean_switch_frac_on10, std_switch_frac_on10, mean_switch_frac_off10, std_switch_frac_off10, "frac")

f_switch_bs = plot_mean_std(mean_switch_frac_on_bs, std_switch_frac_on_bs, mean_switch_frac_off_bs, std_switch_frac_off_bs, "switch")
f_frac_bs = plot_mean_std(mean_switch_frac_on_bs, std_switch_frac_on_bs, mean_switch_frac_off_bs, std_switch_frac_off_bs, "frac")


# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams, [all_switch_rates_on2[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on2)]
[lines!(ax, kdams, [all_switch_rates_off2[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on2)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, kdams, [all_fracs_on2[i][key] for key in kdams], label="on state $i") for i in eachindex(all_fracs_on2)]
[lines!(ax, kdams, [all_fracs_off2[i][key] for key in kdams], label="off state $i") for i in eachindex(all_fracs_off2)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)




# delete!(all_species_mean_on2, 10)

f2_rtca = plot_conc_frac(1, all_species_mean_on2, all_species_mean_off2, all_fracs_on2, all_fracs_off2, "2", logz=true)
f5_rtca = plot_conc_frac(1, all_species_mean_on5, all_species_mean_off5, all_fracs_on5, all_fracs_off5, "5", logz=true)
f10_rtca = plot_conc_frac(1, all_species_mean_on10, all_species_mean_off10, all_fracs_on10, all_fracs_off10, "10", logz=true)
f_bs_rtca = plot_conc_frac(1, all_species_mean_on_bs, all_species_mean_off_bs, all_fracs_on_bs, all_fracs_off_bs, "bs", logz=true)

f2_rh = plot_conc_frac(3, all_species_mean_on2, all_species_mean_off2, all_fracs_on2, all_fracs_off2, "2")
f5_rh = plot_conc_frac(3, all_species_mean_on5, all_species_mean_off5, all_fracs_on5, all_fracs_off5, "5")
f10_rh = plot_conc_frac(3, all_species_mean_on10, all_species_mean_off10, all_fracs_on10, all_fracs_off10, "10")
f_bs_rh = plot_conc_frac(3, all_species_mean_on_bs, all_species_mean_off_bs, all_fracs_on_bs, all_fracs_off_bs, "bs")

f2_rm_a = plot_conc_frac(2, all_species_mean_on2, all_species_mean_off2, all_fracs_on2, all_fracs_off2, "2")
f5_rm_a = plot_conc_frac(2, all_species_mean_on5, all_species_mean_off5, all_fracs_on5, all_fracs_off5, "5")
f10_rm_a = plot_conc_frac(2, all_species_mean_on10, all_species_mean_off10, all_fracs_on10, all_fracs_off10, "10")
f_bs_rm_a = plot_conc_frac(2, all_species_mean_on_bs, all_species_mean_off_bs, all_fracs_on_bs, all_fracs_off_bs, "bs")

display(GLMakie.Screen(), f5)

all_species_mean_on2[5]