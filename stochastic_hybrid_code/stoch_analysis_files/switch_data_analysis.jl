using JLD2, InteractiveViz, GLMakie, Statistics

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))

# @load "/Users/s2257179/Desktop/saved_variables/data_thresh_2.jld2" all_start_indices2 all_stop_indices2 all_switch_rates_on2 all_switch_rates_off2 all_fracs_on2 all_fracs_off2 all_species_mean_on2 all_species_mean_off2
# @load "/Users/s2257179/Desktop/saved_variables/data_thresh_5.jld2" all_start_indices5 all_stop_indices5 all_switch_rates_on5 all_switch_rates_off5 all_fracs_on5 all_fracs_off5 all_species_mean_on5 all_species_mean_off5
# @load "/Users/s2257179/Desktop/saved_variables/data_thresh_10.jld2" all_start_indices10 all_stop_indices10 all_switch_rates_on10 all_switch_rates_off10 all_fracs_on10 all_fracs_off10 all_species_mean_on10 all_species_mean_off10
# @load "/Users/s2257179/Desktop/saved_variables/data_thresh_bs.jld2" all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs thresholds_rtca thresholds_rtcb

@load "/Users/s2257179/Desktop/saved_variables/data_thresh_all.jld2" indices switch_rates fracs species_mean thresholds_bs

# sorting results for plotting
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

switch_rates["on"]
fracs

# calculate the average and std switch rate but first need to rearrange data so they are in a kdam dict 20 kdam values and in each value there should be 19 numbers 
mean_switch_frac = Dict("on"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())), 
                        "off"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())))
std_switch_frac = Dict("on"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())), 
                        "off"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())))

for onoff in ["on", "off"]
    for thresh in ["2", "5", "10", "bs"]
        for kdam in kdams
            rates_switch = [switch_rates[onoff][thresh][i][kdam] for i in eachindex(switch_rates[onoff][thresh])]
            mean_switch_frac[onoff][thresh]["switch"][kdam] = mean(rates_switch)
            std_switch_frac[onoff][thresh]["switch"][kdam] = std(rates_switch)
            
            fracs1 = [fracs[onoff][thresh][i][kdam] for i in eachindex(fracs[onoff][thresh])]
            mean_switch_frac[onoff][thresh]["frac"][kdam] = mean(fracs1)
            std_switch_frac[onoff][thresh]["frac"][kdam] = std(fracs1)
        end
    end
end

f_switch2 = plot_mean_std("2", "switch")
f_frac2 = plot_mean_std("2", "frac")

f_switch5 = plot_mean_std("5", "switch")
f_frac5 = plot_mean_std("5", "frac")

f_switch10 = plot_mean_std("10", "switch")
f_frac10 = plot_mean_std("10", "frac")

f_switch_bs = plot_mean_std("bs", "switch")
f_frac_bs = plot_mean_std("bs", "frac")

switch_all = plot_mean_std(["2", "5", "10", "bs"], "switch")
fracs_all = plot_mean_std(["2", "5", "10", "bs"], "frac")

display(GLMakie.Screen(), f)


# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams, [switch_rates["on"]["2"][i][key] for key in kdams], label="on→off $i") for i in eachindex(switch_rates["on"]["2"])]
[lines!(ax, kdams, [switch_rates["off"]["2"][i][key] for key in kdams], label="on→off $i") for i in eachindex(switch_rates["off"]["2"])]
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, kdams, [fracs["on"]["2"][i][key] for key in kdams], label="on state $i") for i in eachindex(fracs["on"]["2"])]
[lines!(ax, kdams, [fracs["off"]["2"][i][key] for key in kdams], label="off state $i") for i in eachindex(fracs["off"]["2"])]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)

# plotting dots of average concentrations with fraction of time in state as colour
f2_rtca = plot_conc_frac(1, "2", logz=true)
f5_rtca = plot_conc_frac(1, "5", logz=true)
f10_rtca = plot_conc_frac(1, "10", logz=true)
f_bs_rtca = plot_conc_frac(1, "bs", logz=true)

f2_rh = plot_conc_frac(3, "2", logz=true)
f5_rh = plot_conc_frac(3, "5", logz=true)
f10_rh = plot_conc_frac(3, "10", logz=true)
f_bs_rh = plot_conc_frac(3, "bs", logz=true)

f2_rm_a = plot_conc_frac(2, "2", logz=true)
f5_rm_a = plot_conc_frac(2, "5", logz=true)
f10_rm_a = plot_conc_frac(2, "10", logz=true)
f_bs_rm_a = plot_conc_frac(2, "bs", logz=true)

display(GLMakie.Screen(), f5)

