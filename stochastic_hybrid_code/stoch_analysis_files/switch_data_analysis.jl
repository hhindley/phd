using JLD2, InteractiveViz, GLMakie, Statistics

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))
# kdam data
@load "/Users/s2257179/Desktop/saved_variables/data_thresh_all.jld2" indices switch_rates fracs species_mean thresholds_bs
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
mean_switch_frac, std_switch_frac = calc_mean_std_vars(switch_rates, fracs, kdams)

f_switch2 = plot_mean_std(mean_switch_frac, std_switch_frac, "2", "switch")
f_frac2 = plot_mean_std(mean_switch_frac, std_switch_frac, "2", "frac")

f_switch5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "switch")
f_frac5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "frac")

f_switch10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "switch")
f_frac10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "frac")

f_switch_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "switch")
f_frac_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "frac")

switch_all = plot_mean_std(mean_switch_frac, std_switch_frac, ["2", "5", "10", "bs"], "switch")
fracs_all = plot_mean_std(mean_switch_frac, std_switch_frac, ["2", "5", "10", "bs"], "frac")

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
f2_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "2", logz=true)
f5_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "5", logz=true)
f10_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "10", logz=true)
f_bs_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "bs", logz=true)

f2_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "2", logz=true)
f5_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "5", logz=true)
f10_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "10", logz=true)
f_bs_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "bs", logz=true)

f2_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "2", logz=true)
f5_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "5", logz=true)
f10_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "10", logz=true)
f_bs_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "bs", logz=true)

save("/Users/s2257179/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Documents/rtc/stochastic/plots/analysis/switching/thresh_bs/conc_frac_rtca_bs.png", f_bs_rtca)


