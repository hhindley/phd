using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

high_kdam = true
if high_kdam
    @load "/Users/s2257179/Desktop/saved_variables/0210/high_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs
else
    @load "/Users/s2257179/Desktop/saved_variables/0210/low_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs
end

# kdam data
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
mean_switch_frac, std_switch_frac = calc_mean_std_vars(switch_rates, fracs, kdams)

f_switch2 = plot_mean_std(mean_switch_frac, std_switch_frac, "2", "switch", tosave=true)
f_frac2 = plot_mean_std(mean_switch_frac, std_switch_frac, "2", "frac", tosave=true)

f_switch5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "switch", tosave=true)
f_frac5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "frac", tosave=true)

f_switch10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "switch", tosave=true)
f_frac10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "frac", tosave=true)

f_switch_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "switch", tosave=true)
f_frac_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "frac", tosave=true)

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
logz = true
includeoff = true
tosave = true

f2_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "2", logz=logz, includeoff=includeoff, tosave=tosave)
f5_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "5", logz=logz, includeoff=includeoff, tosave=tosave)
f10_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "10", logz=logz, includeoff=includeoff, tosave=tosave)
f_bs_rtca = plot_conc_frac(species_mean, fracs, kdams, 1, "bs", logz=logz, includeoff=includeoff, tosave=tosave)

f2_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "2", logz=logz, includeoff=includeoff, tosave=tosave)
f5_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "5", logz=logz, includeoff=includeoff, tosave=tosave)
f10_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "10", logz=logz, includeoff=includeoff, tosave=tosave)
f_bs_rh = plot_conc_frac(species_mean, fracs, kdams, 3, "bs", logz=logz, includeoff=includeoff, tosave=tosave)

f2_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "2", logz=logz, includeoff=includeoff, tosave=tosave)
f5_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "5", logz=logz, includeoff=includeoff, tosave=tosave)
f10_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "10", logz=logz, includeoff=includeoff, tosave=tosave)
f_bs_rm_a = plot_conc_frac(species_mean, fracs, kdams, 2, "bs", logz=logz, includeoff=includeoff, tosave=tosave)







@load "/Users/s2257179/Desktop/saved_variables/0210/high_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs
species_mean_high = species_mean; switch_rates_high = switch_rates; fracs_high = fracs;
@load "/Users/s2257179/Desktop/saved_variables/0210/low_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs
species_mean_low = species_mean; switch_rates_low = switch_rates; fracs_low = fracs;

mean_switch_frac_high, std_switch_frac_high = calc_mean_std_vars(switch_rates_high, fracs_high, kdams)
mean_switch_frac_low, std_switch_frac_low = calc_mean_std_vars(switch_rates_low, fracs_low, kdams)

ordered_keys = sort(collect(keys(mean_switch_frac_high["on"]["2"]["frac"])))
ordered_fracs_high = [mean_switch_frac_high["on"]["2"]["frac"][key] for key in ordered_keys]
ordered_fracs_low = [mean_switch_frac_low["on"]["2"]["frac"][key] for key in ordered_keys]
log_fracs_high = log.(ordered_fracs_high)
log_fracs_low = log.(ordered_fracs_low)
combined_log_fracs = vcat(log_fracs_low, log_fracs_high)
min_log_frac = minimum(combined_log_fracs)
max_log_frac = maximum(combined_log_fracs)
normalized_log_fracs_high = (log_fracs_high .- min_log_frac) ./ (max_log_frac_high - min_log_frac)
normalized_log_fracs_low = (log_fracs_low .- min_log_frac) ./ (max_log_frac_low - min_log_frac)
colors_high = [get(ColorSchemes.rainbow_bgyr_35_85_c72_n256, frac) for frac in normalized_log_fracs_high]
colors_low = [get(ColorSchemes.rainbow_bgyr_35_85_c72_n256, frac) for frac in normalized_log_fracs_low]


df_concs_on_low = DataFrame(
    "data" => 
        [species_mean_low["on"]["2"][i][kdams[j]][1] for j in eachindex(kdams) for i in eachindex(species_mean_low["on"]["2"])],
    "kdam" => repeat(kdams, inner=length(species_mean_low["on"]["2"])),
    "group" => repeat(range(1,20,length=20), inner=length(species_mean_low["on"]["2"])),
    "color" => repeat(colors_low, inner=length(species_mean_low["on"]["2"]))
)

df_concs_on_high = DataFrame(
    "data" => 
        [species_mean_high["on"]["2"][i][kdams[j]][1] for j in eachindex(kdams) for i in eachindex(species_mean_high["on"]["2"])],
    "kdam" => repeat(kdams, inner=length(species_mean_high["on"]["2"])),
    "group" => repeat(range(1,20,length=20), inner=length(species_mean_high["on"]["2"])),
    "color" => repeat(colors_high, inner=length(species_mean_high["on"]["2"]))
)


f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min⁻¹)", ylabel="RtcA in on state (μM)", title="Hysteresis experiement")
violin!(ax, df_concs_on_low.kdam, df_concs_on_low.data, side=:left, color=df_concs_on_low.color)
violin!(ax, df_concs_on_high.kdam, df_concs_on_high.data, side=:right, color=df_concs_on_high.color)


f = Figure() 
ax = Axis(f[1,1], xlabel="Damage rate (min⁻¹)", ylabel="RtcA in on state (μM)", title="Hysteresis experiement")
boxplot!(ax, df_concs_on_low.group, df_concs_on_low.data, color=df_concs_on_low.color)
boxplot!(ax, df_concs_on_high.kdam, df_concs_on_high.data, side=:right, color=df_concs_on_high.color)












using ModelingToolkit, DifferentialEquations, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))


br = get_br(rtc_model, ssvals_rtc, params_rtc)
df = create_br_df(br)

lines!(df.kdam, df.rtca, linewidth=3)
