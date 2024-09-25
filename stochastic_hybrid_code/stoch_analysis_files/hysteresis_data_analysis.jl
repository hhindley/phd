using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

# hysteresis data
@load "/Users/s2257179/Desktop/saved_variables/hysteresis_high.jld2" indices switch_rates fracs species_mean thresholds_bs
indices_high = indices; switch_rates_high = switch_rates; fracs_high = fracs; species_mean_high = species_mean; thresholds_bs_high = thresholds_bs;
@load "/Users/s2257179/Desktop/saved_variables/hysteresis_low.jld2" indices switch_rates fracs species_mean thresholds_bs
indices_low = indices; switch_rates_low = switch_rates; fracs_low = fracs; species_mean_low = species_mean; thresholds_bs_low = thresholds_bs;

kdams_high = [1.5, 0.8]
mean_switch_frac_high, std_switch_frac_high = calc_mean_std_vars(switch_rates_high, fracs_high, kdams_high)
kdams_low = [0.01, 0.8]
mean_switch_frac_low, std_switch_frac_low = calc_mean_std_vars(switch_rates_low, fracs_low, kdams_low)

f_switch2_high = plot_mean_std(mean_switch_frac_high, std_switch_frac_high, "2", "switch")
f_switch2_low = plot_mean_std(mean_switch_frac_low, std_switch_frac_low, "2", "switch")

f_frac2 = plot_mean_std(mean_switch_frac, std_switch_frac, "2", "frac")

f_switch5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "switch")
f_frac5 = plot_mean_std(mean_switch_frac, std_switch_frac, "5", "frac")

f_switch10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "switch")
f_frac10 = plot_mean_std(mean_switch_frac, std_switch_frac, "10", "frac")

f_switch_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "switch")
f_frac_bs = plot_mean_std(mean_switch_frac, std_switch_frac, "bs", "frac")

switch_all = plot_mean_std(mean_switch_frac, std_switch_frac, ["2", "5", "10", "bs"], "switch")
fracs_all = plot_mean_std(mean_switch_frac, std_switch_frac, ["2", "5", "10", "bs"], "frac")



# tot_mean_species = Dict("high"=>Dict("on"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict()), "off"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict())),
#                         "low"=>Dict("on"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict()), "off"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict())))
# tot_std_species = Dict("on"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict()), 
#                         "off"=>Dict("2"=>Dict(), "5"=>Dict(), "10"=>Dict(), "bs"=>Dict()))
# tot_mean_species["on"]["2"]

# species_mean_high["on"]["2"][1][0.8][1]
# for onoff in ["on", "off"]
#     for thresh in ["2", "5", "10", "bs"]
#         for kdam in eachindex(kdams_high)
#             concs_high = [species_mean_high[onoff][thresh][i][kdams_high[kdam]][1] for i in eachindex(species_mean_high[onoff][thresh])]
#             tot_mean_species["high"][onoff][thresh][kdams_high[kdam]] = mean(concs_high)
#             # tot_std_species[onoff][thresh][kdam] = std(rates_switch)
#             concs_low = [species_mean_low[onoff][thresh][i][kdams_low[kdam]][1] for i in eachindex(species_mean_low[onoff][thresh])]
#             tot_mean_species["low"][onoff][thresh][kdams_low[kdam]] = mean(concs_low)
#         end
#     end
# end

# tot_mean_species["high"]["on"]["2"]
# tot_mean_species["low"]["on"]["2"]
# vcat(kdams_low, reverse(kdams_high))

# lows = [tot_mean_species["low"]["on"]["2"][kdam] for kdam in kdams_low]

# highs = [tot_mean_species["high"]["on"]["2"][kdam] for kdam in kdams_high]


# f = Figure()
# ax = Axis(f[1,1], xticks=([0.01, 0.8, 0.8, 1.5], ["0.01", "0.8", "0.8", "1.5"]))
# barplot!(ax, vcat(kdams_low, reverse(kdams_high)), vcat(lows, highs))


# f = Figure()
# ax = Axis(f[1,1])
# [scatter!(ax, kdams_high, [species_mean_high["on"]["2"][i][key][1] for key in kdams_high], color=:blue) for i in eachindex(species_mean_high["on"]["2"])]
# [scatter!(ax, kdams_low, [species_mean_low["on"]["2"][i][key][1] for key in kdams_low], color=:red) for i in eachindex(species_mean_low["on"]["2"])]

df_concs_on = DataFrame(
    "data" => vcat(
        [species_mean_low["on"]["2"][i][0.01][1] for i in eachindex(species_mean_low["on"]["2"])],
        [species_mean_low["on"]["2"][i][0.8][1] for i in eachindex(species_mean_low["on"]["2"])],
        [species_mean_high["on"]["2"][i][0.8][1] for i in eachindex(species_mean_high["on"]["2"])],
        [species_mean_high["on"]["2"][i][1.5][1] for i in eachindex(species_mean_high["on"]["2"])]
    ),
    "group" => repeat([1, 2, 3, 4], inner=length(species_mean_low["on"]["2"]))
)
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "low_0.8", "high_0.8", "1.5"]), xlabel="kdam", ylabel="RtcA in on state (μM)", title="Hysteresis experiement")
boxplot!(ax, df_concs_on.group, df_concs_on.data)

df_concs_off = DataFrame(
    "data" => vcat(
        [species_mean_low["off"]["2"][i][0.01][1] for i in eachindex(species_mean_low["off"]["2"])],
        [species_mean_low["off"]["2"][i][0.8][1] for i in eachindex(species_mean_low["off"]["2"])],
        [species_mean_high["off"]["2"][i][0.8][1] for i in eachindex(species_mean_high["off"]["2"])],
        [species_mean_high["off"]["2"][i][1.5][1] for i in eachindex(species_mean_high["off"]["2"])]
    ),
    "group" => repeat([1, 2, 3, 4], inner=length(species_mean_low["off"]["2"]))
)
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "low_0.8", "high_0.8", "1.5"]), xlabel="kdam", ylabel="RtcA in off state (μM)", title="Hysteresis experiement")
boxplot!(ax, df_concs_off.group, df_concs_off.data)

df_fracs_on = DataFrame(
    "data" => vcat(
        [fracs_low["on"]["2"][i][0.01] for i in eachindex(fracs_low["on"]["2"])],
        [fracs_low["on"]["2"][i][0.8] for i in eachindex(fracs_low["on"]["2"])],
        [fracs_high["on"]["2"][i][0.8] for i in eachindex(fracs_high["on"]["2"])],
        [fracs_high["on"]["2"][i][1.5] for i in eachindex(fracs_high["on"]["2"])]
    ),
    "group" => repeat([1, 2, 3, 4], inner=length(fracs_low["on"]["2"]))
)
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "low_0.8", "high_0.8", "1.5"]), xlabel="kdam", ylabel="Fraction of time in on state", title="Hysteresis experiement")
boxplot!(ax, df_fracs_on.group, df_fracs_on.data)

df_fracs_off = DataFrame(
    "data" => vcat(
        [fracs_low["off"]["2"][i][0.01] for i in eachindex(fracs_low["off"]["2"])],
        [fracs_low["off"]["2"][i][0.8] for i in eachindex(fracs_low["off"]["2"])],
        [fracs_high["off"]["2"][i][0.8] for i in eachindex(fracs_high["off"]["2"])],
        [fracs_high["off"]["2"][i][1.5] for i in eachindex(fracs_high["off"]["2"])]
    ),
    "group" => repeat([1, 2, 3, 4], inner=length(fracs_low["off"]["2"]))
)
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "0.8", "0.8", "1.5"]), xlabel="Damage rate (min⁻¹)", ylabel="Fraction of time in off state", title="Hysteresis experiement")
boxplot!(ax, df_fracs_off.group, df_fracs_off.data)

df_means_concs_on = combine(groupby(df_concs_on, :group), :data => mean => :mean_data)
df_means_concs_on.color = ["red", "purple", "purple", "blue"]
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "0.8", "0.8", "1.5"]), xlabel="Damage rate (min⁻¹)", ylabel="average RtcA in on state (μM)", title="Hysteresis experiement")
barplot!(ax, df_means_concs_on.group, df_means_concs_on.mean_data, color=df_means_concs_on.color)

df_means_concs_off = combine(groupby(df_concs_off, :group), :data => mean => :mean_data)
df_means_concs_off.color = ["red", "purple", "purple", "blue"]
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "0.8", "0.8", "1.5"]), xlabel="Damage rate (min⁻¹)", ylabel="average RtcA in off state (μM)", title="Hysteresis experiement")
barplot!(ax, df_means_concs_off.group, df_means_concs_off.mean_data, color=df_means_concs_off.color)

df_means_fracs_on = combine(groupby(df_fracs_on, :group), :data => mean => :mean_data)
df_means_fracs_on.color = ["red", "purple", "purple", "blue"]
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "0.8", "0.8", "1.5"]), xlabel="Damage rate (min⁻¹)", ylabel="Fraction of time in on state", title="Hysteresis experiement")
barplot!(ax, df_means_fracs_on.group, df_means_fracs_on.mean_data, color=df_means_fracs_on.color)

df_means_fracs_off = combine(groupby(df_fracs_off, :group), :data => mean => :mean_data)
df_means_fracs_off.color = ["red", "purple", "purple", "blue"]
f = Figure()
ax = Axis(f[1,1], xticks=(1:4, ["0.01", "0.8", "0.8", "1.5"]), xlabel="Damage rate (min⁻¹)", ylabel="Fraction of time in off state", title="Hysteresis experiement")
barplot!(ax, df_means_fracs_off.group, df_means_fracs_off.mean_data, color=df_means_fracs_off.color)






# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams_high, [switch_rates_high["on"]["2"][i][key] for key in kdams_high], label="on→off $i") for i in eachindex(switch_rates_high["on"]["2"])]
[lines!(ax, kdams_high, [switch_rates_high["off"]["2"][i][key] for key in kdams_high], label="on→off $i") for i in eachindex(switch_rates_high["off"]["2"])]
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, kdams, [fracs["on"]["2"][i][key] for key in kdams], label="on state $i") for i in eachindex(fracs["on"]["2"])]
[lines!(ax, kdams, [fracs["off"]["2"][i][key] for key in kdams], label="off state $i") for i in eachindex(fracs["off"]["2"])]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)

# plotting dots of average concentrations with fraction of time in state as colour
f2_rtca = plot_conc_frac(species_mean_high, fracs_high, kdams_high, 1, "2", logz=true)
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


