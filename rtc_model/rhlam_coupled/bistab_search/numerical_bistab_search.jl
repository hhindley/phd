using JLD2
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/bistab_search/setup.jl"))

lamkin_vals = []
inc = []
dec = []

for i in lam_c_vals
    for j in kin_c_vals
        for wab in wab_vals
            test_params = deepcopy(params_rtc1)
            test_params[lam_c] = i
            test_params[kin_c] = j
            test_params[ω_ab] = wab
            push!(lamkin_vals, (i,j,wab))
            res_ss = numerical_bistability_analysis(test, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
            res1_ss = numerical_bistability_analysis(test, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)
            push!(inc, res_ss)
            push!(dec, res1_ss)
        end
    end
end

@save "/home/hollie_hindley/phd/rtc_model/rhlam_coupled/bistab_search/numerical_bistab_search.jld2" lamkin_vals inc dec 

# all_res = vcat(inc, dec)
# inc_lamkin_vals = [("inc",lamkin_vals[i]...) for i in eachindex(lamkin_vals)]
# dec_lamkin_vals = [("dec",lamkin_vals[i]...) for i in eachindex(lamkin_vals)]
# all_lamkinvals = vcat(inc_lamkin_vals, dec_lamkin_vals)
# kdam_inc = fill(kdam_range, length(inc))
# kdam_dec = fill(reverse(kdam_range), length(dec))
# kdams = vcat(kdam_inc, kdam_dec)
# plot([scatter(x=kdams[i], y=all_res[i], name="lam_c, kin_c, wab = $(all_lamkinvals[i])") for i in eachindex(all_res)])
