using JLD2
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/bistab_search/setup.jl"))


lamkin_ssvals = []
lamkin_vals = []

for i in lam_c_vals
    for j in kin_c_vals
        for wab in wab_vals
            test_params = deepcopy(params_rtc1)
            test_params[lam_c] = i
            test_params[kin_c] = j
            test_params[ω_ab] = wab
            push!(lamkin_vals, (i,j,wab))
            res = var_param(test, kdam, test_params, kdam_range, ssvals_rtc)
            push!(lamkin_ssvals, res.rtca)
        end
    end
end


@save "/home/hollie_hindley/phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search.jld2" lamkin_ssvals lamkin_vals


# @load "/home/hollie_hindley/phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search.jld2" lamkin_ssvals lamkin_vals


# plot([scatter(x=kdam_range, y=lamkin_ssvals[i], name="lam_c, kin_c, wab = $(lamkin_vals[i])") for i in eachindex(lamkin_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))
