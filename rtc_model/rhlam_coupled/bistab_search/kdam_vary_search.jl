using JLD2
# using GLMakie, InteractiveViz
model = "lamkin_coupled"
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))
lam_c_vals = 10 .^range(log10(1e-8),log10(0.1), length=5) #8e-7
kin_c_vals = 10 .^range(log10(1e-6),log10(0.1), length=5) #1.5e-5
wab_vals = 10 .^range(log10(1e-6),log10(0.1), length=5) #1.3e-5
kdam_range = range(0,100,length=20)

lamkin_ssvals = []
lamkin_vals = []
iteration_num = 0
println("starting loop")
for i in lam_c_vals
    for j in kin_c_vals
        for wab in wab_vals
            iteration_num += 1
            println("iteration: $iteration_num")
            test_params = deepcopy(params_rtc1)
            test_params[lam_c] = i
            test_params[kin_c] = j
            test_params[ω_ab] = wab
            push!(lamkin_vals, (i,j,wab))
            ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
            res = var_param(lamkin_coupled, kdam, test_params, kdam_range, ssvals_rtc_test)
            push!(lamkin_ssvals, res.rtca)
        end
    end
end

println("saving variables")
@save joinpath(homedir(), "phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search.jld2") lamkin_ssvals lamkin_vals
println("finished!")

# @load joinpath(homedir(),"phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search.jld2") lamkin_ssvals lamkin_vals


# plot([scatter(x=kdam_range, y=lamkin_ssvals[i], name="lam_c, kin_c, wab = $(lamkin_vals[i])") for i in eachindex(lamkin_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))




# function bistab_check(ssvals)
#     m, m_ind = findmax(ssvals)
#     for i in ssvals[m_ind:end]
#         if ((i-m)/m*100) < -50
#             return true
#         end
#     end
#     return false
# end
# indices=[]
# for i in eachindex(lamkin_ssvals)
#     if bistab_check(lamkin_ssvals[i])
#         push!(indices, i)
#     end
# end
# indices

# f = Figure()
# ax = Axis(f[1, 1])
# [ilines!(ax, kdam_range, lamkin_ssvals[i]) for i in indices]

# lamkin_vals[indices]

# lamkin_ssvals[indices[4]]
# f = Figure()
# ax = Axis(f[1, 1])
# lines!(ax, kdam_range, lamkin_ssvals[indices[4]])

# lamkin_vals[indices[5]]

# lam_c_val = lamkin_vals[indices[5]][1]
# kin_c_val = lamkin_vals[indices[5]][2]
# ω_ab_val = lamkin_vals[indices[5]][3]

