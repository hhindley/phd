using JLD2, GLMakie, InteractiveViz

# decide which model you want to test and set it as variable name
model = "lamkin_coupled"
# model = "lam_coupled"

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
            res_ss = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
            res1_ss = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)
            push!(inc, res_ss)
            push!(dec, res1_ss)
        end
    end
end

@save "/home/hollie_hindley/phd/rtc_model/rhlam_coupled/bistab_search/numerical_bistab_search_smaller.jld2" lamkin_vals inc dec 


@load joinpath(homedir(), "phd/rtc_model/rhlam_coupled/bistab_search/numerical_bistab_search_smaller.jld2") lamkin_vals inc dec
all_res = vcat(inc, dec)
inc_lamkin_vals = [("inc",lamkin_vals[i]...) for i in eachindex(lamkin_vals)]
dec_lamkin_vals = [("dec",lamkin_vals[i]...) for i in eachindex(lamkin_vals)]
all_lamkinvals = vcat(inc_lamkin_vals, dec_lamkin_vals)
kdam_inc = fill(kdam_range, length(inc))
kdam_dec = fill(reverse(kdam_range), length(dec))
kdams = vcat(kdam_inc, kdam_dec)
# plot([scatter(x=kdams[i], y=all_res[i], name="lam_c, kin_c, wab = $(all_lamkinvals[i])") for i in eachindex(all_res)])

digit = 4
diff_res = []
for i in eachindex(inc)
    if round.(inc[i], digits=digit) != round.(reverse(dec[i]), digits=digit)
        push!(diff_res, i)
    end
end
diff_res

f = Figure()
ax=Axis(f[1,1])
[lines!(ax, kdams[i], all_res[i]) for i in eachindex(all_res)]

f = Figure()
ax=Axis(f[1,1])
[lines!(ax, kdam_inc[i], inc[i]) for i in diff_res]
[lines!(ax, kdam_dec[i], dec[i]) for i in diff_res]

round.(inc[1], digits=6)[1]
reverse(dec[1])[1]

diff_res

x=9
f = Figure()
ax=Axis(f[1,1])
lines!(ax, kdam_inc[diff_res[x]], inc[diff_res[x]])
lines!(ax, kdam_dec[diff_res[x]], dec[diff_res[x]])


for i in diff_res
    println(i)
    f = Figure()
    ax=Axis(f[1,1])
    lines!(ax, kdam_inc[i], inc[i])
    lines!(ax, kdam_dec[i], dec[i])
    display(f)
    sleep(2)
end

inds = [13,26,41,58,67,68,95]
kdam_range = range(0,60, length=50)
ind=7
lamkin_vals[inds[ind]]
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[inds[ind]][1]
test_params[kin_c] = lamkin_vals[inds[ind]][2]
test_params[ω_ab] = lamkin_vals[inds[ind]][3]
res_ss = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
res1_ss = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)

f = Figure()
ax=Axis(f[1,1])
lines!(ax, kdam_range, res_ss)
lines!(ax, reverse(kdam_range), res1_ss)

interesting_inds = [26,41,58,67,95]

lamkin_vals[26]
lamkin_vals[41]
lamkin_vals[58]
lamkin_vals[67]
lamkin_vals[95]



# 26
kdam_range26 = range(0,20, length=30)
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[26][1]
test_params[kin_c] = lamkin_vals[26][2]
test_params[ω_ab] = lamkin_vals[26][3]
res_ss26 = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range26, kdam; specie=:rtca)
res1_ss26 = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range26), kdam, specie=:rtca)

f26 = Figure()
ax26=Axis(f26[1,1])
lines!(ax26, kdam_range26, res_ss26)
lines!(ax26, reverse(kdam_range26), res1_ss26)

# 41
kdam_range41 = range(0,50, length=30)
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[41][1]
test_params[kin_c] = lamkin_vals[41][2]
test_params[ω_ab] = lamkin_vals[41][3]
res_ss41 = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range41, kdam; specie=:rtca)
res1_ss41 = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range41), kdam, specie=:rtca)

f41 = Figure()
ax41=Axis(f41[1,1])
lines!(ax41, kdam_range41, res_ss41)
lines!(ax41, reverse(kdam_range41), res1_ss41)

# 58
kdam_range58 = range(0,20, length=30)
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[58][1]
test_params[kin_c] = lamkin_vals[58][2]
test_params[ω_ab] = lamkin_vals[58][3]
res_ss58 = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
res1_ss58 = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)

f58 = Figure()
ax58=Axis(f58[1,1])
lines!(ax58, kdam_range58, res_ss58)
lines!(ax58, reverse(kdam_range58), res1_ss58)

# 67
kdam_range67 = range(0,20, length=30)
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[67][1]
test_params[kin_c] = lamkin_vals[67][2]
test_params[ω_ab] = lamkin_vals[67][3]
res_ss67 = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
res1_ss67 = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)

f67 = Figure()
ax67=Axis(f67[1,1])
lines!(ax67, kdam_range67, res_ss67)
lines!(ax67, reverse(kdam_range67), res1_ss67)

# 95
kdam_range95 = range(0,20, length=30)
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[95][1]
test_params[kin_c] = lamkin_vals[95][2]
test_params[ω_ab] = lamkin_vals[95][3]
res_ss95 = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam; specie=:rtca)
res1_ss95 = numerical_bistability_analysis(model, test_params, ssvals_rtc, species_rtc, reverse(kdam_range), kdam, specie=:rtca)

f95 = Figure()
ax95=Axis(f95[1,1])
lines!(ax95, kdam_range95, res_ss95)
lines!(ax95, reverse(kdam_range95), res1_ss95)

