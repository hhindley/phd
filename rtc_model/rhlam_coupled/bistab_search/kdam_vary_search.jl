using JLD2, Dates
println("starting at $(Dates.now())")
# using GLMakie, InteractiveViz
model = "lamkin_coupled"
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))
lam_c_vals = 10 .^range(log10(1e-8),log10(0.1), length=10) #8e-7
kin_c_vals = 10 .^range(log10(1e-6),log10(0.1), length=10) #1.5e-5
wab_vals = 10 .^range(log10(1e-6),log10(0.1), length=10) #1.3e-5
kdam_range = range(0,100,length=20)

lamkin_ssvals = []
lamkin_vals = []

global iteration_num = 0
println("starting loop")
for i in lam_c_vals
    for j in kin_c_vals
        for wab in wab_vals
            global iteration_num += 1
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
# println("finished at $(Dates.now())")

@load joinpath(homedir(),"phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search.jld2") lamkin_ssvals lamkin_vals


# # plot([scatter(x=kdam_range, y=lamkin_ssvals[i], name="lam_c, kin_c, wab = $(lamkin_vals[i])") for i in eachindex(lamkin_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))




function bistab_check(ssvals)
    m, m_ind = findmax(ssvals)
    for i in ssvals[m_ind:end]
        if ((i-m)/m*100) < -50 && m > 0.005
            return true
        end
    end
    return false
end
indices=[]
for i in eachindex(lamkin_ssvals)
    if bistab_check(lamkin_ssvals[i])
        push!(indices, i)
    end
end
indices

f = Figure()
ax = Axis(f[1, 1])
[ilines!(ax, kdam_range, lamkin_ssvals[i]) for i in indices]

lines(kdam_range, lamkin_ssvals[indices[19]])

kdam_range = range(0,100,length=100)
new_res = []
global iteration_num = 0
for i in indices
    global iteration_num += 1
    println("iteration: $iteration_num")
    test_params = deepcopy(params_rtc1)
    test_params[lam_c] = lamkin_vals[i][1]
    test_params[kin_c] = lamkin_vals[i][2]
    test_params[ω_ab] = lamkin_vals[i][3]

    ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
    res = var_param(lamkin_coupled, kdam, test_params, kdam_range, ssvals_rtc_test)
    push!(new_res, res.rtca)
end

for i in eachindex(new_res)
    display(lines(kdam_range, new_res[i]))
    sleep(2)
end
lines(kdam_range, new_res[20])

kdam_range = range(0,200,length=100)
new_res_ss = []
global iteration_num = 0
for i in indices
    global iteration_num += 1
    println("iteration: $iteration_num")
    test_params = deepcopy(params_rtc1)
    test_params[lam_c] = lamkin_vals[i][1]
    test_params[kin_c] = lamkin_vals[i][2]
    test_params[ω_ab] = lamkin_vals[i][3]

    ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
    res_ss = numerical_bistability_analysis(lamkin_coupled, test_params, init_rtc, species_rtc, kdam_range, kdam, specie=:rtca)
    res1_ss = numerical_bistability_analysis(lamkin_coupled, test_params, ssvals_rtc_test, species_rtc, reverse(kdam_range), kdam, specie=:rtca)

    push!(new_res_ss, (res_ss, res1_ss))
end
x=1
lines(kdam_range, new_res_ss[x][1])
lines!(reverse(kdam_range), new_res_ss[x][2])

lams = [lamkin_vals[indices][i][1] for i in eachindex(indices)]
kins = [lamkin_vals[indices][i][2] for i in eachindex(indices)]
wabs = [lamkin_vals[indices][i][3] for i in eachindex(indices)]

f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, 1, 1e-8, color=:black)
scatter!(ax, 1, 0.1, color=:black)
scatter!(ax, fill(1, length(lams)), lams)

f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, 1, 1e-6, color=:black)
scatter!(ax, 1, 0.1, color=:black)
scatter!(ax, fill(1, length(kins)), kins)

f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, 1, 1e-6, color=:black)
scatter!(ax, 1, 0.1, color=:black)
scatter!(ax, fill(1, length(wabs)), wabs)

lamkin_vals[indices][20]
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lamkin_vals[indices][20][1]
test_params[kin_c] = lamkin_vals[indices][20][2]
test_params[ω_ab] = lamkin_vals[indices][20][3]
kdam_range = range(0,500,length=100)
ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
res = var_param(lamkin_coupled, kdam, test_params, kdam_range, ssvals_rtc_test)

lines(kdam_range, res.rh)

res_ss = numerical_bistability_analysis(lamkin_coupled, test_params, init_rtc, species_rtc, kdam_range, kdam)
res1_ss = numerical_bistability_analysis(lamkin_coupled, test_params, ssvals_rtc_test, species_rtc, reverse(kdam_range), kdam)
lines(kdam_range, res_ss.rtca)
lines!(reverse(kdam_range), res1_ss.rtca)


f = Figure()
ax = Axis(f[1, 1])
lines!(ax, kdam_range, lamkin_ssvals[indices[3]])

lamkin_vals[indices[5]]

lam_c_val = lamkin_vals[indices[3]][1]
kin_c_val = lamkin_vals[indices[3]][2]
ω_ab_val = lamkin_vals[indices[3]][3]

kdam_range1 = range(0,200,length=1000)

test_params = deepcopy(params_rtc1)
test_params[lam_c] = lam_c_val
test_params[kin_c] = kin_c_val
test_params[ω_ab] = ω_ab_val
ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
res = var_param(lamkin_coupled, kdam, test_params, kdam_range1, ssvals_rtc_test)

f = Figure()
ax = Axis(f[1, 1])
lines!(ax, kdam_range1, res.rtca)








lam_c_vals = 10 .^range(log10(1e-8),log10(0.1), length=10) #8e-7
kin_c_vals = 10 .^range(log10(1e-6),log10(0.1), length=10) #1.5e-5
wab_vals = 10 .^range(log10(1e-6),log10(0.1), length=10) #1.3e-5
kdam_range = range(0,100,length=20)

lamkin_ssvals1 = []
lamkin_vals1 = []
times_df=[]
global iteration_num = 0
println("starting loop")
for i in lam_c_vals
    for j in kin_c_vals
        for wab in wab_vals
            global iteration_num += 1
            println("iteration: $iteration_num")
            test_params = deepcopy(params_rtc1)
            test_params[lam_c] = i
            test_params[kin_c] = j
            test_params[ω_ab] = wab
            test_params[kdam] = 0
            push!(lamkin_vals1, (i,j,wab))
            solu = sol(lamkin_coupled, init_rtc, tspan, test_params)
            df = create_solu_df(solu, species_rtc)
            # ssvals_rtc_test = steady_states(lamkin_coupled, init_rtc, test_params)
            # res = var_param(lamkin_coupled, kdam, test_params, kdam_range, ssvals_rtc_test)
            push!(lamkin_ssvals1, df.rtca)
            push!(times_df, df.time)
        end
    end
end

times_df[1]
lamkin_ssvals1[1]

f = Figure()
ax = Axis(f[1, 1], xscale=log10)
[lines!(ax, times_df[i], lamkin_ssvals1[i]) for i in eachindex(lamkin_ssvals1)]

lamkin_ssvals1

good_range = []
for i in eachindex(lamkin_ssvals1)
    if lamkin_ssvals1[i][end] < 0.5 && lamkin_ssvals1[i][end] > 0.001
        push!(good_range, i)
    end
end
good_range

f = Figure()
ax = Axis(f[1, 1], xscale=log10)
[lines!(ax, times_df[i], lamkin_ssvals1[i]) for i in (good_range)]



lams=[lamkin_vals1[good_range][i][1] for i in eachindex(lamkin_vals1[good_range])]
minimum(lams)
maximum(lams)

kins=[lamkin_vals1[good_range][i][2] for i in eachindex(lamkin_vals1[good_range])]
minimum(kins)
maximum(kins)

wabs=[lamkin_vals1[good_range][i][3] for i in eachindex(lamkin_vals1[good_range])]
minimum(wabs)
maximum(wabs)


lam_c_range = range(minimum(lams), maximum(lams), length=20)
wab_range = range(minimum(wabs), maximum(wabs), length=20)
kdam_range = range(0,200,length=200)

lam_ssvals = []
lam_vals = []

global iteration_num = 0
println("starting loop")
for i in lam_c_vals
    for wab in wab_vals
        global iteration_num += 1
        println("iteration: $iteration_num")
        test_params = deepcopy(params_rtc1)
        test_params[lam_c] = i
        test_params[ω_ab] = wab
        push!(lam_vals, (i,wab))
        ssvals_rtc_test = steady_states(lam_coupled, init_rtc, test_params)
        res = var_param(lam_coupled, kdam, test_params, kdam_range, ssvals_rtc_test)
        push!(lam_ssvals, res.rtca)
    end
end