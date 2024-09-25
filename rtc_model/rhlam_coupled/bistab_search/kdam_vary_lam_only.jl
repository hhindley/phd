using JLD2, Dates

model = "lam_coupled"
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))

lam_c_vals = 10 .^range(log10(1e-8),log10(0.1), length=100) #8e-7
wab_vals = 10 .^range(log10(1e-6),log10(0.1), length=100) #1.3e-5

lamkin_ssvals1 = []
lamkin_vals1 = []
# times_df=[]
global iteration_num = 0
println("starting loop")
for i in lam_c_vals
    for wab in wab_vals
        global iteration_num += 1
        println("iteration: $iteration_num")
        test_params = deepcopy(params_rtc1)
        test_params[lam_c] = i
        test_params[ω_ab] = wab
        test_params[kdam] = 0
        push!(lamkin_vals1, (i,wab))
        solu = sol(lam_coupled, init_rtc, tspan, test_params)
        df = create_solu_df(solu, species_rtc)
        push!(lamkin_ssvals1, df.rtca)
        # push!(times_df, df.time)
    end
end

good_range = []
for i in eachindex(lamkin_ssvals1)
    if lamkin_ssvals1[i][end] < 0.5 && lamkin_ssvals1[i][end] > 0.001
        push!(good_range, i)
    end
end

# f = Figure()
# ax = Axis(f[1, 1], xscale=log10)
# [lines!(ax, times_df[i], lamkin_ssvals1[i]) for i in (good_range)]

lams=[lamkin_vals1[good_range][i][1] for i in eachindex(lamkin_vals1[good_range])]
wabs=[lamkin_vals1[good_range][i][2] for i in eachindex(lamkin_vals1[good_range])]

lam_c_range = range(minimum(lams), maximum(lams), length=50)
wab_range = range(minimum(wabs), maximum(wabs), length=50)
kdam_range = range(0,200,length=200)

lam_ssvals = []
lam_vals = []

start_time = Dates.now()
println("starting at $start_time")
global iteration_num = 0

println("starting loop")
for i in lam_c_range
    for wab in wab_range
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

println("saving variables")
@save joinpath(homedir(), "phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search_just_lam.jld2") lam_ssvals lam_vals
end_time = Dates.now()
println("finished at $end_time")

duration = (end_time - start_time) / Hour(1)
println("Duration: $(duration) hours")


@load joinpath(homedir(),"phd/rtc_model/rhlam_coupled/bistab_search/kdam_vary_search_just_lam.jld2") lam_ssvals lam_vals


# plot([scatter(x=kdam_range, y=lamkin_ssvals[i], name="lam_c, kin_c, wab = $(lamkin_vals[i])") for i in eachindex(lamkin_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))



using GLMakie, InteractiveViz
function bistab_check(ssvals)
    m, m_ind = findmax(ssvals)
    for i in ssvals[m_ind:end]
        if ((i-m)/m*100) < -50
            return true
        end
    end
    return false
end

indices=[]
for i in eachindex(lam_ssvals)
    if bistab_check(lam_ssvals[i])
        push!(indices, i)
    end
end
indices

f = Figure()
ax = Axis(f[1, 1])
[ilines!(ax, kdam_range, lam_ssvals[i]) for i in indices]


f = Figure()
ax = Axis(f[1, 1])
[ilines!(ax, kdam_range, lam_ssvals[i]) for i in eachindex(lam_ssvals)]

sm = []
for i in eachindex(lam_ssvals)
    if maximum(lam_ssvals[i]) > 0.001 && maximum(lam_ssvals[i]) < 10
        push!(sm, lam_ssvals[i])
    end
end


f = Figure()
ax = Axis(f[1, 1])
[ilines!(ax, kdam_range, sm[i]) for i in eachindex(sm)]
