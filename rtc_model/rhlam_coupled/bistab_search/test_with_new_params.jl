include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))
model = "lamkin_coupled"
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))
lam_c_val = 2.15e-6
kin_c_val = 1e-6
ω_ab_val = 1.3e-5
test_params = deepcopy(params_rtc1)
test_params[lam_c] = lam_c_val
test_params[kin_c] = kin_c_val
test_params[ω_ab] = ω_ab_val
ssvals_rtc_test = steady_states(model, init_rtc, test_params)

kdam_range = range(0,100,length=50)
res = var_param(model, kdam, test_params, kdam_range, ssvals_rtc_test)

lines(kdam_range, res.rtca)

res_n = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam)
res_n1 = numerical_bistability_analysis(model, test_params, ssvals_rtc_test, species_rtc, reverse(kdam_range), kdam)

lines(kdam_range, res_n.rtca)
lines!(reverse(kdam_range), res_n1.rtca)


br = get_br(model, ssvals_rtc_test, test_params; 
            kdam_max=200., 
            ds=1e-4, dsmax=0.0005, dsmin=1e-8, 
            detect_bifurcation=3, 
            n_inversion=2, 
            max_bisection_steps=50, 
            nev=3, 
            max_steps=1000000, 
            θ=0.01
)

df1 = create_br_df(br)
bf = bf_point_df(br)
# kdam1 = findall(x->x==bf.kdam[1],df1.kdam)[1]
# kdam2 = findall(x->x==bf.kdam[2],df1.kdam)[1]
# kdam3 = findall(x->x==bf.kdam[3],df1.kdam)[1]
# kdam4 = findall(x->x==bf.kdam[4],df1.kdam)[1]
# kdam5 = findall(x->x==bf.kdam[5],df1.kdam)[1]
# kdam6 = findall(x->x==bf.kdam[6],df1.kdam)[1]

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="rtca")
lines!(ax, df1.kdam, df1.rtca)
scatter!(ax, [df1.kdam[kdam1]], [df1.rtca[kdam1]])
scatter!(ax, [df1.kdam[kdam2]], [df1.rtca[kdam2]])
scatter!(ax, [df1.kdam[kdam3]], [df1.rtca[kdam3]])
scatter!(ax, [df1.kdam[kdam4]], [df1.rtca[kdam4]])
scatter!(ax, [df1.kdam[kdam5]], [df1.rtca[kdam5]])
scatter!(ax, [df1.kdam[kdam6]], [df1.rtca[kdam6]])

lamkin_vals[indices]

new_params = []
new_ssvals = []
for i in lamkin_vals[indices]
    test_params = deepcopy(params_rtc1)
    test_params[lam_c] = i[1]
    test_params[kin_c] = i[2]
    test_params[ω_ab] = i[3]
    ssvals_rtc_test = steady_states(model, init_rtc, test_params)
    push!(new_params, test_params)
    push!(new_ssvals, ssvals_rtc_test)
end
lamkin_vals[indices][20]
x=20
br = get_br(model, new_ssvals[x], new_params[x]; 
            kdam_max=500., 
            ds=0.0001, dsmax=0.005, dsmin=0.000001, 
            detect_bifurcation=3, 
            n_inversion=2, 
            max_bisection_steps=20, 
            nev=2, 
            max_steps=1000000, 
            θ=0.5,
            tol_stability=1e-12,
            tol_bisection_eigenvalue=1e-20
)

df1 = create_br_df(br)
lines(df1.kdam, df1.rtca)
bf = bf_point_df(br)

kdam1 = findall(x->x==bf.kdam[1],df1.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df1.kdam)[1]

scatter!([df1.kdam[kdam1]], [df1.rtca[kdam1]])
scatter!([df1.kdam[kdam2]], [df1.rtca[kdam2]])















model = "lam_coupled" # don't see "bistability" in this version so there is a difference 
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))
lam_c_val = 2.15e-6
kin_c_val = 1e-6
ω_ab_val = 1.3e-5

test_params = deepcopy(params_rtc1)
test_params[lam_c] = lam_c_val
test_params[kin_c] = kin_c_val
test_params[ω_ab] = ω_ab_val

res = var_param(model, kdam, test_params, kdam_range, ssvals_rtc_test)

lines(kdam_range, res.rtca)

res_n = numerical_bistability_analysis(model, test_params, init_rtc, species_rtc, kdam_range, kdam)
res_n1 = numerical_bistability_analysis(model, test_params, ssvals_rtc_test, species_rtc, kdam_range1, kdam)

lines(kdam_range, res_n.rtca)
lines!(kdam_range, res_n1.rtca)