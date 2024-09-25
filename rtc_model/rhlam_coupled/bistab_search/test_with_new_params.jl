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


function bifurcation(mdoel, init_vals, param_vals; p_max=200., ds=1e-4, dsmax=0.0005, dsmin=1e-6, detect_bifurcation=3, n_inversion=2, max_bisection_steps=50, nev=3, max_steps=1000000, θ=0.5)
    prob_bf = ODEProblem(mdoel, init_vals, tspan, param_vals; jac=true)
    odefun = prob_bf.f
    F = (u,p) -> odefun(u,p,0)
    J = (u,p) -> odefun.jac(u,p,0)
    par_tm = prob_bf.p[1]
    # id_kdam = indexof(kdam, parameters(model))
    id_kdam = indexof(0.0, par_tm)
    # Bifurcation Problem

    prob = BifurcationProblem(F, prob_bf.u0, (par_tm), (@lens _[id_kdam]); J=J,
    record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(p_min = 0., p_max = p_max, ds = ds, #a=0.1,
    dsmax = dsmax, dsmin = dsmin, # 0.15
    # options to detect bifurcations
    detect_bifurcation = detect_bifurcation, n_inversion = n_inversion, max_bisection_steps = max_bisection_steps, #3,2,10
    # number of eigenvalues
    nev = nev, 
    # maximum number of continuation steps
    max_steps = max_steps,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
    # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
    # only using parameters that make a difference to solution
    # continuation of equilibria
    norminf(x) = norm(x, Inf)

    br = continuation(prob, PALC(θ=θ), opts_br; plot = false, bothside=true, normC = norminf)
    return br
end

br = get_br(model, ssvals_rtc_test, test_params; p_max=200., ds=1e-4, dsmax=0.0005, dsmin=1e-8, detect_bifurcation=3, n_inversion=2, max_bisection_steps=50, nev=3, max_steps=1000000, θ=0.01)

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