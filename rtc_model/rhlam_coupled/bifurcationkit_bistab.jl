include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))


br = get_br(test, ssvals_rtc, params_rtc1, 100.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
p_conc = plot([scatter(x=df.kdam, y=df.rtcb, name="RtcB", line=attr(color=:blue)), 
               scatter(x=df.kdam, y=df.rtca, name="RtcA", line=attr(color=:red))], 
               Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))

br = get_br(test, init_on, params_rtc1, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
p_conc = plot(scatter(x=df.kdam, y=df.rtca, name="RtcA"), 
        Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))

br = get_br(test, init_off, params_rtc1, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
p_conc = plot(scatter(x=df.kdam, y=df.rtca, name="RtcA"), 
        Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))
        





prob_bf = ODEProblem(lamkin_coupled, ssvals_rtc_test, tspan, test_params; jac=true)
odefun = prob_bf.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)
par_tm = prob_bf.p[1]
# id_kdam = indexof(kdam, parameters(model))
id_kdam = indexof(0.0, par_tm)
# Bifurcation Problem

prob = BifurcationProblem(F, prob_bf.u0, (par_tm), (@lens _[id_kdam]); J=J,
record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
opts_br = ContinuationPar(p_min = 0., p_max = 200., ds = 1e-4, #a=0.1,
dsmax = 0.0005, dsmin = 1e-6, # 0.15
# options to detect bifurcations
detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 50, #3,2,10
# number of eigenvalues
nev = 3, 
# maximum number of continuation steps
max_steps = 1000000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
# tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
# only using parameters that make a difference to solution
# continuation of equilibria
norminf(x) = norm(x, Inf)


br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)

df1 = create_br_df(br)
bf = bf_point_df(br)
kdam1 = findall(x->x==bf.kdam[1],df1.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df1.kdam)[1]
kdam3 = findall(x->x==bf.kdam[3],df1.kdam)[1]
kdam4 = findall(x->x==bf.kdam[4],df1.kdam)[1]
kdam5 = findall(x->x==bf.kdam[5],df1.kdam)[1]
kdam6 = findall(x->x==bf.kdam[6],df1.kdam)[1]

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="rtca")
lines!(ax, df1.kdam, df1.rtca)
scatter!(ax, [df1.kdam[kdam1]], [df1.rtca[kdam1]])
scatter!(ax, [df1.kdam[kdam2]], [df1.rtca[kdam2]])
scatter!(ax, [df1.kdam[kdam3]], [df1.rtca[kdam3]])
scatter!(ax, [df1.kdam[kdam4]], [df1.rtca[kdam4]])
scatter!(ax, [df1.kdam[kdam5]], [df1.rtca[kdam5]])
scatter!(ax, [df1.kdam[kdam6]], [df1.rtca[kdam6]])

# plot(scatter(x=df1.kdam, y=df1.rtca, name="RtcA"), 
#       Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))