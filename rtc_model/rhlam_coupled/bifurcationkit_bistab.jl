include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

br = get_br(test, ssvals_rtc, params_rtc1, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
p_conc = plot([scatter(x=df.kdam, y=df.rtcb, name="RtcB", line=attr(color=:blue)), scatter(x=df.kdam, y=df.rtca, name="RtcA", line=attr(color=:red))], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))


prob = ODEProblem(test, ssvals_rtc, tspan, params_rtc1; jac=true)
odefun = prob.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)
par_tm = prob.p[1]
# id_kdam = indexof(kdam, parameters(model))
id_kdam = indexof(0.0, par_tm)
# Bifurcation Problem

prob = BifurcationProblem(F, prob.u0, (par_tm), (@lens _[id_kdam]); J=J,
record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
opts_br = ContinuationPar(p_min = 0., p_max = 1., ds = 0.001, #a=0.1,
dsmax = 0.15, dsmin = 0.0001, # 0.15
# options to detect bifurcations
detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 20, #3,2,10
# number of eigenvalues
nev = 2, 
# maximum number of continuation steps
max_steps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
# tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
# only using parameters that make a difference to solution
# continuation of equilibria
br = continuation(prob, PALC(θ=0.1), opts_br; plot = false, bothside=true, normC = norminf)

df1 = create_br_df(br)

plot([scatter(x=df1.kdam, y=df1.rtca, name="RtcA", line=attr(color=:red)), scatter(x=df1.kdam, y=df1.rtcb, name="RtcB", line=attr(color=:blue))], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))