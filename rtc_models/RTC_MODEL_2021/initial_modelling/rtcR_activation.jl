using DifferentialEquations, Plots, DataFrames, BenchmarkTools, Measures

include("rtcR_models.jl");
include("params.jl");
include("functions.jl");

prob = ODEProblem(mass_cons_hill_mm!, init, tspan, params);
sol = solve(prob, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

psol = plot(sol, labels=species, title="new params", titlefontsize=10);# title="RtcR activation/open complex formation"); #,legendfontsize=14,xtickfontsize=14,ytickfontsize=14,titlefontsize=18);

plot(psol)

# plotlyjs()

b, c, d, p = separate_species!(species2, sol);
a_tot = (b/6)+c+((k3*d.*c./(k4.+katp.+k3.*d)));
c_tot = c+((k3*d.*c./(k4.+katp.+k3.*d)));
mass_cons = plot([sol.t, sol.t], [a_tot, c_tot], title="Mass conservation in RtcR activation/open complex formation",labels=["atot" "ctot"]);

plot(mass_cons)


