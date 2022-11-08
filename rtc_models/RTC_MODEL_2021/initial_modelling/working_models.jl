using DifferentialEquations, Plots, DataFrames, BenchmarkTools, Measures

include("rtcR_models.jl");
include("params.jl");
include("functions.jl");

prob1 = ODEProblem(full_model!, init_full, tspan, params_full);
sol1 = solve(prob1, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

prob2 = ODEProblem(hill!, init_full, tspan, params_full);
sol2 = solve(prob2, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

prob3 = ODEProblem(mass_cons_a!, init2, tspan, params);
sol3 = solve(prob3, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

prob4 = ODEProblem(mass_cons_a_with_hill!, init2, tspan, params);
sol4 = solve(prob4, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

prob5 = ODEProblem(mass_cons_hill_mm!, init, tspan, params);
sol5 = solve(prob5, alg=Rodas4(), abstol=1e-6, reltol=1e-6);

prob6 = ODEProblem(simple_odes_mm_e_mass_cons_a!, init, tspan, params);
sol6 = solve(prob6, alg=Rodas4(), abstol=1e-6, reltol=1e-6);


p1 = plot(sol1, labels=species_full, title="1-- simple odes",legendfontsize=16,xtickfontsize=14,ytickfontsize=14,titlefontsize=18); #, ylim=(0,12));
p2 = plot(sol2, labels=species_full, title="2-- 1 + hill",legendfontsize=16,xtickfontsize=14,ytickfontsize=14,titlefontsize=18); #, ylim=(0,12));
p3 = plot(sol3, labels=species3, title="3-- 1 + mass conservation for a",legendfontsize=16,xtickfontsize=14,ytickfontsize=14,titlefontsize=18);
p4 = plot(sol4, labels=species3, title="4-- 2 + 3",legendfontsize=16,xtickfontsize=14,ytickfontsize=14,titlefontsize=18);
p5 = plot(sol5, labels=species, title="with Hill function for formation of c/Ra"#="5-- 4 + mm for e"=#,legendfontsize=14,xtickfontsize=14,ytickfontsize=14,titlefontsize=18);
p6 = plot(sol6, labels=species, title="without Hill function"#="6-- mass cons, mm, no hill"=#,legendfontsize=14,xtickfontsize=14,ytickfontsize=14,titlefontsize=18);

# plotlyjs()

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(2200,2000), margin=5mm)

plot(p6, p5, layout=(2,1), size=(1200,1000),margin=5mm)

a, b, c, d, e, p = separate_species2!(species_full2, sol1);
a_tot_1 = a+c+e;
c_tot_1 = c+e;
p1_mc = plot([sol1.t,sol1.t], [a_tot_1,c_tot_1], title="1 - simple odes");

a2, b2, c2, d2, e2, p2 = separate_species2!(species_full2, sol2);
a_tot_2 = a2+c2+e2;
c_tot_2 = c2+e2;
p2_mc = plot([sol2.t, sol2.t], [a_tot_2, c_tot_2], title="2 - 1 + hill");

b3, c3, d3, e3, p3 = separate_species3!(species4, sol3);
a_tot_3 = (b3/6)+c3+e3;
c_tot_3 = c3+e3;
p3_mc = plot([sol3.t, sol3.t], [a_tot_3, c_tot_3], title="3 - 1 + mass conservation for a");

b4, c4, d4, e4, p4 = separate_species3!(species4, sol4);
a_tot_4 = (b4/6)+c4+e4;
c_tot_4 = c4+e4;
p4_mc = plot([sol4.t, sol4.t], [a_tot_4, c_tot_4], title="4 - 3 + 2");

b5, c5, d5, p5 = separate_species!(species2, sol5);
a_tot_5 = (b5/6)+c5+((k3*d5.*c5./(k4.+katp.+k3.*d5)));
c_tot_5 = c5+((k3*d5.*c5./(k4.+katp.+k3.*d5)));
p5_mc = plot([sol5.t, sol5.t], [a_tot_5, c_tot_5], title="with Hill function for formation of c/Ra"#="5-- 4 + mm for e"=#);

b6, c6, d6, p6 = separate_species!(species2, sol6);
a_tot_6 = (b6/6)+c6+((k3*d6.*c6./(k4.+katp.+(k3.*d6))));
c_tot_6 = c6+((k3*d6.*c6./(k4.+katp.+(k3.*d6))));
p6_mc = plot([sol6.t, sol6.t], [a_tot_6, c_tot_6], title="without Hill function"#="6-- mass cons, mm, no hill"=#);

plot(p1_mc, p2_mc, p3_mc, p4_mc, p5_mc, p6_mc, layout=(6,1), size=(600,1200), labels=["a tot" "c tot"], margin=5mm)

plot(p6_mc, p5_mc, layout=(2,1), size=(1200,1000), labels=["a tot" "c tot"], margin=5mm)

plot(c5, p5)