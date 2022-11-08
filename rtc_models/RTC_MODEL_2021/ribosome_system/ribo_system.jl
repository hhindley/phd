using Plots, DifferentialEquations, DataFrames

include("/home/holliehindley/rtc_model/rtc_expression/MWC_functions.jl")
include("ribo_functions.jl")
include("ribo_params.jl")

# plotlyjs()

sol = simple_solve!(ribo_odes!, init, tspan, params)

# sol_alg = simple_solve!(ribo_odes_alg!, init_alg, tspan, params_alg)

sol_alg = simple_solve!(ribo_all_alg!, new_init, tspan, new_params)

sol_mix1 = simple_solve!(ribo_mix1!, init_mix1, tspan, params_mix1)

sol_mix2 = simple_solve!(ribo_mix2!, init_mix2, tspan, params_mix2)

p1 = plot(sol, labels=["RtcA" "RtcB" "RdRtcA" "RtRtcB" "Rd" "Rt" "Rh"], title="all ODEs");
p2 = plot(sol_alg, labels=["RtcA" "RtcB" "Rd" "Rh"], title="all algebraic");
p3 = plot(sol_mix1, labels=["RtcA" "RtcB" "Rd" "Rt" "Rh"], title="Rh ODE");
p4 = plot(sol_mix2, labels=["RtcA" "RtcB" "Rd" "Rt" "RdRtcA" "RtRtcB"], title="Ints ODE");

plot(p1, p2, p3, p4, size=(800,600))

rtca, rtcb, rdrtca, rtrtcb, rd, rt, rh = split(species_full, sol)
rtca1, rtcb1, rdrtca1, rtrtcb1, rd1, rt1, rh1 = split_all_alg(new_species, sol_alg)
rtca2, rtcb2, rdrtca2, rtrtcb2, rd2, rt2, rh2 = split_mix1(species_mix1, sol_mix1)
rtca3, rtcb3, rdrtca3, rtrtcb3, rd3, rt3, rh3 = split_mix2(species_mix2, sol_mix2)

p1a = plot(sol.t, [rtca, rtcb, rd, rh, rt], labels=["RtcA" "RtcB" "Rd" "Rh" "Rt"], title="all ODEs", xaxis=(:log10, (1,Inf)));
p2a = plot(sol_alg.t, [rtca1, rtcb1, rd1, rh1, rt1], labels=["RtcA" "RtcB" "Rd" "Rh" "Rt"], title="all algebraic", xaxis=(:log10, (1,Inf)));
p3a = plot(sol_mix1.t, [rtca2, rtcb2, rd2, rh2, rt2], labels=["RtcA" "RtcB" "Rd" "Rh" "Rt"], title="Rh ODE", xaxis=(:log10, (1,Inf)));
p4a = plot(sol_mix2.t, [rtca3, rtcb3, rd3, rh3, rt3], labels=["RtcA" "RtcB" "Rd" "Rh" "Rt"], title="Ints ODE", xaxis=(:log10, (1,Inf)));
plot(p1a, p2a, p3a, p4a, size=(800,600))

p1b = plot(sol.t, rt+rd+rh+rdrtca+rtrtcb, title="all ODEs");
p2b = plot(sol_alg.t, rt1+rd1+rh1+rdrtca1+rtrtcb1, title="all algebraic");
p3b = plot(sol_mix1.t, rt2+rd2+rh2+rdrtca2+rtrtcb2, title="Rh ODE");
p4b = plot(sol_mix2.t, rt3+rd3+rh3+rdrtca3+rtrtcb3, title="Ints ODE");
plot(p1b, p2b, p3b, p4b, size=(800,600))

a = plot(sol.t, rtca);
a1 = plot(sol_alg.t, rtca1);
a2 = plot(sol_mix1.t, rtca2);
a3 = plot(sol_mix2.t, rtca3);
plot(a, a1, a2, a3, size=(800,600))

b = plot(sol.t, rtcb);
b1 = plot(sol_alg.t, rtcb1);
b2 = plot(sol_mix1.t, rtcb2);
b3 = plot(sol_mix2.t, rtcb3);
plot(b, b1, b2, b3, size=(800,600))

c = plot(sol.t, rd);
c1 = plot(sol_alg.t, rd1);
c2 = plot(sol_mix1.t, rd2);
c3 = plot(sol_mix2.t, rd3);
plot(c, c1, c2, c3, size=(800,600))

d = plot(sol.t, rh);
d1 = plot(sol_alg.t, rh1);
d2 = plot(sol_mix1.t, rh2);
d3 = plot(sol_mix2.t, rh3);
plot(d, d1, d2, d3, size=(800,600))

e = plot(sol.t, rt);
e1 = plot(sol_alg.t, rt1);
e2 = plot(sol_mix1.t, rt2);
e3 = plot(sol_mix2.t, rt3);
plot(e, e1, e2, e3, size=(800,600))

f = plot(sol.t, rdrtca);
f1 = plot(sol_alg.t, rdrtca1);
f2 = plot(sol_mix1.t, rdrtca2);
f3 = plot(sol_mix2.t, rdrtca3);
plot(f, f1, f2, f3, size=(800,600))

g = plot(sol.t, rtrtcb);
g1 = plot(sol_alg.t, rtrtcb1);
g2 = plot(sol_mix1.t, rtrtcb2);
g3 = plot(sol_mix2.t, rtrtcb3);
plot(g, g1, g2, g3, size=(800,600))

h = plot(sol.t, rdrtca+rtrtcb);
h1 = plot(sol_alg.t, rdrtca1+rtrtcb1);
h2 = plot(sol_mix1.t, rdrtca2+rtrtcb2);
h3 = plot(sol_mix2.t, rdrtca3+rtrtcb3);
plot(h, h1, h2, h3, size=(800,600))

i = plot(sol.t, rt+rd+rh);
i1 = plot(sol_alg.t, rt1+rd1+rh1);
i2 = plot(sol_mix1.t, rt2+rd2+rh2);
i3 = plot(sol_mix2.t, rt3+rd3+rh3);
plot(i, i1, i2, i3, size=(800,600))

j = plot(sol.t, rtrtcb+rd+rh);
j1 = plot(sol_alg.t, rtrtcb1+rd1+rh1);
j2 = plot(sol_mix1.t, rtrtcb2+rd2+rh2);
j3 = plot(sol_mix2.t, rtrtcb3+rd3+rh3);
plot(j, j1, j2, j3, size=(800,600))

k = plot(sol.t, rdrtca+rtca);
k1 = plot(sol_alg.t, rdrtca1+rtca1);
k2 = plot(sol_mix1.t, rdrtca2+rtca2);
k3 = plot(sol_mix2.t, rdrtca3+rtca3);
plot(k, k1, k2, k3, size=(800,600))

l = plot(sol.t, rtrtcb+rtcb);
l1 = plot(sol_alg.t, rtrtcb1+rtcb1);
l2 = plot(sol_mix1.t, rtrtcb2+rtcb2);
l3 = plot(sol_mix2.t, rtrtcb3+rtcb3);
plot(l, l1, l2, l3, size=(800,600))

println(rd+rh+rt+rtrtcb+rdrtca)
println(rd1+rh1+rt1+rtrtcb1+rdrtca1)
println(rd2+rh2+rt2+rtrtcb2+rdrtca2)
println(rd3+rh3+rt3+rtrtcb3+rdrtca3)


plot(sol.t, [rh, rd, rt], labels=["Rh" "Rd" "Rt"], title="Ribosome full ODEs")
plot(sol_alg.t, [rh_alg, rd_alg, rt_alg], labels=["Rh" "Rd" "Rt"], title="Ribosome Rt algebraic")