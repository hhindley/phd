using Plots, DifferentialEquations, DataFrames, Statistics, Measures
include("fullmodel_functions.jl")
include("fullmodel_params.jl")
include("models_new.jl")
# plotlyjs()

# sol_alg = simple_solve!(full_alg!, init_alg, tspan, params_alg) # rh and ints as algebraic - all algebraic 
# p2_sol = plot(sol_alg[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_alg, palette=:tab10, title="2 - all algebraic (ints and rh)");

# sol_full = simple_solve!(full_ODEs!, init_full, tspan, params_full)  # all ribosomes as ODEs and ints as odes - all odes 
# p3_sol = plot(sol_full[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_full, palette=:tab10,  title="5 - new ODEs")#title="3 - all ODEs (ints and rh)")

sol_full = simple_solve!(new_ODEs!, init_new, tspan, params_new)  # all ribosomes as ODEs and ints as odes - all odes 
p3_sol = plot(sol_full[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, palette=:tab10,  title="5 - new ODEs")#title="3 - all ODEs (ints and rh)")




# sol_mix2 = simple_solve!(mix2!, init_mix2, tspan, params_mix2) # ints as odes and rh algebraic
# p4_sol = plot(sol_mix2[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_mix2, palette=:tab10, title="4 - ints as ODEs, rh algebraic");

# sol_mix1 = simple_solve!(mix1!, init_mix1, tspan, params_mix1) #all ribosomes as ODEs, ints algebraic 
# p1_sol = plot(sol_mix1[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_mix1, palette=:tab10, title="1 - all ribo ODEs, ints algebraic");

sol_rtc = simple_solve!(RTC_SYSTEM!, init_rtc, tspan, params_rtc) #all ribosomes as ODEs, ints algebraic 
p5_sol = plot(sol_rtc[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_rtc, palette=:tab10, title="final working model")

# plot(p5_sol, p2_sol, p3_sol, p4_sol, size=(1000,800))
plot(p5_sol, p3_sol, size=(800,400))


# splitting species and looking at solution results
# rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, ribo_d1, ribo_h1, ribo_t1, alpha1, fa1, ra1, v1, tscr1, tlr1, tscr_b1, tlr_b1, rtcr_tscr1, rtcr_tlr1, rdrtca1, rtrtcb1 = split(species_mix1, sol_mix1)
# p1 = plot(sol_mix1.t, [rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_d1, ribo_h1, ribo_t1], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="1 - all ribo ODEs, ints algebraic");

# p_rtot1 = plot(sol_mix1.t, [ribo_d1+ribo_h1+ribo_t1+rdrtca1+rtrtcb1], xaxis=(:log10, (1,Inf)),title="1 - ribo odes, ints alg", legend=false);
# p_rtc1 = plot(sol_mix1.t, xlabel="t", [rtca1, rtcb1, rtcr1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"], title="1 - ribo odes, ints alg");
# p_mrna1 = plot(sol_mix1.t, xlabel="t", [rm_a1, rm_b1, rm_r1], xaxis=(:log10, (1,Inf)), labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], title="1 - ribo odes, ints alg");
# p_ribo1 = plot(sol_mix1.t, xlabel="t", [ribo_d1, ribo_h1, ribo_t1], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"], title="1 - ribo odes, ints alg");
# p_int1 = plot(sol_mix1.t, xlabel="t", [rtrtcb1, rdrtca1], xaxis=(:log10, (1,Inf)),labels=["RtRtcB" "RdRtcA"], title="1 - ribo odes, ints alg");

# rm_a2, rtca2, rm_b2, rtcb2, rm_r2, rtcr2, ribo_d2, ribo_h2, ribo_t2, alpha2, fa2, ra2, v2, tscr2, tlr2, tscr_b2, tlr_b2, rtcr_tscr2, rtcr_tlr2, rdrtca2, rtrtcb2 = split_rh(species_alg, sol_alg)
# p2 = plot(sol_alg.t, [rm_a2, rtca2, rm_b2, rtcb2, rm_r2, rtcr2, rdrtca2, rtrtcb2, ribo_d2, ribo_h2, ribo_t2], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="2 - ints and ribo h as algebraic");

# p_rtot2 = plot(sol_alg.t, [ribo_d2+ribo_h2+ribo_t2+rdrtca2+rtrtcb2], xaxis=(:log10, (1,Inf)),title="2 - ints and ribo h as algebraic", legend=false);
# p_rtc2 = plot(sol_alg.t, xlabel="t", [rtca2, rtcb2, rtcr2], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"], title="2 - ints and ribo h as algebraic");
# p_mrna2 = plot(sol_alg.t, xlabel="t", [rm_a2, rm_b2, rm_r2], xaxis=(:log10, (1,Inf)), labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], title="2 - ints and ribo h as algebraic");
# p_ribo2 = plot(sol_alg.t, xlabel="t", [ribo_d2, ribo_h2, ribo_t2], xaxis=(:log10, (1,Inf)),labels=["Rd" "Rh" "Rt"], title="2 - ints and ribo h as algebraic");
# p_int2 = plot(sol_alg.t, xlabel="t", [rtrtcb2, rdrtca2], xaxis=(:log10, (1,Inf)),labels=["RtRtcB" "RdRtcA"], title="2 - ints and ribo h as algebraic");

rm_a3, rtca3, rm_b3, rtcb3, rm_r3, rtcr3, ribo_d3, ribo_h3, ribo_t3, alpha3, fa3, ra3, v3, tscr3, tlr3, tscr_b3, tlr_b3, rtcr_tscr3, rtcr_tlr3, rdrtca3, rtrtcb3 = split_riboint(species_full, sol_full)
p3 = plot(sol_full.t, [rm_a3, rtca3, rm_b3, rtcb3, rm_r3, rtcr3, rdrtca3, rtrtcb3, ribo_d3, ribo_h3, ribo_t3], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="3 - ribo odes and int odes")

p_rtot3 = plot(sol_full.t, [ribo_d3+ribo_h3+ribo_t3+rdrtca3+rtrtcb3], xaxis=(:log10, (1,Inf)),title="3 - ribo odes and int odes", legend=false);
p_rtc3 = plot(sol_full.t, xlabel="t", [rtca3, rtcb3, rtcr3], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"], title="3 - ribo odes and int odes");
p_mrna3 = plot(sol_full.t, xlabel="t", [rm_a3, rm_b3, rm_r3], xaxis=(:log10, (1,Inf)),labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], title="3 - ribo odes and int odes");
p_ribo3 = plot(sol_full.t, xlabel="t", [ribo_d3, ribo_h3, ribo_t3], xaxis=(:log10, (1,Inf)),labels=["Rd" "Rh" "Rt"], title="3 - ribo odes and int odes");
p_int3 = plot(sol_full.t, xlabel="t", [rtrtcb3, rdrtca3], xaxis=(:log10, (1,Inf)),labels=["RtRtcB" "RdRtcA"], title="3 - ribo odes and int odes");


# rm_a4, rtca4, rm_b4, rtcb4, rm_r4, rtcr4, ribo_d4, ribo_h4, ribo_t4, alpha4, fa4, ra4, v4, tscr4, tlr4, tscr_b4, tlr_b4, rtcr_tscr4, rtcr_tlr4, rdrtca4, rtrtcb4 = split_all(species_mix2, sol_mix2)
# p4 = plot(sol_mix2.t, [rm_a4, rtca4, rm_b4, rtcb4, rm_r4, rtcr4, rdrtca4, rtrtcb4, ribo_d4, ribo_h4, ribo_t4], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="4 - ribo h alg and int odes")

# p_rtot4 = plot(sol_mix2.t, [ribo_d4+ribo_h4+ribo_t4+rdrtca4+rtrtcb4], xaxis=(:log10, (1,Inf)),title="4 - ribo h alg and int odes", legend=false)
# p_rtc4 = plot(sol_mix2.t, xlabel="t", [rtca4, rtcb4, rtcr4], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"], title="4 - ribo h alg and int odes")
# p_mrna4 = plot(sol_mix2.t, xlabel="t", [rm_a4, rm_b4, rm_r4], xaxis=(:log10, (1,Inf)),labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], title="4 - ribo h alg and int odes")
# p_ribo4 = plot(sol_mix2.t, xlabel="t", [ribo_d4, ribo_h4, ribo_t4], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"], title="4 - ribo h alg and int odes")
# p_int4 = plot(sol_mix2.t, xlabel="t", [rtrtcb4, rdrtca4], xaxis=(:log10, (1,Inf)), labels=["RtRtcB" "RdRtcA"], title="4 - ribo h alg and int odes")

rm_a5, rtca5, rm_b5, rtcb5, rm_r5, rtcr5, ribo_d5, ribo_h5, ribo_t5, alpha5, fa5, ra5, v5, tscr5, tlr5, tscr_b5, tlr_b5, rtcr_tscr5, rtcr_tlr5, rdrtca5, rtrtcb5 = split_rh(species_rtc, sol_rtc)
p5 = plot(sol_rtc.t, [rm_a5, rtca5, rm_b5, rtcb5, rm_r5, rtcr5, rdrtca5, rtrtcb5, ribo_d5, ribo_h5, ribo_t5], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="5 - new model")

p_rtot5 = plot(sol_rtc.t, [ribo_d5+ribo_h5+ribo_t5+rdrtca5+rtrtcb5], xaxis=(:log10, (1,Inf)),title="5 - new model", legend=false);
p_rtc5 = plot(sol_rtc.t, xlabel="t", [rtca5, rtcb5, rtcr5], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"], title="5 - new model");
p_mrna5 = plot(sol_rtc.t, xlabel="t", [rm_a5, rm_b5, rm_r5], xaxis=(:log10, (1,Inf)), labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], title="5 - new model");
p_ribo5 = plot(sol_rtc.t, xlabel="t", [ribo_d5, ribo_h5, ribo_t5], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"], title="5 - new model");
p_int5 = plot(sol_rtc.t, xlabel="t", [rtrtcb5, rdrtca5], xaxis=(:log10, (1,Inf)),labels=["RtRtcB" "RdRtcA"], title="5 - new model");

plot(sol_rtc.t, [rdrtca5, rtrtcb5], xaxis=(:log10, (1,Inf)), labels=["RdRtcA" "RtRtcB"])
print(maximum(rdrtca5))
# plot(p5, p2, p3, p4, layout=(2,2), size=(1000,800))
# plot(p_rtc5, p_rtc2, p_rtc3, p_rtc4, layout=(2,2), size=(1000,800))
# plot(p_mrna5, p_mrna2, p_mrna3, p_mrna4, layout=(2,2), size=(1000,800))
# plot(p_ribo5, p_ribo2, p_ribo3, p_ribo4, layout=(2,2), size=(1000,800))
# plot(p_int5, p_int2, p_int3, p_int4, layout=(2,2), size=(1000,800))
# plot(p_rtot5, p_rtot2, p_rtot3, p_rtot4, layout=(2,2), size=(1000,800), plot_title="rtot")

plot(p5, p3, size=(800,400))
plot(p_rtc5, p_rtc3, size=(800,400))
plot(p_mrna5, p_mrna3, size=(800,400))
plot(p_ribo5, p_ribo3, size=(800,400))
plot(p_int5, p_int3, size=(800,400))
plot(p_rtot5, p_rtot3, size=(800,400), plot_title="rtot")

plot(sol_full.t, [ribo_d3+ribo_h3+ribo_t3], xaxis=(:log10, (1,Inf)),title="3 - ribo odes and int odes", labels="Rd+Rt+Rh")

# p1_i = plot(sol_mix1.t, [rtca1+rdrtca1, rtcb1+rtrtcb1], title="1 - ribo odes, ints alg", labels=["RtcA" "RtcB"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)));
# p2_i = plot(sol_alg.t, [rtca2+rdrtca2, rtcb2+rtrtcb2], title="2 - ints and ribo h as algebraic", labels=["RtcA" "RtcB"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)));
p3_i = plot(sol_full.t, [rtca3+rdrtca3, rtcb3+rtrtcb3], title="3 - ribo odes and int odes", labels=["RtcA" "RtcB"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)));
# p4_i = plot(sol_mix2.t, [rtca4+rdrtca4, rtcb4+rtrtcb4], title="4 - ribo h alg and int odes", labels=["RtcA" "RtcB"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)));
p5_i = plot(sol_rtc.t, [rtca5+rdrtca5, rtcb5+rtrtcb5], title="5 - new model", labels=["RtcA" "RtcB"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)));

# plot(p5_i, p2_i, p3_i, p4_i, size=(1000,800), plot_title="RtcA and RtcB total")
plot(p5_i, p3_i, size=(800,400), plot_title="RtcA and RtcB total")

function plot_comp(species1, species2, species3, species4, species)
    plot(sol_rtc.t, species1, xaxis=(:log10, (1,Inf)), labels="5 - new model", title="$species")
    plot!(sol_alg.t, species2, xaxis=(:log10, (1,Inf)), labels="2 - ints and ribo h as algebraic")
    plot!(sol_full.t, species3, xaxis=(:log10, (1,Inf)), labels="3 - ribo odes and int odes")
    plot!(sol_mix2.t, species4, xaxis=(:log10, (1,Inf)), labels="4 - ribo h alg and int odes")
end

# plot_comp(ribo_h5, ribo_h2, ribo_h3, ribo_h4, "ribo_h")
# plot_comp(ribo_t5, ribo_t2, ribo_t3, ribo_t4, "ribo_t")
# plot_comp(ribo_d5, ribo_d2, ribo_d3, ribo_d4, "ribo_d")
# plot_comp(rtca5, rtca2, rtca3, rtca4, "RtcA")
# plot_comp(rtcb5, rtcb2, rtcb3, rtcb4, "RtcB")
# plot_comp(rtcr5, rtcr2, rtcr3, rtcr4, "RtcR")
# plot_comp(rdrtca5, rdrtca2, rdrtca3, rdrtca4, "RdRtcA")
# plot_comp(rtrtcb5, rtrtcb2, rtrtcb3, rtrtcb4, "RtRtcB")
# plot_comp(rm_a5, rm_a2, rm_a3, rm_a4, "RtcA mRNA")
# plot_comp(rm_b5, rm_b2, rm_b3, rm_b4, "RtcB mRNA")
# plot_comp(rm_r5, rm_r2, rm_r3, rm_r4, "RtcR mRNA")
# plot_comp(alpha5, alpha2, alpha3, alpha4, "alpha")
# plot_comp(fa5, fa2, fa3, fa4, "fa")
# plot_comp(ra5, ra2, ra3, ra4, "Ra")
# plot_comp(v5, v2, v3, v4, "v")
# plot_comp(tscr5, tscr2, tscr3, tscr4, "tscr_a")
# plot_comp(tscr_b5, tscr_b2, tscr_b3, tscr_b4, "tscr_b")
# plot_comp(tlr5, tlr2, tlr3, tlr4, "tlr_a")
# plot_comp(tlr_b5, tlr_b2, tlr_b3, tlr_b4, "tlr_b")
# plot_comp(rtcr_tlr5, rtcr_tlr2, rtcr_tlr3, rtcr_tlr4, "rtcr_tlr")

function plot_comp2(species1, species3, species)
    plot(sol_rtc.t, species1, xaxis=(:log10, (1,Inf)), labels="5 - new model", title="$species")
    plot!(sol_full.t, species3, xaxis=(:log10, (1,Inf)), labels="3 - ribo odes and int odes")
end

plot_comp2(ribo_h5, ribo_h3, "ribo_h")
plot_comp2(ribo_t5, ribo_t3, "ribo_t")
plot_comp2(ribo_d5, ribo_d3, "ribo_d")
plot_comp2(rtca5, rtca3, "RtcA")
plot_comp2(rtcb5, rtcb3, "RtcB")
plot_comp2(rtcr5, rtcr3, "RtcR")
plot_comp2(rdrtca5, rdrtca3, "RdRtcA")
plot_comp2(rtrtcb5, rtrtcb3, "RtRtcB")
plot_comp2(rm_a5, rm_a3, "RtcA mRNA")
plot_comp2(rm_b5, rm_b3, "RtcB mRNA")
plot_comp2(rm_r5, rm_r3, "RtcR mRNA")
plot_comp2(alpha5, alpha3, "alpha")
plot_comp2(fa5, fa3, "fa")
plot_comp2(ra5, ra3, "Ra")
plot_comp2(v5, v3, "v")
plot_comp2(tscr5, tscr3, "tscr_a")
plot_comp2(tscr_b5, tscr_b3, "tscr_b")
plot_comp2(tlr5, tlr3, "tlr_a")
plot_comp2(tlr_b5, tlr_b3, "tlr_b")
plot_comp2(rtcr_tlr5, rtcr_tlr3, "rtcr_tlr")







rd_range = range(0, 10, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in full model to ode and see if that works (probs will) 

params_full = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k]
init_full = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_h_0, ribo_d_0, ribo_t_0];
plot1_rd_ribo, plot2_rd_ribo = logplots_for_param("Rd", init_full, 10, rd_range, full_ODEs!, init_full, tspan, params_full, species_full, labels_full)
plot1_rd_ribo
plot2_rd_ribo # all ribosomes as algebraic


function empty_res()
    sols=[]
    sols_t=[]
    rm_a=[]
    rtca=[]
    rm_b=[]
    rtcb=[]
    rm_r=[]
    rtcr=[]
    ribo_d=[]
    ribo_h=[]
    ribo_t=[]
    rdrtca=[]
    rtrtcb=[]
    return sols, sols_t, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, rdrtca, rtrtcb
end

sols, sols_t, trm_a, trtca, trm_b, trtcb, trm_r, trtcr, tribo_d, tribo_h, tribo_t, trdrtca, trtrtcb = empty_res()
for i in rd_range
    init_full[10] = i
    # ribo_d_0 = i
    sol = simple_solve!(full_ODEs!, init_full, tspan, params_full)
    push!(sols, sol)
    push!(sols_t, sol.t)
    solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species_full)
    push!(trm_a, solDF[:, :rm_a])
    push!(trtca, solDF[:, :rtca])
    push!(trm_b, solDF[:, :rm_b])
    push!(trtcb, solDF[:, :rtcb])
    push!(trm_r, solDF[:, :rm_r])
    push!(trtcr, solDF[:, :rtcr])
    push!(tribo_d, solDF[:, :ribo_d])
    push!(tribo_h, solDF[:, :ribo_h])
    push!(tribo_t, solDF[:, :ribo_t])
    push!(trdrtca, solDF[:, :rdrtca])
    push!(trtrtcb, solDF[:, :rtrtcb])
end

p1 = plot(sols[1][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 0");
p2 = plot(sols[2][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 4");
p3 = plot(sols[3][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 8");
p4 = plot(sols[4][2:end], yaxis=(:log10, (1,Inf)),xaxis=(:log10, (1,Inf)), title="Rd = 12");
p5 = plot(sols[5][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 16");
p6 = plot(sols[6][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 20");
plot(p1, p2, p3, p4, p5, p6, size=(800,600), labels=labels_full, plot_title=("Full ODEs"))

t1 = plot(sols_t, trtcr, xaxis=(:log10, (1,Inf)), title="RtcR");
t2 = plot(sols_t, trtca, xaxis=(:log10, (1,Inf)), title="RtcA");
t3 = plot(sols_t, trtcb, xaxis=(:log10, (1,Inf)), title="RtcB");
t4 = plot(sols_t, tribo_d, xaxis=(:log10, (1,Inf)), title="Rd");
t5 = plot(sols_t, tribo_t, xaxis=(:log10, (1,Inf)), title="Rt");
t6 = plot(sols_t, tribo_h, xaxis=(:log10, (1,Inf)), title="Rh");
plot(t1, t2, t3, t4, t5, t6, size=(800,600), labels=["0" "4" "8" "12" "16" "20"], plot_title="Effect of Rd")


params_mix1 = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot]
init_mix1 = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_h_0, ribo_d_0, ribo_t_0];

sols1, sols_t1, trm_a1, trtca1, trm_b1, trtcb1, trm_r1, trtcr1, tribo_d1, tribo_h1, tribo_t1, trdrtca1, trtrtcb1 = empty_res()

for i in rd_range
    init_mix1[8] = i
    sol1 = simple_solve!(mix1!, init_mix1, tspan, params_mix1)
    push!(sols1, sol1)
    push!(sols_t1, sol1.t)
    solDF1 = DataFrame([[j[i] for j in sol1.u] for i=1:length(sol1.u[1])], species_mix1)
    push!(trm_a1, solDF1[:, :rm_a])
    push!(trtca1, solDF1[:, :rtca])
    push!(trm_b1, solDF1[:, :rm_b])
    push!(trtcb1, solDF1[:, :rtcb])
    push!(trm_r1, solDF1[:, :rm_r])
    push!(trtcr1, solDF1[:, :rtcr])
    push!(tribo_d1, solDF1[:, :ribo_d])
    push!(tribo_h1, solDF1[:, :ribo_h])
    push!(tribo_t1, solDF1[:, :ribo_t])
    push!(trtrtcb1, ka_b.*solDF1[:, :ribo_t].*rtcb_tot./(ka_b.*solDF1[:, :ribo_t].+kb_b.+kc_b.*atp))
    push!(trdrtca1, k1_a.*solDF1[:, :ribo_d].*rtca_tot./(k1_a.*solDF1[:, :ribo_d].+k2_a.+k3_a.*atp))
end

p11 = plot(sols1[1][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 0");
p21 = plot(sols1[2][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 4");
p31 = plot(sols1[3][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 8");
p41 = plot(sols1[4][2:end], yaxis=(:log10, (1,Inf)),xaxis=(:log10, (1,Inf)), title="Rd = 12");
p51 = plot(sols1[5][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 16");
p61 = plot(sols1[6][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 20");
plot(p11, p21, p31, p41, p51, p61, size=(800,600), labels=labels_mix1, plot_title="Intermediates as algebraic")

t11 = plot(sols_t1, trtcr1, xaxis=(:log10, (1,Inf)), title="RtcR");
t21 = plot(sols_t1, trtca1, xaxis=(:log10, (1,Inf)), title="RtcA");
t31 = plot(sols_t1, trtcb1, xaxis=(:log10, (1,Inf)), title="RtcB");
t41 = plot(sols_t1, tribo_d1, xaxis=(:log10, (1,Inf)), title="Rd");
t51 = plot(sols_t1, tribo_t1, xaxis=(:log10, (1,Inf)), title="Rt");
t61 = plot(sols_t1, tribo_h1, xaxis=(:log10, (1,Inf)), title="Rh");
plot(t11, t21, t31, t41, t51, t61, size=(800,600), labels=["0" "4" "8" "12" "16" "20"], plot_title="Effect of Rd")




params_mix2 = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot]
init_mix2 = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_d_0, ribo_t_0];

sols2, sols_t2, trm_a2, trtca2, trm_b2, trtcb2, trm_r2, trtcr2, tribo_d2, tribo_h2, tribo_t2, trdrtca2, trtrtcb2 = empty_res()

for i in rd_range
    init_mix2[9] = i
    println(init_mix2[9])

    sol2 = simple_solve!(mix2!, init_mix2, tspan, params_mix2)
    push!(sols2, sol2)
    push!(sols_t2, sol2.t)
    solDF2 = DataFrame([[j[i] for j in sol2.u] for i=1:length(sol2.u[1])], species_mix2)
    push!(trm_a2, solDF2[:, :rm_a])
    push!(trtca2, solDF2[:, :rtca])
    push!(trm_b2, solDF2[:, :rm_b])
    push!(trtcb2, solDF2[:, :rtcb])
    push!(trm_r2, solDF2[:, :rm_r])
    push!(trtcr2, solDF2[:, :rtcr])
    push!(tribo_d2, solDF2[:, :ribo_d])
    push!(tribo_t2, solDF2[:, :ribo_t])
    push!(trdrtca2, solDF2[:, :rtrtcb])
    push!(trtrtcb2, solDF2[:, :rdrtca])
    push!(tribo_h2, ribo_tot.-solDF2[:, :ribo_d].-solDF2[:, :ribo_t].-solDF2[:, :rtrtcb].-solDF2[:, :rdrtca])
end
print(length(sols2))
p12 = plot(sols2[1][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 0");
p22 = plot(sols2[2][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 4");
p32 = plot(sols2[3][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 8");
p42 = plot(sols2[4][2:end], yaxis=(:log10, (1,Inf)),xaxis=(:log10, (1,Inf)), title="Rd = 12");
p52 = plot(sols2[5][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 16");
p62 = plot(sols2[6][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 20");
plot(p12, p22, p32, p42, p52, p62, size=(800,600), labels=labels_mix2, plot_title="Rh as algebraic")

t12 = plot(sols_t2, trtcr2, xaxis=(:log10, (1,Inf)), title="RtcR");
t22 = plot(sols_t2, trtca2, xaxis=(:log10, (1,Inf)), title="RtcA");
t32 = plot(sols_t2, trtcb2, xaxis=(:log10, (1,Inf)), title="RtcB");
t42 = plot(sols_t2, tribo_d2, xaxis=(:log10, (1,Inf)), title="Rd");
t52 = plot(sols_t2, tribo_t2, xaxis=(:log10, (1,Inf)), title="Rt");
t62 = plot(sols_t2, tribo_h2, xaxis=(:log10, (1,Inf)), title="Rh");
plot(t12, t22, t32, t42, t52, t62, size=(800,600), labels=["0" "4" "8" "12" "16" "20"], plot_title="Effect of Rd")


params_alg = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot, ribo_tot]
init_alg = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_t_0];

sols3, sols_t3, trm_a3, trtca3, trm_b3, trtcb3, trm_r3, trtcr3, tribo_d3, tribo_h3, tribo_t3, trdrtca3, trtrtcb3 = empty_res()

for i in rd_range
    init_alg[7] = i
    sol3 = simple_solve!(full_alg!, init_alg, tspan, params_alg)
    push!(sols3, sol3)
    push!(sols_t3, sol3.t)
    solDF3 = DataFrame([[j[i] for j in sol3.u] for i=1:length(sol3.u[1])], species_alg)
    push!(trm_a3, solDF3[:, :rm_a])
    push!(trtca3, solDF3[:, :rtca])
    push!(trm_b3, solDF3[:, :rm_b])
    push!(trtcb3, solDF3[:, :rtcb])
    push!(trm_r3, solDF3[:, :rm_r])
    push!(trtcr3, solDF3[:, :rtcr])
    push!(tribo_d3, solDF3[:, :ribo_d])
    push!(tribo_t3, solDF3[:, :ribo_t])
    push!(trtrtcb3, ka_b.*solDF3[:, :ribo_t].*rtcb_tot./(ka_b.*solDF3[:, :ribo_t].+kb_b.+kc_b.*atp))
    push!(trdrtca3, k1_a.*solDF3[:, :ribo_d].*rtca_tot./(k1_a.*solDF3[:, :ribo_d].+k2_a.+k3_a.*atp))
    push!(tribo_h3, ribo_tot.-solDF3[:, :ribo_d].-solDF3[:, :ribo_t].-(ka_b.*solDF3[:, :ribo_t].*rtcb_tot./(ka_b.*solDF3[:, :ribo_t].+kb_b.+kc_b.*atp)).-(k1_a.*solDF3[:, :ribo_d].*rtca_tot./(k1_a.*solDF3[:, :ribo_d].+k2_a.+k3_a.*atp)))
end

p13 = plot(sols3[1][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 0");
p23 = plot(sols3[2][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 4");
p33 = plot(sols3[3][2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), title="Rd = 8");
p43 = plot(sols3[4][2:end], #=yaxis=(:log10, (1,Inf)),=# xaxis=(:log10, (1,Inf)), title="Rd = 12");
p53 = plot(sols3[5][2:end], #=yaxis=(:log10, (1,Inf)),=# xaxis=(:log10, (1,Inf)), title="Rd = 16");
p63 = plot(sols3[6][2:end], #=yaxis=(:log10, (1,Inf)),=# xaxis=(:log10, (1,Inf)), title="Rd = 20");
plot(p13, p23, p33, p43, p53, p63, size=(800,600), labels=labels_alg, plot_title="All algebraic")

t13 = plot(sols_t3, trtcr3, xaxis=(:log10, (1,Inf)));
t23 = plot(sols_t3, trtca3, xaxis=(:log10, (1,Inf)));
t33 = plot(sols_t3, trtcb3, xaxis=(:log10, (1,Inf)));
t43 = plot(sols_t3, tribo_d3, xaxis=(:log10, (1,Inf)));
t53 = plot(sols_t3, tribo_t3, xaxis=(:log10, (1,Inf)));
t63 = plot(sols_t3, tribo_h3, xaxis=(:log10, (1,Inf)));
plot(t13, t23, t33, t43, t53, t63, size=(800,600), labels=["0" "4" "8" "12" "16" "20"], plot_title="Effect of Rd")




















rh_range = range(0, 20, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in full model to ode and see if that works (probs will) 

params_full = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k]
init_full = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_h_0, ribo_d_0, ribo_t_0];
plot1_rh_ribo, plot2_rh_ribo = plots_for_param("Rh", init_full, 9, rh_range, full_ODEs!, init_full, tspan, params_full, species_full)
# plot1_rd_ribo
plot2_rh_ribo # all ribosomes as algebraic

params_mix2 = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot]
init_mix2 = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_d_0, ribo_t_0];
plot1_rh_ribod, plot2_rh_ribod = plots_for_param_rh("Rh", init_mix2, 7, rh_range, rtc_system_ribos_h!, init_riboh, tspan, params_ribo, species_riboh)
# plot1_rd_riboh
plot2_rh_ribod # rd as algebraic


init = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0];
params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, #=gr,=# d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_t_0#=, ribo_tot=#];
plot1_rh, plot2_rh = plots_for_param("Rh", init, 8, rh_range, rtc_system!, init, tspan, params, species)
# plot1_rd
plot2_rh # rt as algebraic



rt_range = range(0, 20, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in full model to ode and see if that works (probs will) 

init_ribo = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, ribo_t_0];
params_ribo = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0]
plot1_rt_ribo, plot2_rt_ribo = plots_for_param("Rt", init_ribo, 9, rt_range, rtc_system_ribos!, init_ribo, tspan, params_ribo, species_ribo)
# plot1_rd_ribo
plot2_rt_ribo # all ribosomes as algebraic


init_ribod = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_h_0, ribo_t_0];
params_ribo = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0]
plot1_rt_ribod, plot2_rt_ribod = plots_for_param_rh("Rt", init_ribod, 8, rt_range, rtc_system_ribos_h!, init_riboh, tspan, params_ribo, species_riboh)
# plot1_rd_riboh
plot2_rt_ribod # rd as algebraic


init_riboh = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_t_0];
params_ribo = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0]
plot1_rt_riboh, plot2_rt_riboh = plots_for_param_rh("Rt", init_riboh, 8, rt_range, rtc_system_ribos_h!, init_riboh, tspan, params_ribo, species_riboh)
# plot1_rt_riboh
plot2_rt_riboh # rh as algebraic













function split_rh(species, solution) # 2

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]

    n = 6
    rdrtca = k1_a.*ribo_d.*rtca./(k1_a.*ribo_d.+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*ribo_t.*rtcb./(ka_b.*ribo_t.+kb_b.+kc_b.*atp)
   
    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0 #10
    ribo_h = ribo_tot .- ribo_d .- ribo_t .- rdrtca .- rtrtcb

    alpha = ribo_t./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n)+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*rtcr # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel



    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end