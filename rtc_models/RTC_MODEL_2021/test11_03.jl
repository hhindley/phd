using Plots, DifferentialEquations, DataFrames, Statistics, Measures, LSODA
include("full_model/fullmodel_functions.jl")
include("full_model/fullmodel_params.jl")
include("full_model/models_new.jl")

prob = ODEProblem(new_ODEs!, init_new, tspan, params_new)
sol = solve(prob, Rodas4(), abstol=1e-12, reltol=1e-9)
p = plot(sol[2:end], ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="full ODEs")

prob1 = ODEProblem(RTC_SYSTEM!, init_rtc, tspan, params_new)
sol1 = solve(prob1, Rodas4(), abstol=1e-15, reltol=1e-12)
p1 = plot(sol1[2:end], yaxis=(:log10, (1,Inf)),  ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)), labels=labels_rtc, title="reduced model");

# prob1 = ODEProblem(RTC_SYSTEM_rh_alg!, init_alg, tspan, params_alg)
# sol1 = solve(prob1, Rodas4(), abstol=1e-15, reltol=1e-12)
# p1 = plot(sol1[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_rtc, title="reduced model");

prob2 = ODEProblem(RTC_SYSTEM2!, init_rtc, tspan, params_new)
sol2 = solve(prob2, Rodas4(), abstol=1e-15, reltol=1e-12)
p2 = plot(sol2[2:end], yaxis=(:log10, (1,Inf)), ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)), labels=labels_rtc, title="reduced model");

plot(p, p1, p2, layout=(3,1), size=(800,600))

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_full_model(species_full, sol)
# rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_h1, ribo_d1, ribo_t1 = split_reduced_model_rh(species_alg, sol1)

rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_h1, ribo_d1, ribo_t1 = split_reduced_model(species_rtc, sol1)
rm_a2, rtca2, rm_b2, rtcb2, rm_r2, rtcr2, rdrtca2, rtrtcb2, ribo_h2, ribo_d2, ribo_t2 = split_reduced_model(species_rtc, sol2)

# p_full = plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="full ODEs");
# p1_full = plot(sol1.t, [rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_d1, ribo_h1, ribo_t1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="reduced model");
# plot(p_full, p1_full, size=(800,500))

plot(sol.t, [rm_a, rm_b, rm_r], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["mRNA-A-full" "mRNA-B-full" "mRNA-R-full"])
plot!(sol1.t, [rm_a1, rm_b1, rm_r1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["mRNA-A" "mRNA-B" "mRNA-R"])
plot!(sol2.t, [rm_a2, rm_b2, rm_r2], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["mRNA-A" "mRNA-B" "mRNA-R"])

plot(sol.t, [rtca, rtcb, rtcr], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA-full" "RtcB-full" "RtcR-full"])
plot!(sol1.t, [rtca1, rtcb1, rtcr1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"])
plot!(sol2.t, [rtca2, rtcb2, rtcr2], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"])

plot(sol.t, [rdrtca], xaxis=(:log10, (1,Inf)), labels="RdRtcA-full")
plot!(sol1.t, [rdrtca1], xaxis=(:log10, (1,Inf)), labels="RdRtcA-reduced")
plot!(sol2.t, [rdrtca2], xaxis=(:log10, (1,Inf)), labels="RdRtcA")

plot(sol.t, [rtrtcb], xaxis=(:log10, (1,Inf)), labels="RtRtcB-full")
plot!(sol1.t, [rtrtcb1], xaxis=(:log10, (1,Inf)), labels="RtRtcB-reduced")
plot!(sol2.t, [rtrtcb2], xaxis=(:log10, (1,Inf)), labels="RtRtcB")
print(maximum(rtrtcb2[1]))
plot(sol.t, ribo_d, xaxis=(:log10, (1,Inf)), labels="Rd-full")
plot!(sol1.t, ribo_d1, xaxis=(:log10, (1,Inf)), labels="Rd-reduced")
plot!(sol2.t, ribo_d2, xaxis=(:log10, (1,Inf)), labels="Rd")

plot(sol.t, ribo_h, xaxis=(:log10, (1,Inf)), labels="Rh-full")
plot!(sol1.t, ribo_h1, xaxis=(:log10, (1,Inf)), labels="Rh-reduced")
plot!(sol2.t, ribo_h2, xaxis=(:log10, (1,Inf)), labels="Rh")

plot(sol.t, ribo_t, xaxis=(:log10, (1,Inf)), labels="Rt-full")
plot!(sol1.t, ribo_t1, xaxis=(:log10, (1,Inf)), labels="Rt-reduced")
plot!(sol2.t, ribo_t2, xaxis=(:log10, (1,Inf)), labels="Rt")

plot(sol.t, [ribo_d, ribo_h, ribo_t], xaxis=(:log10, (1,Inf)), labels=["Rd-full" "Rh-full" "Rt-full"])
plot!(sol1.t, [ribo_d1, ribo_h1, ribo_t1], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"])
plot!(sol2.t, [ribo_d2, ribo_h2, ribo_t2], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"])

plot(sol.t, ribo_d+ribo_h+ribo_t+rdrtca+rtrtcb, xaxis=(:log10, (1,Inf)), labels="Rtot-full", title="Full ODEs")
plot(sol1.t, (ribo_d1)+ribo_h1+(ribo_t1), xaxis=(:log10, (1,Inf)), labels="Rt+Rd+Rh+RdRtcA+RtRtcB", title="reduced model")
plot(sol2.t, ribo_d2+ribo_h2+ribo_t2, xaxis=(:log10, (1,Inf)), labels="Rtot")

print((ribo_d2[10000]))
print((ribo_h2[10000]))
print((ribo_t2[10000]))
print((rdrtca2[10000]))
print((rtrtcb2[10000]))
print(ribo_d2[10000]+ribo_h2[10000]+ribo_t2[10000]+rdrtca2[10000]+rtrtcb2[10000])

rtot = ribo_d2+ribo_h2+ribo_t2+rdrtca2+rtrtcb2

print(rtot[10000])
rtot_10 = []
for i in rtot
    if i > 10.2
        push!(rtot_10, i)
    end
end

println(length(rtot_10))