using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, Plots #, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")


prob = ODEProblem(rtc_model, init, tspan, params)
solu = solve(prob, Rodas4())

Plots.plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))


function get_curve(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    species = df[:, species]
    return species
end

rm_a = get_curve(solu, :rm_a)
rtca = get_curve(solu, :rtca)
rm_b = get_curve(solu, :rm_b)
rtcb = get_curve(solu, :rtcb)
rm_r = get_curve(solu, :rm_r)
rtcr = get_curve(solu, :rtcr)
rh = get_curve(solu, :rh)
rd = get_curve(solu, :rd)
rt = get_curve(solu, :rt)

Plots.plot(solu.t, rtca, xaxis=(:log10, (0.01,Inf)), labels="RtcA")
Plots.plot(solu.t, rtcb, xaxis=(:log10, (0.01,Inf)), labels="RtcB")
Plots.plot(solu.t, rtcr, xaxis=(:log10, (0.01,Inf)), labels="RtcR")
Plots.plot(solu.t, rm_a, xaxis=(:log10, (0.01,Inf)), labels="rm_a")
Plots.plot(solu.t, rm_b, xaxis=(:log10, (0.01,Inf)), labels="rm_b")
Plots.plot(solu.t, rm_r, xaxis=(:log10, (0.01,Inf)), labels="rm_r")

Plots.plot(solu.t, rt, xaxis=(:log10, (0.01,Inf)), labels="Rt")
Plots.plot(solu.t, rh, xaxis=(:log10, (0.01,Inf)), labels="Rh")
Plots.plot(solu.t, rd, xaxis=(:log10, (0.01,Inf)), labels="Rd")

r_tot = rh+rt+rd
Plots.plot(solu.t, (@.rh/r_tot *100), xaxis=(:log10, (0.01,Inf)), labels="Rh")
Plots.plot!(solu.t, (@.rt/r_tot *100), xaxis=(:log10, (0.01,Inf)), labels="Rt")
Plots.plot!(solu.t, (@.rd/r_tot *100), xaxis=(:log10, (0.01,Inf)), labels="Rd")

v = collect(1:1:100)

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rt, rd, rtca1 = [], [], [], [], [], [], [], [], [], []
for i in v
    k_a = i
    params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, k_a, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
    prob = ODEProblem(rtc_model_new, init, tspan, params)
    solu = solve(prob, Rodas4())
    solDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    push!(rm_a, solDF[end, :rm_a])
    push!(rtca, solDF[end, :rtca])
    push!(rm_b, solDF[end, :rm_b])
    push!(rtcb, solDF[end, :rtcb])
    push!(rm_r, solDF[end, :rm_r])
    push!(rtcr, solDF[end, :rtcr])
    push!(rh, solDF[end, :rh])
    push!(rt, solDF[end, :rt])
    push!(rd, solDF[end, :rd])
    push!(rtca1, (atp*solDF[end, :rtca])/(atp+(ktag*solDF[end, :rd])/k_a))
end

plot(v, rt, markershapes=[:circle])

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rt, rd, rtca1, rtcb1 = [], [], [], [], [], [], [], [], [], [], []
for i in v1
    KM = i
    params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, KM, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
    # println(params[12])
    prob = ODEProblem(rtc_model_new, init, tspan, params)
    solu = solve(prob, Rodas4())
    solDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    push!(rm_a, solDF[end, :rm_a])
    push!(rtca, solDF[end, :rtca])
    push!(rm_b, solDF[end, :rm_b])
    push!(rtcb, solDF[end, :rtcb])
    push!(rm_r, solDF[end, :rm_r])
    push!(rtcr, solDF[end, :rtcr])
    push!(rh, solDF[end, :rh])
    push!(rt, solDF[end, :rt])
    push!(rd, solDF[end, :rd])
    push!(rtca1, (atp*solDF[end, :rtca])/(atp+(ktag*solDF[end, :rd])/k_a))
    push!(rtcb1, (atp*solDF[end, :rtcb])/(atp+(krep*solDF[end, :rt])/k_b))
end
plot(v1, rd, markershapes=[:circle])

v1 = collect(0:100:1000)

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rt, rd, rtca1, rtcb1 = [], [], [], [], [], [], [], [], [], [], []
for i in v1
    ktag = i
    params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, KM, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
    prob = ODEProblem(rtc_model_new, init, tspan, params)
    solu = solve(prob, Rodas4())
    solDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    push!(rm_a, solDF[end, :rm_a])
    push!(rtca, solDF[end, :rtca])
    push!(rm_b, solDF[end, :rm_b])
    push!(rtcb, solDF[end, :rtcb])
    push!(rm_r, solDF[end, :rm_r])
    push!(rtcr, solDF[end, :rtcr])
    push!(rh, solDF[end, :rh])
    push!(rt, solDF[end, :rt])
    push!(rd, solDF[end, :rd])
    push!(rtca1, (atp*solDF[end, :rtca])/(atp+(ktag*solDF[end, :rd])/k_a))
    push!(rtcb1, (atp*solDF[end, :rtcb])/(atp+(krep*solDF[end, :rt])/k_b))
end
plot(v1, rt, markershapes=[:circle])


rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rt, rd, rtca1, rtcb1 = [], [], [], [], [], [], [], [], [], [], []
for i in v1
    krep = i
    params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, k_a, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
    prob = ODEProblem(rtc_model_new, init, tspan, params)
    solu = solve(prob, Rodas4())
    solDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    push!(rm_a, solDF[end, :rm_a])
    push!(rtca, solDF[end, :rtca])
    push!(rm_b, solDF[end, :rm_b])
    push!(rtcb, solDF[end, :rtcb])
    push!(rm_r, solDF[end, :rm_r])
    push!(rtcr, solDF[end, :rtcr])
    push!(rh, solDF[end, :rh])
    push!(rt, solDF[end, :rt])
    push!(rd, solDF[end, :rd])
    push!(rtca1, (atp*solDF[end, :rtca])/(atp+(ktag*solDF[end, :rd])/k_a))
    push!(rtcb1, (atp*solDF[end, :rtcb])/(atp+(krep*solDF[end, :rt])/k_b))
end
plot(v1, rh, markershapes=[:circle])

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rt, rd, rtca1, rtcb1 = [], [], [], [], [], [], [], [], [], [], []
for i in collect(0.1:10:1000)
    atp = i
    params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, k_a, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
    prob = ODEProblem(rtc_model_new, init, tspan, params)
    solu = solve(prob, Rodas4())
    solDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    push!(rm_a, solDF[end, :rm_a])
    push!(rtca, solDF[end, :rtca])
    push!(rm_b, solDF[end, :rm_b])
    push!(rtcb, solDF[end, :rtcb])
    push!(rm_r, solDF[end, :rm_r])
    push!(rtcr, solDF[end, :rtcr])
    push!(rh, solDF[end, :rh])
    push!(rt, solDF[end, :rt])
    push!(rd, solDF[end, :rd])
    push!(rtca1, (atp*solDF[end, :rtca])/(atp+(ktag*solDF[end, :rd])/k_a))
    push!(rtcb1, (atp*solDF[end, :rtcb])/(atp+(krep*solDF[end, :rt])/k_b))
end
plot(collect(0.1:10:1000), rtca)#, markershapes=[:circle])