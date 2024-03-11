using BlackBoxOptim, Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")

# set time span and how many time points to solve at 
tspan = (0, 100)
t = exp10.(range(-3,2,15))
pushfirst!(t, 0)

sol_syn = sol_with_t(rtc_model, init, tspan, params, t)

# plot solution 
Plots.plot(sol_syn[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))

rtca_syn = get_curve(sol_syn, :rtca)
rtcr_syn = get_curve(sol_syn, :rtcr)
rd_syn = get_curve(sol_syn, :rd)
rh_syn = get_curve(sol_syn, :rh)

function rtc_bo1(x)
    solu = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, x[3], ktag, kdeg, kin, atp]), t)
    rtca = get_curve(solu, :rtca)
    rtcr = get_curve(solu, :rtcr)
    obj = []
    for (i,j) in zip(rtca, rtca_syn)
        append!(obj, abs2(i-j))
    end
    for (i,j) in zip(rtcr, rtcr_syn)
        append!(obj, abs2(i-j))
    end
    return sum(obj)

        # kdam 
    for (i,j) in zip(rh, rh_syn)
        append!(obj, abs2(i-j))
    end
    for (i,j) in zip(rd, rd_syn)
        append!(obj, abs2(i-j))
    end
end

res = bboptimize(rtc_bo1; SearchRange = [(0, 10.0), (0, 10.0), (0, 1.0)], NumDimensions = 3, MaxSteps = 3500)

best_candidate(res)
best_fitness(res)