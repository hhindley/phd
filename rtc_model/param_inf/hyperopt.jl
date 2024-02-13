using Hyperopt, Plots
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

# set time span and how many time points to solve at 
tspan = (0, 100)
t = exp10.(range(-3,2,15))
pushfirst!(t, 0)

sol_syn = sol(rtc_model, init, tspan, params, t)

# plot solution 
Plots.plot(sol_syn[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))

rtca_syn = get_curve(sol_syn, :rtca)
rtcr_syn = get_curve(sol_syn, :rtcr)

function rtc_bo(ω_ab)#, ω_r)
    solu = sol(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]), t)
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
end



ho = @hyperopt for i=1000,
    sampler = BOHB(), # This is default if none provided
    ω_ab = LinRange(1,7, 1000)
    # ω_r = LinRange(1,7, 1000)
print(i, "\t", ω_ab, "\t")#, ω_r, "\t")
@show  rtc_bo(ω_ab)
end

Plots.plot(ho, size=(1000,800))

print(ho.history)
print(ho.results)
print(ho.candidates)
print(ho.params)
print(ho.minimum) # minumum of function rtc_bo

print(ho.minimizer) # best params

printmin(ho)