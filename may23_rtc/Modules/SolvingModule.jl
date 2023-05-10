# using StaticArrays, LabelledArrays, DifferentialEquations

module SolvingModule
export sol
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl")
end


# using .Solving, .Models
# include("/home/holliehindley/phd/may23_rtc/parameters/init.jl")
# include("/home/holliehindley/phd/may23_rtc/parameters/params.jl")
# initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];
# params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

# # prob = ODEProblem(rtc_model, initial, tspan, params)
# # solu = solve(prob, Rodas4())


# sol(rtc_model, initial, tspan, params)