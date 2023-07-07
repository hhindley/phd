using Plots, Printf, Measures
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")


params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014) 	

params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014) 	

br = get_br(rtc_mod, params1, initial, 3.)

plot(br)




plots = []
for i in (10 .^ range(-4, stop=0, length = 5))
    for j in (10 .^ range(-2, stop=1, length = 5))
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = j, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:ω_r])
        push!(plots, plot(get_br(rtc_mod, params1, initial, 3.), title="ω_ab = $(@sprintf "%.2E" j), ω_r = $(@sprintf "%.2E" i)"))
    end
end

p_all = plot(plots[1], plots[2], plots[3], plots[4], plots[5], plots[6],
plots[7], plots[8], plots[9], plots[10], plots[11], plots[12],
plots[13], plots[14], plots[15], plots[16], plots[17], plots[18],
plots[19], plots[20], plots[21], plots[22], plots[23], plots[24], plots[25],
 margin=1mm, size=(2000,2000), layout=(5,5))

plots_br = []
for i in range(1,stop=5000,length=25)
    params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = i, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014) 	
    
    push!(plots_br, plot(get_br(rtc_mod, params1, initial, 3.), title="atp=$i"))
end

plot(plots_br[1], plots_br[2], plots_br[3], plots_br[4], plots_br[5], plots_br[6],
plots_br[7], plots_br[8], plots_br[9], plots_br[10], plots_br[11], plots_br[12],
plots_br[13], plots_br[14], plots_br[15], plots_br[16], plots_br[17], plots_br[18],
plots_br[19], plots_br[20], plots_br[21], plots_br[22], plots_br[23], plots_br[24], plots_br[25],
 margin=1mm, size=(2000,2000), layout=(5,5))