include("/home/holliehindley/phd/rtc_model/parameters/params.jl")

rh = 11.29 #75 # conc of ribosomes in exponential phase 
thr_t = 10 #30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 0.00175
kdeg_trna = 0.00001

trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]

init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5

params_trna = @LArray [L, c, kr*12, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg_trna, kin_trna, atp, nA, nB, nR, lam, kc, k_diss, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :kc, :k_diss, :rh, :thr_t)

params_trna_bf = (L = L, c = c, kr = kr*12, Vmax_init = Vmax_init, Km_init = Km_init,
ω_ab = ω_ab, ω_r = ω_r, θtscr = θtscr, g_max = g_max, θtlr = θtlr, km_a = km_a, km_b = km_b, 
d = d, krep = krep, kdam = kdam, ktag = ktag, kdeg = kdeg_trna, kin = kin_trna, atp = atp, 
na = nA, nb = nB, nr = nR, lam = lam, kc = kc, k_diss = k_diss, rh = rh, thr_t = thr_t)

trna_init = 135.5
lam*trna_init/((g_max*atp/(θtlr+atp))*(trna_init/(thr_t+trna_init)))