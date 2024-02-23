include("/home/holliehindley/phd/rtc_model/parameters/params.jl")

rh_val = 30#75 # conc of ribosomes in exponential phase 
thr_t_val = 4#25/90#100 #30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna_val = 0.00175
kdeg_trna_val = 0.00001

k_inhib2_val_trna = 0.0025
inhib_val_trna = 0.1

k_inhib1_val_trna = 0.1
k_inhib1a_val_trna = 0.03
k_inhib1b_val_trna = 0.05

k_trna_inhib_vals = [k_inhib1a_val_trna, k_inhib1b_val_trna, k_inhib1_val_trna]

species_trna = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
species_trna_inhib = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt, :rtc_i]

init_trna = [rtc_trna_model.rm_a=>0.0,rtc_trna_model.rtca=>0.0,rtc_trna_model.rm_b=>0.0,rtc_trna_model.rtcb=>0.0,rtc_trna_model.rm_r=>0.0,rtc_trna_model.rtcr=>0.0,rtc_trna_model.trna=>135.5,rtc_trna_model.rd=>0.0,rtc_trna_model.rt=>0.0] # tRNA initial conc = 135.5

tspan = (0, 1e9);

params_trna = Dict(L=>L_val, c=>c_val, kr=>kr_val*12, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_trna_val, kin=>kin_trna_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, rh=>rh_val, thr_t=>thr_t_val)

solu = sol(rtc_trna_model, init_trna, tspan, params_trna)
ssvals_trna = get_all_ssvals(solu, species_trna)


# trna_init = 135.5
# lam*trna_init/((g_max*atp/(θtlr+atp))*(trna_init/(thr_t+trna_init)))

# kin_trna*(g_max*atp/(θtlr+atp))