include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")

rh_val = 30#75 # conc of ribosomes in exponential phase 
thr_t_val = 4#25/90#100 #30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna_val = 0.00175
kdeg_trna_val = 0.00001

k_inhib2_val_trna = 0.0025
inhib_val_trna = 0.1

k_inhib1_val_trna = 0.5
k_inhib1a_val_trna = 0.15
k_inhib1b_val_trna = 0.25

k_trna_inhib_vals = [k_inhib1a_val_trna, k_inhib1b_val_trna, k_inhib1_val_trna]

# species_trna = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
# species_trna_inhib = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt, :rtc_i]




# trna_init = 135.5
# lam*trna_init/((g_max*atp/(θtlr+atp))*(trna_init/(thr_t+trna_init)))

# kin_trna*(g_max*atp/(θtlr+atp))