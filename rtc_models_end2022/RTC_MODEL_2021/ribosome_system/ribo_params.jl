k1 = 0.1
k2 = 0.1
k3 = 1
ka = 0.1
kb = 0.1
kc = 1
k = 0.1
atp = 10

rtca_0 = 1
rtcb_0 = 1
rdrtca_0 = 0
rtrtcb_0 = 0
ribo_d_0 = 0
ribo_t_0 = 0
ribo_h_0 = 10

params = [k1, k2, k3, ka, kb, kc, k, atp]
init = [rtca_0, rtcb_0, rdrtca_0, rtrtcb_0, ribo_d_0, ribo_t_0, ribo_h_0]
tspan = (0, 1000)

species = [:rtca, :rtcb, :ribo_d, :ribo_h]
species_full = [:rtca, :rtcb, :rdrtca, :rtrtcb, :ribo_d, :ribo_t, :ribo_h]


rtot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0
params_alg = [k1, k2, k3, ka, kb, kc, k, atp, rdrtca_0, rtrtcb_0, rtot]
init_alg = [rtca_0, rtcb_0, ribo_d_0, ribo_h_0]

rtca_tot = rtca_0 + rdrtca_0
rtcb_tot = rtcb_0 + rtrtcb_0
ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0 

new_params = [k1, k2, k3, ka, kb, kc, k, atp, rtca_tot, rtcb_tot, ribo_tot]
new_species = [:rtca, :rtcb, :ribo_d, :ribo_t]
new_init = [rtca_0, rtcb_0, ribo_d_0, ribo_t_0]

params_mix1 = [k1, k2, k3, ka, kb, kc, k, atp, rtca_tot, rtcb_tot]
params_mix2 = [k1, k2, k3, ka, kb, kc, k, atp, ribo_tot]
init_mix1 = [rtca_0, rtcb_0, ribo_d_0, ribo_t_0, ribo_h_0]
init_mix2 = [rtca_0, rtcb_0, ribo_d_0, ribo_t_0, rdrtca_0, rtrtcb_0]
species_mix1 = [:rtca, :rtcb, :ribo_d, :ribo_t, :ribo_h]
species_mix2 = [:rtca, :rtcb, :ribo_d, :ribo_t, :rdrtca, :rtrtcb]