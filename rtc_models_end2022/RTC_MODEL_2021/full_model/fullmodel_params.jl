
k1 = 0.1; 
k2 = 0.1;
k3 = 0.1;
k4 = 0.1;
k5 = 10;
sigma = 10# 10;
atp = 10;
# rtcr = 10; #10;
# rt = 10; # 10;
kr = 10;
# global kt = 1000;
L = 100; # t0/r0
c = 0.01; #kr/kt
w_rtc = 4.14 # 5 # 10;    # max. transcription rate 
theta_rtc = 20 # 4.38 # 5000;  # transcription threshold
max = 4;            # max. translation rate 
thr = 20 #7 # 5000;        # translation threshold 
# r = 10;
d = 0.01;

#ribosomes
k1_a = 0.1; #0
k2_a = 0.1; #0
k3_a = 1; #0
ka_b = 0.1; #0
kb_b = 0.1; #0
kc_b = 1; #0
k = 0.05; #0
d1 = 0.001;
# rtc
rm_a_0 = 0;
rtca_0 = 1;
rm_b_0 = 0;
rtcb_0 = 1;
rm_r_0 = 0;
rtcr_0 = 0;

#ribosomes 
rdrtca_0 = 0;
rtrtcb_0 = 0;
ribo_d_0 = 0;
ribo_t_0 = 0;
ribo_h_0 = 10;

ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0
rtca_tot = rtca_0 + rdrtca_0
rtcb_tot = rtcb_0 + rtrtcb_0

# gr = 0.01*ribo_h

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_t_0#=ribo_tot=#];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_d_0, ribo_t_0];

init = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0];

tspan = (0,1e9)

labels_ode = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rd" "Rt"]

labels = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rd" "Rh"]
labels_ribo = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rh" "Rd" "Rt"]
labels_riboh = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rd" "Rt"]
labels_ribod = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rh" "Rt"]
labels_int = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rh" "Rd" "Rt"]


params_mix1b = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtcb_tot]
init_mix1b = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, ribo_h_0, ribo_d_0, ribo_t_0];
labels_mix1b = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "Rh" "Rd" "Rt"]


params_mix1 = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot]
init_mix1 = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_h_0, ribo_d_0, ribo_t_0];
labels_mix1 = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rh" "Rd" "Rt"]
species_mix1 = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :ribo_h, :ribo_d, :ribo_t]

params_mix2 = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot]
init_mix2 = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_d_0, ribo_t_0];
labels_mix2 = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rd" "Rt"]
species_mix2 = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rdrtca, :rtrtcb, :ribo_d, :ribo_t]

params_full = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k]
init_full = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_h_0, ribo_d_0, ribo_t_0];
labels_full = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rh" "Rd" "Rt"]
species_full = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rdrtca, :rtrtcb, :ribo_h, :ribo_d, :ribo_t]

params_alg = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot]
init_alg = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_t_0];
labels_alg = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rd" "Rt"]
species_alg = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :ribo_d, :ribo_t]

labels_all = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rd" "Rh" "Rt"]


params_rtc = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot]
init_rtc = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_h_0, ribo_d_0, ribo_t_0]
labels_rtc = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "Rh" "Rd" "Rt"]
species_rtc = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :ribo_h, :ribo_d, :ribo_t]

params_new = [kr, L, c, sigma, atp, w_rtc, theta_rtc, max, thr, d, k1, k2, k3, k4, k5, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, d1]
init_new = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rdrtca_0, rtrtcb_0, ribo_h_0, ribo_d_0, ribo_t_0]
labels_new = ["mRNA-RtcA" "RtcA" "mRNA-RtcB" "RtcB" "mRNA-RtcR" "RtcR" "RdRtcA" "RtRtcB" "Rh" "Rd" "Rt"]

