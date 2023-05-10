using Plots, Measures, Interpolations

include("MWC_functions.jl");

#ranges and labels
rtcr_range = range(0, 100, length=100);
rt_range = range(0, 100, length=100);

c_range = [0, 0.01, 0.1, 0.15, 0.2];
c_labels = ["c=0" "c=0.01" "c=0.1" "c=0.15" "c=0.2"];

kr_range = [1, 5, 10, 25, 50];
kr_labels = ["kr=1" "kr=5" "kr=10" "kr=25" "kr=50"];

rtcr_range1 = [1, 5, 10, 15, 20];
rtcr_labels = ["RtcR=1" "RtcR=5" "RtcR=10" "RtcR=15" "RtcR=20"];

rt_range1 = [1, 5, 10, 15, 20];
rt_labels = ["Rt=1" "Rt=5" "Rt=10" "Rt=15" "Rt=20"];

L_range = exp10.(range(log10(0.1), log10(10000), length=5));
L_labels = ["L=0.1" "L=1.78" "L=31.62" "L=562.34" "L=10000"];

sigma_range = [0, 1, 5, 10, 20];
sigma_labels = ["\\sigma=0" "\\sigma=1" "\\sigma=5" "\\sigma=10" "\\sigma=20"];

atp_range = [0, 1, 5, 10, 50];
atp_labels = ["ATP=0" "ATP=1" "ATP=5" "ATP=10" "ATP=50"];

k_range = [0, 0.01, 0.1, 1, 10];
k_labels = ["k=0" "k=0.01" "k=0.1" "k=1" "k=10"];

#parameters
k1 = 0.1; 
k2 = 0.1;
k3 = 0.1;
k4 = 0.1;
k5 = 10;
sigma = 10;
atp = 10;
rtcr = 10; #10;
rt = 10; # 10;
kr = 10;
kt = 1000;
L = 100; # t0/r0
c = 0.01; #kr/kt
n = 6;

#plots
#changing rt or rtcR
v_rtcr = [vary_all(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp) for i in rtcr_range];
p_rtcr = plotting1(rtcr_range, v_rtcr, "Varying RtcR", "RtcR", "\\sigma*");

v_rt = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp) for i in rt_range];
p_rt = plotting1(rt_range, v_rt, "Varying Rt", "Rt", "\\sigma*");

# comparision to mm method 
p_rtcr_mm = plot_rtc_comp1(mm_rtcr, rtcr_range, "Varying RtcR with MM method", "RtcR", "\\sigma*");
p_rt_mm = plot_rtc_comp1(mm_rt, rt_range, "Varying Rt with MM method", "Rt", "\\sigma*");

#changing rt/rtcR and another parameter 
v_rtcr_c = [vary_all(rt, kr, n, L, j, i, k1, k2, k3, k4, k5, sigma, atp) for j in c_range for i in rtcr_range];
p_rtcr_c = plotting2(rtcr_range, v_rtcr_c, "Varying RtcR with c", "RtcR", "\\sigma*", c_labels);

v_rt_c = [vary_all(i, kr, n, L, j, rtcr, k1, k2, k3, k4, k5, sigma, atp) for j in c_range for i in rt_range];
p_rt_c = plotting2(rt_range, v_rt_c, "Varying Rt with c", "Rt", "\\sigma*", c_labels);

v_rtcr_kr = [vary_all(rt, j, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp) for j in kr_range for i in rtcr_range];
p_rtcr_kr = plotting2(rtcr_range, v_rtcr_kr, "Varying RtcR with kr", "RtcR", "\\sigma*", kr_labels);

v_rt_kr = [vary_all(i, j, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp) for j in kr_range for i in rt_range];
p_rt_kr = plotting2(rt_range, v_rt_kr, "Varying Rt with kr", "Rt", "\\sigma*", kr_labels);

v_rtcr_L = [vary_all(rt, kr, n, j, c, i, k1, k2, k3, k4, k5, sigma, atp) for j in L_range for i in rtcr_range];
p_rtcr_L = plotting2(rtcr_range, v_rtcr_L, "Varying RtcR with L", "RtcR", "\\sigma*", L_labels);

v_rt_L = [vary_all(i, kr, n, j, c, rtcr, k1, k2, k3, k4, k5, sigma, atp) for j in L_range for i in rt_range];
p_rt_L = plotting2(rt_range, v_rt_L, "Varying Rt with L", "Rt", "\\sigma*", L_labels);

v_rtcr_sigma = [vary_all(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, j, atp) for j in sigma_range for i in rtcr_range];
p_rtcr_sigma = plotting2(rtcr_range, v_rtcr_sigma, "Varying RtcR with \\sigma", "RtcR", "\\sigma*", sigma_labels);

v_rt_sigma = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, j, atp) for j in sigma_range for i in rt_range];
p_rt_sigma = plotting2(rt_range, v_rt_sigma, "Varying Rt with \\sigma", "Rt", "\\sigma*", sigma_labels);

v_rtcr_atp = [vary_all(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, j) for j in atp_range for i in rtcr_range];
p_rtcr_atp = plotting2(rtcr_range, v_rtcr_atp, "Varying RtcR with ATP", "RtcR", "\\sigma*", atp_labels);

v_rt_atp = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, j) for j in atp_range for i in rt_range];
p_rt_atp = plotting2(rt_range, v_rt_atp, "Varying Rt with ATP", "Rt", "\\sigma*", atp_labels);


#rt vs rtcR
v_rtcr_rt = [vary_all(j, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp) for j in rt_range1 for i in rtcr_range];
p_rtcr_rt = plotting2(rtcr_range, v_rtcr_rt, "Varying RtcR with Rt", "RtcR", "\\sigma*", rt_labels);

v_rt_rtcr = [vary_all(i, kr, n, L, c, j, k1, k2, k3, k4, k5, sigma, atp) for j in rtcr_range1 for i in rt_range];
p_rt_rtcr = plotting2(rt_range, v_rt_rtcr, "Varying Rt with RtcR", "Rt", "\\sigma*", rtcr_labels);


#changing rt/rtcR and rates 
v_rtcr_k1 = [vary_all(rt, kr, n, L, c, i, j, k2, k3, k4, k5, sigma, atp) for j in k_range for i in rtcr_range];
p_rtcr_k1 = plotting2(rtcr_range, v_rtcr_k1, "Varying RtcR with k1", "RtcR", "\\sigma*", k_labels);

v_rt_k1 = [vary_all(i, kr, n, L, c, rtcr, j, k2, k3, k4, k5, sigma, atp) for j in k_range for i in rt_range];
p_rt_k1 = plotting2(rt_range, v_rt_k1, "Varying Rt with k1", "Rt", "\\sigma*", k_labels);

v_rtcr_k2 = [vary_all(rt, kr, n, L, c, i, k1, j, k3, k4, k5, sigma, atp) for j in k_range for i in rtcr_range];
p_rtcr_k2 = plotting2(rtcr_range, v_rtcr_k2, "Varying RtcR with k2", "RtcR", "\\sigma*", k_labels);

v_rt_k2 = [vary_all(i, kr, n, L, c, rtcr, k1, j, k3, k4, k5, sigma, atp) for j in k_range for i in rt_range];
p_rt_k2 = plotting2(rt_range, v_rt_k2, "Varying Rt with k2", "Rt", "\\sigma*", k_labels);

v_rtcr_k3 = [vary_all(rt, kr, n, L, c, i, k1, k2, j, k4, k5, sigma, atp) for j in k_range for i in rtcr_range];
p_rtcr_k3 = plotting2(rtcr_range, v_rtcr_k3, "Varying RtcR with k3", "RtcR", "\\sigma*", k_labels);

v_rt_k3 = [vary_all(i, kr, n, L, c, rtcr, k1, k2, j, k4, k5, sigma, atp) for j in k_range for i in rt_range];
p_rt_k3 = plotting2(rt_range, v_rt_k3, "Varying Rt with k3", "Rt", "\\sigma*", k_labels);

v_rtcr_k4 = [vary_all(rt, kr, n, L, c, i, k1, k2, k3, j, k5, sigma, atp) for j in k_range for i in rtcr_range];
p_rtcr_k4 = plotting2(rtcr_range, v_rtcr_k4, "Varying RtcR with k4", "RtcR", "\\sigma*", k_labels);

v_rt_k4 = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, j, k5, sigma, atp) for j in k_range for i in rt_range];
p_rt_k4 = plotting2(rt_range, v_rt_k4, "Varying Rt with k4", "Rt", "\\sigma*", k_labels);

v_rtcr_k5 = [vary_all(rt, kr, n, L, c, i, k1, k2, k3, k4, j, sigma, atp) for j in k_range for i in rtcr_range];
p_rtcr_k5 = plotting2(rtcr_range, v_rtcr_k5, "Varying RtcR with k5", "RtcR", "\\sigma*", k_labels);

v_rt_k5 = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, k4, j, sigma, atp) for j in k_range for i in rt_range];
p_rt_k5 = plotting2(rt_range, v_rt_k5, "Varying Rt with k5", "Rt", "\\sigma*", k_labels);


#comparison plots 
plot(p_rtcr, p_rt, layout=(2,1), size=(600,600))

plot(p_rtcr, p_rt, p_rtcr_mm, p_rt_mm, layout=(2,2), size=(900, 600), margin=2mm)

plot(p_rtcr_c, p_rtcr_kr, p_rtcr_L, p_rtcr_sigma, p_rtcr_atp, p_rtcr_rt, layout=(3,2), size=(1000, 800), margin=2mm)

plot(p_rt_c, p_rt_kr, p_rt_L, p_rt_sigma, p_rt_atp, p_rt_rtcr, layout=(3,2), size=(1000, 800), margin=2mm)

plot(p_rtcr, p_rtcr_k1, p_rtcr_k2, p_rtcr_k3, p_rtcr_k4, p_rtcr_k5, layout=(3,2), size=(1000,800), margin=2mm)

plot(p_rt, p_rt_k1, p_rt_k2, p_rt_k3, p_rt_k4, p_rt_k5, layout=(3,2), size=(1000,800), margin=2mm)

# addition of transcription     
w_rtc = 4.14;
theta_rtc = 20;
w_rtc_range = [0, 1, 5, 10, 20];
w_rtc_labels = ["0" "1" "5" "10" "20"];
theta_rtc_range = [0, 1, 5, 10, 20];
theta_rtc_labels = ["0" "1" "5" "10" "20"];

c_range1 = range(0, 0.5, length=100);
L_range1 = range(0, 10000, length=100);
k_range1 = range(0, 100, length=100);
atp_range1 = range(0, 20000, length=100);
sigma_range1 = range(0, 50, length=100);


v_rtcr_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in rtcr_range];
p_rtcr_tr = plotting1(rtcr_range, v_rtcr_tr, "Varying RtcR on transcription", "RtcR", "\$RtcAB\$");

v_rt_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in rt_range];
p_rt_tr = plotting1(rt_range, v_rt_tr, "Varying Rt", "Rt", "\$RtcAB\$");

v_kr_tr = [vary_all_transcription(rt, i, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in rtcr_range];
p_kr_tr = plotting1(rtcr_range, v_kr_tr, "Varying kr on transcription", "kr", "\$RtcAB\$");

v_L_tr = [vary_all_transcription(rt, kr, n, i, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in L_range1];
p_L_tr = plotting1(L_range1, v_L_tr, "Varying L on transcription", "L", "\$RtcAB\$");

v_c_tr = [vary_all_transcription(rt, kr, n, L, i, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in c_range1];
p_c_tr = plotting1(c_range1, v_c_tr, "Varying c on transcription", "c", "\$RtcAB\$");

v_k1_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, i, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in k_range1];
p_k1_tr = plotting1(k_range1, v_k1_tr, "Varying k1 on transcription", "k1", "\$RtcAB\$");

v_k2_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, i, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in k_range1];
p_k2_tr = plotting1(k_range1, v_k2_tr, "Varying k2 on transcription", "k2", "\$RtcAB\$");

v_k3_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, i, k4, k5, sigma, atp, w_rtc, theta_rtc) for i in k_range1];
p_k3_tr = plotting1(k_range1, v_k3_tr, "Varying k3 on transcription", "k3", "\$RtcAB\$");

v_k4_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, i, k5, sigma, atp, w_rtc, theta_rtc) for i in k_range1];
p_k4_tr = plotting1(k_range1, v_k4_tr, "Varying k4 on transcription", "k4", "\$RtcAB\$");

v_k5_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, i, sigma, atp, w_rtc, theta_rtc) for i in k_range1];
p_k5_tr = plotting1(k_range1, v_k5_tr, "Varying k5 on transcription", "k5", "\$RtcAB\$");

v_atp_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, i, w_rtc, theta_rtc) for i in atp_range1];
p_atp_tr = plotting1(atp_range1, v_atp_tr, "Varying ATP on transcription", "ATP", "\$RtcAB\$");

v_sigma_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, i, atp, w_rtc, theta_rtc) for i in sigma_range1];
p_sigma_tr = plotting1(sigma_range1, v_sigma_tr, "Varying \\sigma on transcription", "\\sigma", "\$RtcAB\$");

v_wrtc_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, i, theta_rtc) for i in sigma_range1];
p_wrtc_tr = plotting1(sigma_range1, v_wrtc_tr, "Varying w_rtc on transcription", "w_rtc", "\$RtcAB\$");

v_trtc_tr = [vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, i) for i in sigma_range1];
p_trtc_tr = plotting1(sigma_range1, v_trtc_tr, "Varying theta_rtc on transcription", "theta_rtc", "\$RtcAB\$");




v_rtcr_wrtc = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp, j, theta_rtc) for j in w_rtc_range for i in rtcr_range];
p_rtcr_wrtc = plotting2(rtcr_range, v_rtcr_wrtc, "Vary RtcR and W_rtc for transcription", "RtcR", "\$RtcAB\$", w_rtc_labels);

v_rt_wrtc = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, j, theta_rtc) for j in w_rtc_range for i in rt_range];
p_rt_wrtc = plotting2(rt_range, v_rt_wrtc, "Vary Rt and W_rtc for transcription", "Rt", "\$RtcAB\$", w_rtc_labels);

v_rtcr_trtc = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, j) for j in theta_rtc_range for i in rtcr_range];
p_rtcr_trtc = plotting2(rtcr_range, v_rtcr_trtc, "Vary RtcR and theta_rtc for transcription", "RtcR", "\$RtcAB\$", theta_rtc_labels);

v_rt_trtc = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, j) for j in theta_rtc_range for i in rt_range];
p_rt_trtc = plotting2(rt_range, v_rt_trtc, "Vary Rt and theta_rtc for transcription", "Rt", "\$RtcAB\$", theta_rtc_labels);



v_rtcr_c_tr = [vary_all_transcription(rt, kr, n, L, j, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in c_range for i in rtcr_range];
p_rtcr_c_tr = plotting2(rtcr_range, v_rtcr_c, "Varying RtcR with c", "RtcR", "\$RtcAB\$", c_labels);

v_rt_c_tr = [vary_all_transcription(i, kr, n, L, j, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in c_range for i in rt_range];
p_rt_c_tr = plotting2(rt_range, v_rt_c, "Varying Rt with c", "Rt", "\$RtcAB\$", c_labels);

v_rtcr_kr_tr = [vary_all_transcription(rt, j, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in kr_range for i in rtcr_range];
p_rtcr_kr_tr = plotting2(rtcr_range, v_rtcr_kr, "Varying RtcR with kr", "RtcR", "\$RtcAB\$", kr_labels);

v_rt_kr_tr = [vary_all_transcription(i, j, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in kr_range for i in rt_range];
p_rt_kr_tr = plotting2(rt_range, v_rt_kr, "Varying Rt with kr", "Rt", "\$RtcAB\$", kr_labels);

v_rtcr_L_tr = [vary_all_transcription(rt, kr, n, j, c, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in L_range for i in rtcr_range];
p_rtcr_L_tr = plotting2(rtcr_range, v_rtcr_L, "Varying RtcR with L", "RtcR", "\$RtcAB\$", L_labels);

v_rt_L_tr = [vary_all_transcription(i, kr, n, j, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in L_range for i in rt_range];
p_rt_L_tr = plotting2(rt_range, v_rt_L, "Varying Rt with L", "Rt", "\$RtcAB\$", L_labels);

v_rtcr_sigma_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, j, atp, w_rtc, theta_rtc) for j in sigma_range for i in rtcr_range];
p_rtcr_sigma_tr = plotting2(rtcr_range, v_rtcr_sigma, "Varying RtcR with \\sigma", "RtcR", "\$RtcAB\$", sigma_labels);

v_rt_sigma_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, j, atp, w_rtc, theta_rtc) for j in sigma_range for i in rt_range];
p_rt_sigma_tr = plotting2(rt_range, v_rt_sigma, "Varying Rt with \\sigma", "Rt", "\$RtcAB\$", sigma_labels);

v_rtcr_atp_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, j, w_rtc, theta_rtc) for j in atp_range for i in rtcr_range];
p_rtcr_atp_tr = plotting2(rtcr_range, v_rtcr_atp, "Varying RtcR with ATP", "RtcR", "\$RtcAB\$", atp_labels);

v_rt_atp_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, j, w_rtc, theta_rtc) for j in atp_range for i in rt_range];
p_rt_atp_tr = plotting2(rt_range, v_rt_atp, "Varying Rt with ATP", "Rt", "\$RtcAB\$", atp_labels);


v_rtcr_rt_tr = [vary_all_transcription(j, kr, n, L, c, i, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in rt_range1 for i in rtcr_range];
p_rtcr_rt_tr = plotting2(rtcr_range, v_rtcr_rt, "Varying RtcR with Rt", "RtcR", "\$RtcAB\$", rt_labels);

v_rt_rtcr_tr = [vary_all_transcription(i, kr, n, L, c, j, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in rtcr_range1 for i in rt_range];
p_rt_rtcr_tr = plotting2(rt_range, v_rt_rtcr, "Varying Rt with RtcR", "Rt", "\$RtcAB\$", rtcr_labels);


v_rtcr_k1_tr = [vary_all_transcription(rt, kr, n, L, c, i, j, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rtcr_range];
p_rtcr_k1_tr = plotting2(rtcr_range, v_rtcr_k1, "Varying RtcR with k1", "RtcR", "\$RtcAB\$", k_labels);

v_rt_k1_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, j, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rt_range];
p_rt_k1_tr = plotting2(rt_range, v_rt_k1, "Varying Rt with k1", "Rt", "\$RtcAB\$", k_labels);

v_rtcr_k2_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, j, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rtcr_range];
p_rtcr_k2_tr = plotting2(rtcr_range, v_rtcr_k2, "Varying RtcR with k2", "RtcR", "\$RtcAB\$", k_labels);

v_rt_k2_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, j, k3, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rt_range];
p_rt_k2_tr = plotting2(rt_range, v_rt_k2, "Varying Rt with k2", "Rt", "\$RtcAB\$", k_labels);

v_rtcr_k3_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, j, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rtcr_range];
p_rtcr_k3_tr = plotting2(rtcr_range, v_rtcr_k3, "Varying RtcR with k3", "RtcR", "\$RtcAB\$", k_labels);

v_rt_k3_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, j, k4, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rt_range];
p_rt_k3_tr = plotting2(rt_range, v_rt_k3, "Varying Rt with k3", "Rt", "\$RtcAB\$", k_labels);

v_rtcr_k4_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, j, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rtcr_range];
p_rtcr_k4_tr = plotting2(rtcr_range, v_rtcr_k4, "Varying RtcR with k4", "RtcR", "\$RtcAB\$", k_labels);

v_rt_k4_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, j, k5, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rt_range];
p_rt_k4_tr = plotting2(rt_range, v_rt_k4, "Varying Rt with k4", "Rt", "\$RtcAB\$", k_labels);

v_rtcr_k5_tr = [vary_all_transcription(rt, kr, n, L, c, i, k1, k2, k3, k4, j, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rtcr_range];
p_rtcr_k5_tr = plotting2(rtcr_range, v_rtcr_k5, "Varying RtcR with k5", "RtcR", "\$RtcAB\$", k_labels);

v_rt_k5_tr = [vary_all_transcription(i, kr, n, L, c, rtcr, k1, k2, k3, k4, j, sigma, atp, w_rtc, theta_rtc) for j in k_range for i in rt_range];
p_rt_k5_tr = plotting2(rt_range, v_rt_k5, "Varying Rt with k5", "Rt", "\$RtcAB\$", k_labels);


plot(p_rtcr_tr, p_rt_tr, p_rtcr_wrtc, p_rt_wrtc, p_rtcr_trtc, p_rt_trtc, layout=(3,2), size=(1000, 800), margin=2mm)

plot(p_rtcr_tr, p_rt_tr, p_kr_tr, p_L_tr, p_c_tr, p_atp_tr, p_sigma_tr, p_wrtc_tr, layout=(4,2), size=(1000, 800), margin=2mm)

plot(p_trtc_tr, p_k1_tr, p_k2_tr, p_k3_tr, p_k4_tr, p_k5_tr, layout=(3,2), size=(1000, 800), margin=2mm)

plot(p_rtcr_c_tr, p_rtcr_kr_tr, p_rtcr_L_tr, p_rtcr_sigma_tr, p_rtcr_atp_tr, p_rtcr_rt_tr, p_rtcr_wrtc, p_rtcr_trtc, layout=(4,2), size=(1000, 800), margin=2mm)

plot(p_rt_c_tr, p_rt_kr_tr, p_rt_L_tr, p_rt_sigma_tr, p_rt_atp_tr, p_rt_rtcr_tr, p_rt_wrtc, p_rt_trtc, layout=(4,2), size=(1000, 800), margin=2mm)

plot(p_rtcr_tr, p_rtcr_k1_tr, p_rtcr_k2_tr, p_rtcr_k3_tr, p_rtcr_k4_tr, p_rtcr_k5_tr, layout=(3,2), size=(1000,800), margin=2mm)

plot(p_rt_tr, p_rt_k1_tr, p_rt_k2_tr, p_rt_k3_tr, p_rt_k4_tr, p_rt_k5_tr, layout=(3,2), size=(1000,800), margin=2mm)


# checking fraction of active RtcR 
include("parameters_init_rtcrab_expression.jl");

rt_fa = [vary_all_fraction_active_rtcr(i, kr, L, c) for i in rt_range];
p_rt_fa = plotting1(rt_range, rt_fa, "Varying Rt on fraction of active RtcR", "Rt", "Fraction active RtcR");

kr_fa = [vary_all_fraction_active_rtcr(rt, i, L, c) for i in kr_range];
p_kr_fa = plotting1(kr_range, kr_fa, "Varying kr on fraction of active RtcR", "kr", "Fraction active RtcR");

L_fa = [vary_all_fraction_active_rtcr(rt, kr, i, c) for i in L_range];
p_L_fa = plotting1(L_range, L_fa, "Varying L on fraction of active RtcR", "L", "Fraction active RtcR");

c_fa = [vary_all_fraction_active_rtcr(rt, kr, L, i) for i in c_range];
p_c_fa = plotting1(c_range, c_fa, "Varying c on fraction of active RtcR", "c", "Fraction active RtcR");

plot(p_rt_fa, p_kr_fa, p_L_fa, p_c_fa, size=(1000,800))

# varying ATP and checking the Km value for using the correct parameters in the ODE model 
# plotlyjs()
atp1_range = range(0, 20000, length=1000);

oc =  [vary_all(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, i) for i in atp1_range];
# oc = var_atp(atp1_range)
plotting1(atp1_range, oc, "ATP on \\sigma* formation", "ATP", "Open complex")

# oc1 = var_atp_sigma(atp1_range, sigma_range)
oc1 = [vary_all(i, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, j, i) for j in sigma_range for i in atp1_range];

oc1 = reshape(oc1, (1000,5));
oc1 = [d[:] for d in eachcol(oc1)];

plot(atp1_range, oc1, title="ATP and \\sigma on \\sigma* formation", xlabel="ATP", ylabel="Open complex", labels=sigma_labels)


oc1_sigma0_hm = oc1[1][end]/2
oc1_sigma1_hm = oc1[2][end]/2
oc1_sigma5_hm = oc1[3][end]/2
oc1_sigma10_hm = oc1[4][end]/2
oc1_sigma20_hm = oc1[5][end]/2

itp = LinearInterpolation(oc1[1], atp1_range)
itp(oc1_sigma0_hm)

itp2 = LinearInterpolation(oc1[2], atp1_range)
itp2(oc1_sigma1_hm)

itp3 = LinearInterpolation(oc1[3], atp1_range)
itp3(oc1_sigma5_hm)

itp4 = LinearInterpolation(oc1[4], atp1_range)
itp4(oc1_sigma10_hm)

itp5 = LinearInterpolation(oc1[5], atp1_range)
itp5(oc1_sigma20_hm)

