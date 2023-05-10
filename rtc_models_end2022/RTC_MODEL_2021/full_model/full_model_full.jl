using Plots, DifferentialEquations, DataFrames, LSODA, Statistics
include("fullmodel_functions.jl")
include("fullmodel_params.jl")

# plotlyjs()

sol_ode = simple_solve!(rtc_system_fullodes!, init_ode, tspan, params)

plot(sol_ode, labels=labels_ode, palette=:tab10)
plot(sol_ode[2:end], yaxis=(:log10, (1,Inf)), labels=labels, palette=:tab10)


# full ode model 
rm_a_ode, rtca_ode, rm_b_ode, rtcb_ode, rm_r_ode, rtcr_ode, ribo_d_ode, ribo_h_ode, ribo_t_ode, alpha_ode, fa_ode, ra_ode, v_ode, tscr_ode, tlr_ode, tscr_b_ode, tlr_b_ode, rtcr_tscr_ode, rtcr_tlr_ode, rdrtca_ode, rtrtcb_ode = split(species_ode, sol_ode)

p_rtc_ode = plot(sol_ode.t, xlabel="t", [rtca_ode, rtcb_ode, rtcr_ode], labels=["RtcA" "RtcB" "RtcR"]);
p_mrna_ode = plot(sol_ode.t, xlabel="t", [rm_a_ode, rm_b_ode, rm_r_ode], labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"]);
p_ribo_ode = plot(sol_ode.t, xlabel="t", [ribo_d_ode, ribo_h_ode, ribo_t_ode], labels=["Rd" "Rh" "Rt"], title="Ribosomes")

p_fartcr_ode = plot(sol_ode.t, xlabel="t", fa_ode, title="Fraction of active RtcR", legend=false);
p_rartcr_ode = plot(sol_ode.t, xlabel="t", ra_ode, title="Active RtcR", legend=false);
p_trsc_i_ode = plot(sol_ode.t, xlabel="t", v_ode, title="Rate of transcription initiation", legend=false);
p_trsc_ode = plot(sol_ode.t, xlabel="t", [tscr_ode, tscr_b_ode], title="Transcription rates", labels=["RtcA" "RtcB"]);
p_trsl_ode = plot(sol_ode.t, xlabel="t", [tlr_ode, tlr_b_ode, rtcr_tlr_ode], title="Translation rates", labels=["RtcA" "RtcB" "RtcR"];)
p_int_ode = plot(sol_ode.t, xlabel="t", [rdrtca_ode, rtrtcb_ode], title="Intermediates", labels=["RdRtcA" "RtRtcB"]);

p_sol_ode = plot(sol_ode);

plot(p_sol_ode, p_rtc_ode, p_mrna_ode, p_ribo_ode, p_fartcr_ode, p_rartcr_ode, p_trsc_ode, p_trsc_i_ode, p_trsl_ode, p_int_ode, size=(800,1000))

plot(sol.t, ribo_d, xlabel="t", title="Damaged ribosomes", legend=false)
plot(sol.t, ribo_t, xlabel="t", title="Tagged ribosomes", legend=false)
plot(sol.t, ribo_h, xlabel="t", title="Healthy ribosomes", legend=false)


# checking difference in full ode system and using mm laws to get rid of intermediate odes
function checking_difference_with_mm_and_odes()
    all = [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb]

    all_max = []
    for i in all
        push!(all_max, (maximum(i)))
    end

    println(typeof(all_max)) 

    all_ode = [rm_a_ode, rtca_ode, rm_b_ode, rtcb_ode, rm_r_ode, rtcr_ode, ribo_d_ode, ribo_h_ode, ribo_t_ode, alpha_ode, fa_ode, ra_ode, v_ode, tscr_ode, tlr_ode, tscr_b_ode, tlr_b_ode, rtcr_tscr_ode, rtcr_tlr_ode, rdrtca_ode, rtrtcb_ode]
    all_max_ode = []
    for j in all_ode
        push!(all_max_ode, (maximum(j)))
    end

    println(all_max_ode)

    dif = all_max_ode.-all_max

    names = ["rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "ribo_d", "ribo_h", "ribo_t", "alpha", "fa", "ra", "v", "tscr", "tlr", "tscr_b", "tlr_b", "rtcr_tscr", "rtcr_tlr", "rdrtca", "rtrtcb"]

    df = DataFrame(species=names, all_max=all_max, all_max_ode=all_max_ode, difference=dif)
end
println(checking_difference_with_mm_and_odes())

# varying parameters

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
atp_range = range(0, 10, length=6) 
plot1_atp, plot2_atp = plots_for_param("ATP", params, 5, atp_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_atp
plot2_atp

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
kr_range = range(0, 10, length=6) 
plot1_kr, plot2_kr = plots_for_param("kr", params, 1, kr_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_kr
plot2_kr

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
L_range = range(0, 10, length=6) 
plot1_L, plot2_L = plots_for_param("L", params, 2, L_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_L
plot2_L

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
c_range = range(0, 5, length=6)
plot1_c, plot2_c = plots_for_param("c", params, 3, c_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_c
plot2_c

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
sigma_range = range(0, 10, length=6)
plot1_sigma, plot2_sigma = plots_for_param("\\sigma", params, 4, sigma_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_sigma
plot2_sigma

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k1_range = range(0, 10, length=6)
plot1_k1, plot2_k1 = plots_for_param("k1", params, 6, k1_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k1
plot2_k1

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k2_range = range(0, 10, length=6)
plot1_k2, plot2_k2 = plots_for_param("k2", params, 7, k2_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k2
plot2_k2

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k3_range = range(0, 10, length=6)
plot1_k3, plot2_k3 = plots_for_param("k3", params, 8, k3_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k3
plot2_k3

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k4_range = range(0, 10, length=6)
plot1_k4, plot2_k4 = plots_for_param("k4", params, 9, k4_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k4
plot2_k4

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k5_range = range(0, 10, length=6)
plot1_k5, plot2_k5 = plots_for_param("k5", params, 10, k5_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k5
plot2_k5

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
wrtc_range = range(0, 10, length=6)
plot1_wrtc, plot2_wrtc = plots_for_param("w_rtc", params, 11, wrtc_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_wrtc
plot2_wrtc

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
thrtc_range = range(0, 10, length=6)
plot1_thrtc, plot2_thrtc = plots_for_param("thrtc", params, 12, thrtc_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_thrtc
plot2_thrtc

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
max_range = range(0, 10, length=6)
plot1_max, plot2_max = plots_for_param("max", params, 13, max_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_max
plot2_max

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
thr_range = range(0, 10, length=6)
plot1_thr, plot2_thr = plots_for_param("thr", params, 14, thr_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_thr
plot2_thr

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
d_range = range(0, 10, length=6)
plot1_d, plot2_d = plots_for_param("d", params, 16, d_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_d
plot2_d

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k1a_range = range(0, 10, length=6)
plot1_k1a, plot2_k1a = plots_for_param("k1a", params, 17, k1a_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k1a
plot2_k1a

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k2a_range = range(0, 10, length=6)
plot1_k2a, plot2_k2a = plots_for_param("k2a", params, 18, k2a_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k2a
plot2_k2a

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k3a_range = range(0, 10, length=6)
plot1_k3a, plot2_k3a = plots_for_param("k3a", params, 19, k3a_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k3a
plot2_k3a

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
kab_range = range(0, 10, length=6)
plot1_kab, plot2_kab = plots_for_param("kab", params, 20, kab_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_kab
plot2_kab

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
kbb_range = range(0, 10, length=6)
plot1_kbb, plot2_kbb = plots_for_param("kbb", params, 21, kbb_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_kbb
plot2_kbb

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
kcb_range = range(0, 10, length=6)
plot1_kcb, plot2_kcb = plots_for_param("kcb", params, 22, kcb_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_kcb
plot2_kcb

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
k_range = range(0, 4, length=6)
plot1_k, plot2_k = plots_for_param("k", params, 23, k_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_k
plot2_k

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
rh_range = range(0, 20, length=6) # over 10 and start seeing negative values?
plot1_rh, plot2_rh = plots_for_param("Rh", init_ode, 8, rh_range, rtc_system!, init_ode, tspan, params, species_ode)
plot1_rh
plot2_rh

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
rtca_range = range(0, 20, length=6) # over 10 and start seeing negative values?
plot1_rtca, plot2_rtca = plots_for_param("RtcA", init_ode, 2, rtca_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_rtca
plot2_rtca

params = [kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init_ode = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0, rdrtca_0, rtrtcb_0];
rtcb_range = range(0, 20, length=6) # over 10 and start seeing negative values?
plot1_rtcb, plot2_rtcb = plots_for_param("RtcB", init_ode, 4, rtcb_range, rtc_system_fullodes!, init_ode, tspan, params, species_ode)
plot1_rtcb
plot2_rtcb






atp_long = range(0, 100, length=100)

rm_a_res1, rtca_res1, rm_b_res1, rtcb_res1, rm_r_res1, rtcr_res1, ribo_d_res1, ribo_h_res1, ribo_t_res1, rdrtca_res1, rtrtcb_res1 = vary_param_xaxis(atp_long, params, 5, species)

plot(atp_long, [rtca_res, rtcb_res])
plot(atp_long, [rtcr_res])
plot(atp_long, [ribo_d_res, ribo_h_res, ribo_t_res])
plot(atp_long, ribo_d_res)

