using DifferentialEquations, Plots, DataFrames, Measures, Dierckx

include("MWC_functions.jl");
include("parameters_init_rtcrab_expression.jl");
# plotlyjs()

# res_rm_r, res_rtc_r, res_rm_r_r, res_rtcr_r = vary_param_for_rtcRab_model!(r_range, 16, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr)
# res_rm_atp, res_rtc_atp, res_rm_r_atp, res_rtcr_atp = vary_param_for_rtcRab_model!(atp_range, 6, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr)

# plot(r_range, res_rtc_r)
# plot(atp_range, res_rtc_atp)
# plot(r_range, [res_rtc_r, res_rtc_atp], ylim=(1e5,1e10), yscale=(:log10), legend=:right)
# plot(r_range, [res_rtc_r, res_rtc_atp], legend=:right)


sol = simple_solve!(rtc_expression!, init, tspan, params);
plot(sol, labels=["mRNA_AB" "RtcAB"])

sol_rtcr = simple_solve!(rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr);
plot(sol_rtcr, labels=["mRNA_AB" "RtcAB" "mRNA_R" "RtcR"], title="Model solution")
plot(sol_rtcr, ylim=(1,1e6), yscale=(:log10), ylabel=("log(concentration)"), labels=["mRNA_AB" "RtcAB" "mRNA_R" "RtcR"], title="Model solution")


# varying parameters for model for expression of rtcab 
vary_rtcr = vary_param_for_rtcab_model_plot!(rtcr_range, 5, rtc_expression!, init, tspan, params, "RtcR", "RtcBA expression");
vary_rt = vary_param_for_rtcab_model_plot!(rt_range, 1, rtc_expression!, init, tspan, params, "Rt", "RtcBA expression");
vary_sigma = vary_param_for_rtcab_model_plot!(sigma_range, 6, rtc_expression!, init, tspan, params, "\\sigma", "RtcBA expression");
vary_atp = vary_param_for_rtcab_model_plot!(atp_range, 7, rtc_expression!, init, tspan, params, "ATP", "RtcBA expression");
vary_kr = vary_param_for_rtcab_model_plot!(kr_range, 2, rtc_expression!, init, tspan, params, "kr", "RtcBA expression");
vary_L = vary_param_for_rtcab_model_plot!(L_range, 3, rtc_expression!, init, tspan, params, "L", "RtcBA expression");
vary_c = vary_param_for_rtcab_model_plot!(c_range, 4, rtc_expression!, init, tspan, params, "c", "RtcBA expression");
vary_wrtc = vary_param_for_rtcab_model_plot!(wrtc_range, 13, rtc_expression!, init, tspan, params, "Max transcription rate of \$RtcBA\$", "RtcBA expression");
vary_thetartc = vary_param_for_rtcab_model_plot!(thetartc_range, 14, rtc_expression!, init, tspan, params, "Transcription threshold for \$RtcBA\$", "RtcBA expression");
vary_max = vary_param_for_rtcab_model_plot!(max_range, 15, rtc_expression!, init, tspan, params, "Max translation rate of RtcBA", "RtcBA expression");
vary_thr = vary_param_for_rtcab_model_plot!(thr_range, 16, rtc_expression!, init, tspan, params, "Translation threshold for RtcBA", "RtcBA expression");
vary_r = vary_param_for_rtcab_model_plot!(r_range, 17, rtc_expression!, init, tspan, params, "r", "RtcBA expression");

plot(vary_rtcr, vary_rt, vary_sigma, vary_atp, vary_kr, vary_L, vary_c, vary_wrtc, vary_thetartc, vary_max, vary_thr, vary_r, layout=(4,3), size=(1000,800), margin=3mm)

#rtcr addition

vary_rt_full = vary_param_for_rtcRab_plot!(rt_range, 1, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "Rt", "RtcBA expression");
vary_sigma_full = vary_param_for_rtcRab_plot!(sigma_range, 5, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "\\sigma", "RtcBA expression");
vary_atp_full = vary_param_for_rtcRab_plot!(atp_range, 6, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "ATP", "RtcBA expression");
vary_kr_full = vary_param_for_rtcRab_plot!(kr_range, 2, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "kr", "RtcBA expression");
vary_L_full = vary_param_for_rtcRab_plot!(L_range, 3, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "L", "RtcBA expression");
vary_c_full = vary_param_for_rtcRab_plot!(c_range, 4, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "c", "RtcBA expression");
vary_wrtc_full = vary_param_for_rtcRab_plot!(wrtc_range, 12, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "Max transcription rate of \$RtcBA\$", "RtcBA expression");
vary_thetartc_full = vary_param_for_rtcRab_plot!(thetartc_range, 13, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "Transcription threshold for \$RtcBA\$", "RtcBA expression");
vary_max_full = vary_param_for_rtcRab_plot!(max_range, 14, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "Max translation rate of RtcBA", "RtcBA expression");
vary_thr_full = vary_param_for_rtcRab_plot!(thr_range, 15, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "Translation threshold for RtcBA", "RtcBA expression");
vary_r_full = vary_param_for_rtcRab_plot!(r_range, 16, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "r", "RtcBA expression");

plot(vary_rt_full, vary_rt_full, vary_sigma_full, vary_atp_full, vary_kr_full, vary_L_full, vary_c_full, vary_wrtc_full, vary_thetartc_full, vary_max_full, vary_thr_full, vary_r_full, layout=(4,3), size=(1000,800), margin=3mm, legend=false)


vary_k1_full = vary_param_for_rtcRab_plot!(k_range, 7, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "k1", "RtcBA expression");
vary_k2_full = vary_param_for_rtcRab_plot!(k_range, 8, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "k2", "RtcBA expression");
vary_k3_full = vary_param_for_rtcRab_plot!(k_range, 9, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "k3", "RtcBA expression");
vary_k4_full = vary_param_for_rtcRab_plot!(k_range, 10, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "k4", "RtcBA expression");
vary_k5_full = vary_param_for_rtcRab_plot!(k_range, 11, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr, "k5", "RtcBA expression");

plot(vary_k1_full, vary_k2_full, vary_k3_full, vary_k4_full, vary_k5_full, layout=(5,1), size=(600,1000), margin=5mm)

# ATP
atp1_range = range(0, 100, length=5)
atp_range1 = range(0.1, 100, length=1000);

p_intercept_plot_atp = log_intercept_plot!(atp_range1, 6, "ATP", "0.01-10");
p_vary_params_solutions_atp, p_vary_params_species_atp, p_vary_param_long_sol_atp, p_vary_param_long_sol_log_atp = param_explore_plots!(atp1_range, 6, "ATP", "0.01-100", atp_range1);

plot(p_vary_params_solutions_atp)
plot(p_vary_params_species_atp)
plot(p_intercept_plot_atp)
plot(p_vary_param_long_sol_atp)
plot(p_vary_param_long_sol_log_atp)



# Rt 
rt1_range = range(0, 50, length=5)  
rt_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_rt = log_intercept_plot!(rt_range1, 1, "Rt", "0.01-100");
p_vary_params_solutions_rt, p_vary_params_species_rt, p_vary_param_long_sol_rt, p_vary_param_long_sol_log_rt = param_explore_plots!(rt1_range, 1, "Rt", "0.01-100", rt_range1);

plot(p_vary_params_solutions_rt)
plot(p_vary_params_species_rt)
# plot(p_intercept_plot_rt)
plot(p_vary_param_long_sol_rt)
plot(p_vary_param_long_sol_log_rt)


# kr 
kr1_range = range(10, 1000, length=5)
kr_range1 = range(0.01, 1000, length=10000);

# p_intercept_plot_kr = log_intercept_plot!(rt_range1, 2, "Rt", "0.01-100");
p_vary_params_solutions_kr, p_vary_params_species_kr, p_vary_param_long_sol_kr, p_vary_param_long_sol_log_kr = param_explore_plots!(kr1_range, 2, "kr", "0.01-1000", kr_range1);

plot(p_vary_params_solutions_kr)
plot(p_vary_params_species_kr)
# plot(p_intercept_plot_kr)
plot(p_vary_param_long_sol_kr)
plot(p_vary_param_long_sol_log_kr)


# sigma 
sigma1_range = range(0, 1, length=5)
sigma_range1 = range(0.01, 10, length=1000);

p_intercept_plot_sigma = log_intercept_plot!(sigma_range1, 5, "\\sigma", "0.01-10");
p_vary_params_solutions_sigma, p_vary_params_species_sigma, p_vary_param_long_sol_sigma, p_vary_param_long_sol_log_sigma = param_explore_plots!(sigma1_range, 5, "\\sigma", "0.01-10", sigma_range1);

plot(p_vary_params_solutions_sigma)
plot(p_vary_params_species_sigma)
plot(p_intercept_plot_sigma)
plot(p_vary_param_long_sol_sigma)
plot(p_vary_param_long_sol_log_sigma)


# r
r1_range = range(0, 10, length=5)
r_range1 = range(0.01, 10, length=1000);

p_intercept_plot_r = log_intercept_plot!(r_range1, 16, "r", "0.01-10");
p_vary_params_solutions_r, p_vary_params_species_r, p_vary_param_long_sol_r, p_vary_param_long_sol_log_r = param_explore_plots!(r1_range, 16, "r", "0.01-100", r_range1);

plot(p_vary_params_solutions_r)
plot(p_vary_params_species_r)
plot(p_intercept_plot_r)
plot(p_vary_param_long_sol_r)
plot(p_vary_param_long_sol_log_r)

# L
L1_range = range(0, 100, length=5)
L_range1 = range(0.01, 1000, length=1000);

# p_intercept_plot_L = log_intercept_plot!(L_range1, 3, "L", "0.01-100");
p_vary_params_solutions_L, p_vary_params_species_L, p_vary_param_long_sol_L, p_vary_param_long_sol_log_L = param_explore_plots!(L1_range, 3, "L", "0.01-1000", L_range1);

plot(p_vary_params_solutions_L)
plot(p_vary_params_species_L)
# plot(p_intercept_plot_L)
plot(p_vary_param_long_sol_L)
plot(p_vary_param_long_sol_log_L)


# c
c1_range = range(0, 1, length=5)
c_range1 = range(0.01, 10, length=1000);

p_intercept_plot_c = log_intercept_plot!(c_range1, 4, "c", "0.01-10");
p_vary_params_solutions_c, p_vary_params_species_c, p_vary_param_long_sol_c, p_vary_param_long_sol_log_c = param_explore_plots!(c1_range, 4, "c", "0.01-10", c_range1);

plot(p_vary_params_solutions_c)
plot(p_vary_params_species_c)
plot(p_intercept_plot_c)
plot(p_vary_param_long_sol_c)
plot(p_vary_param_long_sol_log_c)


# max transcription rate
wrtc1_range = range(0, 0.1, length=5)
wrtc_range1 = range(0.01, 10, length=1000);

p_intercept_plot_wrtc = log_intercept_plot!(wrtc_range1, 12, "wrtc", "0.01-10");
p_vary_params_solutions_wrtc, p_vary_params_species_wrtc, p_vary_param_long_sol_wrtc, p_vary_param_long_sol_log_wrtc = param_explore_plots!(wrtc1_range, 12, "wrtc", "0.01-10", wrtc_range1);

plot(p_vary_params_solutions_wrtc)
plot(p_vary_params_species_wrtc)
plot(p_intercept_plot_wrtc)
plot(p_vary_param_long_sol_wrtc)
plot(p_vary_param_long_sol_log_wrtc)


# transcription threshold 
thrtc1_range = range(0, 100, length=5)
thrtc_range1 = range(0.01, 100, length=1000);

# p_intercept_plot_thrtc = log_intercept_plot!(thrtc_range1, 13, "thrtc", "0.01-100");
p_vary_params_solutions_thrtc, p_vary_params_species_thrtc, p_vary_param_long_sol_thrtc, p_vary_param_long_sol_log_thrtc = param_explore_plots!(thrtc1_range, 13, "thrtc", "0.01-100", thrtc_range1);

plot(p_vary_params_solutions_thrtc)
plot(p_vary_params_species_thrtc)
# plot(p_intercept_plot_thrtc)
plot(p_vary_param_long_sol_thrtc)
plot(p_vary_param_long_sol_log_thrtc)


# max translation rate 
maxtl1_range = range(0, 100, length=5)
maxtl_range1 = range(0.01, 10, length=1000);

p_intercept_plot_maxtl = log_intercept_plot!(maxtl_range1, 14, "maxtl", "0.01-10");
p_vary_params_solutions_maxtl, p_vary_params_species_maxtl, p_vary_param_long_sol_maxtl, p_vary_param_long_sol_log_maxtl = param_explore_plots!(maxtl1_range, 14, "maxtl", "0.01-10", maxtl_range1);

plot(p_vary_params_solutions_maxtl)
plot(p_vary_params_species_maxtl)
plot(p_intercept_plot_maxtl)
plot(p_vary_param_long_sol_maxtl)
plot(p_vary_param_long_sol_log_maxtl)


# translation threshold 
thrtl1_range = range(0, 100, length=5)
thrtl_range1 = range(0.01, 1000, length=1000);

# p_intercept_plot_thrtl = log_intercept_plot!(thrtl_range1, 15, "thrtl", "0.01-100");
p_vary_params_solutions_thrtl, p_vary_params_species_thrtl, p_vary_param_long_sol_thrtl, p_vary_param_long_sol_log_thrtl = param_explore_plots!(thrtl1_range, 15, "thrtl", "0.01-100", thrtl_range1);

plot(p_vary_params_solutions_thrtl)
plot(p_vary_params_species_thrtl)
# plot(p_intercept_plot_thrtl)
plot(p_vary_param_long_sol_thrtl)
plot(p_vary_param_long_sol_log_thrtl)

# k1
k11_range = range(0, 100, length=5)
k1_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_k1 = log_intercept_plot!(k1_range1, 7, "k1", "0.01-100");
p_vary_params_solutions_k1, p_vary_params_species_k1, p_vary_param_long_sol_k1, p_vary_param_long_sol_log_k1 = param_explore_plots!(k11_range, 7, "k1", "0.01-100", k1_range1);

plot(p_vary_params_solutions_k1)
plot(p_vary_params_species_k1)
# plot(p_intercept_plot_k1)
plot(p_vary_param_long_sol_k1)
plot(p_vary_param_long_sol_log_k1)


#k2
k21_range = range(0, 100, length=5)
k2_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_k2 = log_intercept_plot!(k2_range1, 8, "k2", "0.01-100");
p_vary_params_solutions_k2, p_vary_params_species_k2, p_vary_param_long_sol_k2, p_vary_param_long_sol_log_k2 = param_explore_plots!(k21_range, 8, "k2", "0.01-100", k2_range1);

plot(p_vary_params_solutions_k2)
plot(p_vary_params_species_k2)
# plot(p_intercept_plot_k2)
plot(p_vary_param_long_sol_k2)
plot(p_vary_param_long_sol_log_k2)


#k3 
k31_range = range(0, 100, length=5)
k3_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_k3 = log_intercept_plot!(k3_range1, 9, "k3", "0.01-100");
p_vary_params_solutions_k3, p_vary_params_species_k3, p_vary_param_long_sol_k3, p_vary_param_long_sol_log_k3 = param_explore_plots!(k31_range, 9, "k3", "0.01-100", k3_range1);

plot(p_vary_params_solutions_k3)
plot(p_vary_params_species_k3)
# plot(p_intercept_plot_k3)
plot(p_vary_param_long_sol_k3)
plot(p_vary_param_long_sol_log_k3)


#k4
k41_range = range(0, 100, length=5)
k4_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_k4 = log_intercept_plot!(k4_range1, 10, "k4", "0.01-100");
p_vary_params_solutions_k4, p_vary_params_species_k4, p_vary_param_long_sol_k4, p_vary_param_long_sol_log_k4 = param_explore_plots!(k41_range, 10, "k4", "0.01-100", k4_range1);

plot(p_vary_params_solutions_k4)
plot(p_vary_params_species_k4)
# plot(p_intercept_plot_k4)
plot(p_vary_param_long_sol_k4)
plot(p_vary_param_long_sol_log_k4)


#k5 
k51_range = range(0, 100, length=5)
k5_range1 = range(0.1, 100, length=1000);

# p_intercept_plot_k5 = log_intercept_plot!(k5_range1, 11, "k5", "0.01-100");
p_vary_params_solutions_k5, p_vary_params_species_k5, p_vary_param_long_sol_k5, p_vary_param_long_sol_log_k5 = param_explore_plots!(k51_range, 11, "k5", "0.01-100", k5_range1);

plot(p_vary_params_solutions_k5)
plot(p_vary_params_species_k5)
# plot(p_intercept_plot_k5)
plot(p_vary_param_long_sol_k5)
plot(p_vary_param_long_sol_log_k5)





# all equations solved

sol_rtcr = simple_solve!(rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr);
plot(sol_rtcr)
solDF = DataFrame([[j[i] for j in sol_rtcr.u] for i=1:length(sol_rtcr.u[1])], species_rtcr);
rtcR = solDF[:, :rtcr]
mrnaAB = solDF[:, :rm]
rtcAB = solDF[:, :rtc]
mrnaR = solDF[:, :rm_r]

rtcr_plot = plot(sol_rtcr.t, rtcR, ylabel="RtcR", xlabel="t");

n = 6
alpha = rt/kr
fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
ra = fa*rtcR # amount of active RtcR
v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
tlel = max*atp/(thr+atp) # translation elongation rate
tlr = r*mrnaAB*tlel # translation rate so formation of Rtc 

# rtcr 
rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
rtcr_tlel = max*atp/(thr+atp)
rtcr_tlr = r*mrnaR*rtcr_tlel


plot(sol_rtcr.t, [rtcR, mrnaAB, rtcAB, mrnaR, ra, v, tscr, tlr, rtcr_tlr], ylim=(0,2500), ylabel="rate/conc.", xlabel="t", labels=["RtcR" "AB-mRNA" "RtcAB" "R-mRNA" "Ra" "v" "tscr" "tlr" "rtcr_tlr"], palette=:tab10)
plot(sol_rtcr.t[2:end], [rtcR[2:end], mrnaAB[2:end], rtcAB[2:end], mrnaR[2:end], ra[2:end], v[2:end], tscr[2:end], tlr[2:end], rtcr_tlr[2:end]], yaxis=(:log10, (1,Inf)), ylabel="log(rate/conc.)", xlabel="t", labels=["RtcR" "AB-mRNA" "RtcAB" "R-mRNA" "Ra" "v" "tscr" "tlr" "rtcr_tlr"], palette=:tab10)






