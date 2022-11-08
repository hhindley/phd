using Plots, DifferentialEquations, DataFrames
include("MWC_functions.jl")

plotlyjs()


function rtc_system!(dydt, initial, params, t)

    kr, L, c, rtcr, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_hdt = zeros(length(dydt))

    n = 6

    ribo_t = ribo_tot - ribo_d - ribo_h #10 

    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    # ribosomes 
    rdrtca = rtca_0 + rdrtca_0 - rtca #0
    rtrtcb = rtcb_0 + rtrtcb_0 - rtcb #0

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - gr*rtca - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - gr*rtcb - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[8] = kc_b*atp*rtrtcb - k*ribo_h 
end




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
w_rtc = 4.14 # 5 # 10;    # max. transcription rate 
theta_rtc = 20 # 4.38 # 5000;  # transcription threshold
max = 4;            # max. translation rate 
thr = 20 #7 # 5000;        # translation threshold 
# r = 10;
d = 0.01;

#ribosomes
k1_a = 0.1; #0
k2_a = 0.1; #0
k3_a = 0.1; #0
ka_b = 0.1; #0
kb_b = 0.1; #0
kc_b = 1; #0
k = 1; #0

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

ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 #10
gr = 0.01*ribo_tot

params = [kr, L, c, rtcr, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_tot];
init = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, ribo_d_0, ribo_h_0];

tspan = (0,200)

sol = simple_solve!(rtc_system!, init, tspan, params)

labels = ["RtcA-mRNA" "RtcA" "RtcB-mRNA" "RtcB" "RtcR-mRNA" "RtcR" "Ribo_d" "Ribo_h"]
plot(sol, labels=labels, palette=:tab10)

species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :ribo_d, :ribo_h]

function split_species(species)
    solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_h = solDF[:, :ribo_h]

    ribo_t = ribo_tot .- ribo_d .- ribo_h 
    alpha = ribo_t./kr
    fa_1 = (alpha.+1).^6 # fraction of active RtcR
    fa_2 = L*(c.*alpha.+1).^6 +(alpha.+1).^6
    fa = fa_1./fa_2
    ra = fa.*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel

    # ribosomes 
    rdrtca = rtca_0 + rdrtca_0 .- rtca
    rtrtcb = rtcb_0 + rtrtcb_0 .- rtcb
return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_species(species)

println(maximum(rtca))

plot(sol.t, [rtca, rtcb], labels=["RtcA" "RtcB"])
plot(sol.t, [rm_a, rm_b], labels=["RtcA-mRNA" "RtcB-mRNA"])
plot(sol.t, rtcr)
plot(sol.t, [ribo_d, ribo_h, ribo_t], labels=["Rd" "Rh" "Rt"])

plot(sol.t, fa)
plot(sol.t, ra)
plot(sol.t, v)
plot(sol.t, [tscr, tscr_b])
plot(sol.t, [tlr, tlr_b, rtcr_tlr])
plot(sol.t, [rdrtca, rtrtcb])

println(rtca)




