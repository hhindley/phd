tspan = (0, 200);
global k1 = 0.1; 
global k2 = 0.1;
global k3 = 0.1;
global k4 = 0.1;
global k5 = 10;
global sigma = 10;
global atp = 10;
global rtcr = 10; #10;
global rt = 10; # 10;
global kr = 10;
global kt = 1000;
global L = 100; # t0/r0
global c = 0.01; #kr/kt
global w_rtc = 4.14 # 5 # 10;    # max. transcription rate 
global theta_rtc = 20 # 4.38 # 5000;  # transcription threshold
global max = 4;            # max. translation rate 
global thr = 20 #7 # 5000;        # translation threshold 
global r = 10;
global gr = 0.01*r
global d = 0.01;

rm_0 = 0;
rtc_0 = 0;
rm_r_0 = 0;
rtcr_0 = 0; 

params = [rt, kr, L, c, rtcr, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, r, gr, d];
params_rtcr = [rt, kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, r, gr, d];
init = [rm_0, rtc_0];
init_rtcr = [rm_0, rtc_0, rm_r_0, rtcr_0];

species = [:rm, :rtc];
species_rtcr = [:rm, :rtc, :rm_r, :rtcr];
species_rtcr_labels = ["RtcAB-mRNA" "RtcAB" "RtcR-mRNA" "RtcR"]

rtcr_range = range(0, 100, length=100);
rt_range = range(0, 100, length=100);
sigma_range = range(0, 50, length=100);
atp_range = range(0, 1000, length=100);
kr_range = range(1, 100, length=100);
L_range = exp10.(range(log10(1), log10(50), length=100));
c_range = range(0, 5, length=100);
wrtc_range = range(0, 100, length=100);
thetartc_range = range(0, 3000, length=100);
max_range = range(0, 50, length=100);
thr_range = range(0, 3000, length=100);
r_range = range(0, 50, length=100);
k_range = range(0, 100, length=100);
