rtca_0 = 0#0.00894; 
rtcb_0 = 0#0.0216; 
rh_0 = 11.29; #69.56; #69.4
rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
rm_a_0 = 0; 
rm_b_0 = 0; 
rm_r_0 = 0#0.0131#0.04 # 0; 
rd_0 = 0; 
rt_0 = 0;

init_rtc = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]

init_inhib = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, 0]


# rtc_init_sig = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, 0]
# species_rtc_sig = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :sig]
