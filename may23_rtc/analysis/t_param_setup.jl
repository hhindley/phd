csv_atp = DataFrame(CSV.File("/home/holliehindley/phd/data/atp_for_rtcmodel.csv"))

csv_atp.atp = csv_atp.atp/5

lam_t = QuadraticInterpolation(csv_atp."gr",csv_atp."t")
atp_t = QuadraticInterpolation(csv_atp."atp",csv_atp."t")

# p_atp = plot(scatter(x=csv_atp."t", y=atp_t), Layout(xaxis_type="", xaxis_range=(0,1500)));
# p_lam = plot(scatter(x=csv_atp."t", y=lam_t), Layout(xaxis_type="", xaxis_range=(0,1500)));

rh = 11.29
kin_model = @. lam_t*rh/g_max
kin_t = QuadraticInterpolation(kin_model, csv_atp."t") 
# p_kin = plot(scatter(x=csv_atp."t", y=kin_t), Layout(xaxis_type="", xaxis_range=(0,1500)));

# p = plot_time_vars(lam_t, atp_t, kin_t)
