function set_time_vars(file)
    csv = DataFrame(CSV.File(file))

    csv.atp = csv.atp/5

    kin = @. csv.gr*11.29/(2.0923*csv.atp/(255.73+csv.atp))

    lam_t = QuadraticInterpolation(csv_atp."gr",csv_atp."t")
    atp_t = QuadraticInterpolation(csv_atp."atp",csv_atp."t")

    # set kin curve  
    # rh = 11.29; g_max = 2.0923;
    # kin_model = @. lam_t*11.29/2.0923
    # kin_model = @. lam_t*11.29/(2.0923*atp_t/(255.73+atp_t))
    kin_t = QuadraticInterpolation(kin, csv_atp."t") 

    return csv_atp.t, atp_t, lam_t, kin_t
end