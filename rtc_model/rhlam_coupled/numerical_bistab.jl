include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

kdam_range1 = reverse(kdam_range)
res = numerical_bistability_analysis(test, params_rtc1, init_rtc, :rtca, species_rtc, kdam_range, kdam)
res1 = numerical_bistability_analysis(test, params_rtc1, init_rtc, :rtca, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res), scatter(x=kdam_range1, y=res1)])


res_ss = numerical_bistability_analysis(test, params_rtc1, ssvals_rtc, :rtca, species_rtc, kdam_range, kdam)
res1_ss = numerical_bistability_analysis(test, params_rtc1, ssvals_rtc, :rtca, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res_ss), scatter(x=kdam_range1, y=res1_ss)])
