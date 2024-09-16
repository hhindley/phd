include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

kdam_range = range(0.0, 1.5, length=50)
kdam_range1 = reverse(kdam_range)
res = numerical_bistability_analysis(test, params_rtc1, init_rtc, species_rtc, kdam_range, kdam)
res1 = numerical_bistability_analysis(test, params_rtc1, init_rtc, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res.rtca), scatter(x=kdam_range1, y=res1.rtca)])

plot([scatter(x=kdam_range1, y=res1[!,s], name="$s") for s in species_rtc])


res_ss = numerical_bistability_analysis(test, params_rtc1, ssvals_rtc, :rtca, species_rtc, kdam_range, kdam)
res1_ss = numerical_bistability_analysis(test, params_rtc1, ssvals_rtc, :rtca, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res_ss), scatter(x=kdam_range1, y=res1_ss)])



init_on = [res[1,s] for s in species_rtc]
init_off = [res1[1,s] for s in species_rtc]

# check stability of branches 

param_dict_dam = deepcopy(params_rtc1)

unstable, eigenvals = calc_unstable_points(kdam_range, res, param_dict_dam)
unstable1, eigenvals1 = calc_unstable_points(kdam_range1, res1, param_dict_dam)


p = scatter(x=kdam_range, y=res.rtca)
p1 = scatter(x=kdam_range1, y=res1.rtca)

pu = [scatter(x=[kdam_range1[i]], y=[res1[i,:rtca]], mode="markers") for i in unstable1]

plot([p, p1, pu...])


s = :rm_a

eigs_real_on = eigs_to_df(eigenvals)
eigs_imag_on = eigs_to_df(eigenvals; real1=false)
plot(scatter(x=eigs_real_on[!,s], y=eigs_imag_on[!,s], mode="markers"), Layout(xaxis_title="Real", yaxis_title="Imaginary"))

eigs_real_off = eigs_to_df(eigenvals1)
eigs_imag_off = eigs_to_df(eigenvals1; real1=false)
plot(scatter(x=eigs_real_off[!,s], y=eigs_imag_off[!,s], mode="markers"), Layout(xaxis_title="Real", yaxis_title="Imaginary"))

plot([scatter(x=eigs_real_on[!,s], y=eigs_imag_on[!,s], mode="markers"), scatter(x=eigs_real_off[!,s], y=eigs_imag_off[!,s], mode="markers")], Layout(xaxis_title="Real", yaxis_title="Imaginary"))

s = :rm_r
plot([scatter(x=kdam_range, y=eigs_real_on[!,s], mode="markers"), scatter(x=kdam_range, y=eigs_real_off[!,s], mode="markers")], Layout(xaxis_title="kdam", yaxis_title="Eigenvalue"))
