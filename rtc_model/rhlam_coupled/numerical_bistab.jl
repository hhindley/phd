include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

kdam_range = range(0.0, 5, length=100)
kdam_range1 = reverse(kdam_range)
res = numerical_bistability_analysis(test, params_rtc1, init_rtc, species_rtc, kdam_range, kdam)
res1 = numerical_bistability_analysis(test, params_rtc1, init_rtc, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res.rtca), scatter(x=kdam_range1, y=res1.rtca)])

plot([scatter(x=kdam_range1, y=res1[!,s], name="$s") for s in species_rtc])


res_ss = numerical_bistability_analysis(test, params_rtc1, init_rtc, species_rtc, kdam_range, kdam)
res1_ss = numerical_bistability_analysis(test, params_rtc1, ssvals_rtc, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res_ss.rtca), scatter(x=kdam_range1, y=res1_ss.rtca)])



init_on = [res[1,s] for s in species_rtc]
init_off = [res1[1,s] for s in species_rtc]

# check stability of branches 

param_dict_dam = deepcopy(params_rtc1)

unstable, eigenvals_on = calc_unstable_points(kdam_range, res, param_dict_dam)
unstable1, eigenvals_off = calc_unstable_points(kdam_range1, res1, param_dict_dam)


p = scatter(x=kdam_range, y=res.rtca)
p1 = scatter(x=kdam_range1, y=res1.rtca)

pu = [scatter(x=[kdam_range1[i]], y=[res1[i,:rtca]], mode="markers") for i in unstable1]

plot([p, p1, pu...], Layout(xaxis_title="kdam", yaxis_title="rtca"))

kdam_full = [fill(i, 9) for i in kdam_range]
kdam_full = vcat(kdam_full...)

real_eigs_on = [real.(e) for e in eigenvals_on]
real_eigs_on = vcat(real_eigs_on...)
imag_eigs_on = [imag.(e) for e in eigenvals_on]
imag_eigs_on = vcat(imag_eigs_on...)
plot(scatter(x=kdam_full,y=real_eigs_on,mode="markers"), Layout(xaxis_title="kdam",yaxis_title="Real part of eigenvalues", title="on branch"))
plot(scatter(x=real_eigs_on,y=imag_eigs_on,mode="markers"), Layout(xaxis_title="Real part of eigenvalues",yaxis_title="Imaginary part of eigenvalues", title="on branch"))

real_eigs_off = [real.(e) for e in eigenvals_off]
real_eigs_off = vcat(real_eigs_off...)
imag_eigs_off = [imag.(e) for e in eigenvals_off]
imag_eigs_off = vcat(imag_eigs_off...)
plot(scatter(x=kdam_full,y=real_eigs_off,mode="markers"), Layout(xaxis_title="kdam",yaxis_title="Real part of eigenvalues", title="off branch"))
plot(scatter(x=real_eigs_off,y=imag_eigs_off,mode="markers"), Layout(xaxis_title="Real part of eigenvalues",yaxis_title="Imaginary part of eigenvalues", title="off branch"))



typeof(eigenvals_on[2])

complex_on = [e for e in eachindex(eigenvals_on) if typeof(eigenvals_on[e]) == Vector{ComplexF64}]

complex_off = [e for e in eachindex(eigenvals_off) if typeof(eigenvals_off[e]) == Vector{ComplexF64}]


p = scatter(x=kdam_range, y=res.rtca)
p1 = scatter(x=kdam_range1, y=res1.rtca)

pu = [scatter(x=[kdam_range[i]], y=[res[i,:rtca]], mode="markers") for i in complex_on]
pu1 = [scatter(x=[kdam_range1[i]], y=[res1[i,:rtca]], mode="markers") for i in complex_off]

plot([p; p1; pu; pu1], Layout(xaxis_title="kdam", yaxis_title="rtca"))
