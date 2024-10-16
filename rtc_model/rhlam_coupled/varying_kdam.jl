include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]


df_ssvals = var_param(test, kdam, params_rtc1, kdam_range, ssvals_rtc)
plot(scatter(x=kdam_range, y=df_ssvals.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

df_init = var_param(test, kdam, params_rtc1, kdam_range, init_rtc)
plot(scatter(x=kdam_range, y=df_init.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

df_on = var_param(test, kdam, params_rtc1, kdam_range, init_on)
plot(scatter(x=kdam_range, y=df_on.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

df_off = var_param(test, kdam, params_rtc1, kdam_range, init_off)
plot(scatter(x=kdam_range, y=df_off.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

ssvals_diff = deepcopy(ssvals_rtc)
ssvals_diff[2] = 0
df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
plot(scatter(x=kdam_range, y=df_diff.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))


L_range = [10, 100, 1000, 10000, 100000]
c_range = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

L_ssvals = []
for i in L_range
    L_params = deepcopy(params_rtc1)
    L_params[L] = i
    push!(L_ssvals, var_param(test, kdam, L_params, kdam_range, ssvals_rtc))
end
plot([scatter(x=kdam_range, y=i.rtca, name="L = $(L_range[j])", line=attr(color=colours[j])) for (i,j) in zip(L_ssvals, 1:length(L_range))], Layout(xaxis_title="kdam", yaxis_title="rtca"))

c_ssvals = []
for i in c_range
    c_params = deepcopy(params_rtc1)
    c_params[c] = i
    push!(c_ssvals, var_param(test, kdam, c_params, kdam_range, ssvals_rtc))
end
plot([scatter(x=kdam_range, y=i.rtca, name="c = $(c_range[j])", line=attr(color=colours[j])) for (i,j) in zip(c_ssvals, 1:length(c_range))], Layout(xaxis_title="kdam", yaxis_title="rtca"))

Lc_ssvals = []
Lc_vals = []
for i in L_range
    for j in c_range
        Lc_params = deepcopy(params_rtc1)
        Lc_params[L] = i
        Lc_params[c] = j
        push!(Lc_vals, (i,j))
        push!(Lc_ssvals, var_param(test, kdam, Lc_params, kdam_range, ssvals_rtc))
    end
end

plot([scatter(x=kdam_range, y=Lc_ssvals[i].rtca, name="L,c = $(Lc_vals[i])") for i in eachindex(Lc_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))




lam_c_val = range(1e-8,1, length=20) #8e-7
kin_c_val = range(1e-8,1, length=20) #1.5e-5
kdam_range = range(0,100,length=20)
kdam_range1 = reverse(kdam_range)
lamkin_ssvals = []
lamkin_ssvals_vals = []
inc = []
dec = []
for i in lam_c_val
    for j in kin_c_val
        lamkin_ssvals_params = deepcopy(params_rtc1)
        lamkin_ssvals_params[lam_c] = i
        lamkin_ssvals_params[kin_c] = j
        push!(lamkin_ssvals_vals, (i,j))
        # push!(lamkin_ssvals, var_param(test, kdam, lamkin_ssvals_params, kdam_range, ssvals_rtc))
        res_ss = numerical_bistability_analysis(test, lamkin_ssvals_params, init_rtc, species_rtc, kdam_range, kdam)
        res1_ss = numerical_bistability_analysis(test, lamkin_ssvals_params, ssvals_rtc, species_rtc, kdam_range1, kdam)
        push!(inc, res_ss)
        push!(dec, res1_ss)

    end
end
plot([[scatter(x=kdam_range, y=inc[i].rtca, name="lam_c, kin_c = $(lamkin_ssvals_vals[i])") for i in eachindex(lamkin_ssvals_vals)],
[scatter(x=kdam_range, y=inc[i].rtca, name="lam_c, kin_c = $(lamkin_ssvals_vals[i])") for i in eachindex(lamkin_ssvals_vals)]],
 Layout(xaxis_title="kdam", yaxis_title="rtca"))


plot([scatter(x=kdam_range, y=lamkin_ssvals[i].rtca, name="lam_c, kin_c = $(lamkin_ssvals_vals[i])") for i in eachindex(lamkin_ssvals_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))
