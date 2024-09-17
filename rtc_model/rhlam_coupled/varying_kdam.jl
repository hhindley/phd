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
