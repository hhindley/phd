include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

res=[]
minus_num = 20
for i in eachindex(ssvals_rtc)
    ssvals_diff = deepcopy(ssvals_rtc)
    if ssvals_diff[i]-minus_num < 0
        ssvals_diff[i] = 0
    else 
        ssvals_diff[i] = ssvals_diff[i]-minus_num
    end
    println(ssvals_diff)
    df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
    push!(res, df_diff)
end
plot([scatter(x=kdam_range, y=res[i].rtca, name="$(species_rtc[i]) - $minus_num") for i in eachindex(res)])

res1=[]
minus_num = 20
for i in eachindex(ssvals_rtc)
    ssvals_diff = deepcopy(ssvals_rtc)
    if ssvals_diff[i]-minus_num < 0 
        ssvals_diff[i] = 0
    else 
        ssvals_diff[i] = ssvals_diff[i]-minus_num
    end
    if i+1 < length(ssvals_diff)
        if ssvals_diff[i+1]-minus_num < 0 
            ssvals_diff[i+1] = 0
        else 
            ssvals_diff[i+1] = ssvals_diff[i+1]-minus_num
        end
        if ssvals_diff[i+2]-minus_num < 0 
            ssvals_diff[i+2] = 0
        else 
            ssvals_diff[i+2] = ssvals_diff[i+2]-minus_num
        end
        println(ssvals_diff)
        df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
        push!(res1, df_diff)
    end
end
res1
plot([scatter(x=kdam_range, y=res1[i], name="$(species_rtc[i]) - $minus_num") for i in eachindex(res1)])



res1 = []
combs = []
minus_num = 2
num_conditions = length(ssvals_rtc)
for k in eachindex(ssvals_rtc)
    println(k)
    for comb in combinations(1:num_conditions, k)
        # println(comb)
        ssvals_diff = deepcopy(ssvals_rtc)
        for i in comb
            # println(i)
            if ssvals_diff[i] - minus_num < 0
                ssvals_diff[i] = 0
            else
                ssvals_diff[i] = ssvals_diff[i] - minus_num
            end
        end

        df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)

        if df_diff.rtca[end] != round(df_ssvals.rtca[end], digits=6)
            # println(ssvals_diff)
            push!(res1, df_diff)
            push!(combs, comb)
        end
        
    end
end

p = plot([scatter(x=kdam_range, y=res1[i].rtca, name="change $([species_rtc[s] for s in combs[i]])") for i in eachindex(res1)], Layout(xaxis_title="kdam",yaxis_title="rtca",title="$([round(ssvals_rtc[i], digits=6) for i in eachindex(ssvals_rtc)]) - $minus_num"))
open("/Users/s2257179/Desktop/minus2.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
