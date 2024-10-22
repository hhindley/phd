function remove_missing(df_rtca, df_times)
    res = Dict()
    times_res = Dict()
    for i in eachindex(names(df_rtca))
        res[parse(Float64, names(df_rtca)[i])] = collect(skipmissing(df_rtca[!,names(df_rtca)[i]]))
        times_res[parse(Float64, names(df_rtca)[i])] = collect(skipmissing(df_times[!,names(df_times)[i]]))
    end
    return res, times_res 
end

function log_results(res)
    res_log = Dict()
    for kdam in kdams
        res_log[kdam] = log.(res[kdam].+1)
    end
    return res_log 
end

function determine_state(res, times; threshold=2)
    res_on = Dict([kdam=>Float64[] for kdam in kdams])
    res_off = Dict([kdam=>Float64[] for kdam in kdams])
    res_on_times = Dict([kdam=>Float64[] for kdam in kdams])
    res_off_times = Dict([kdam=>Float64[] for kdam in kdams])
    for kdam in kdams
        for i in eachindex(res[kdam])
            if res[kdam][i] >= threshold
                push!(res_on[kdam], res[kdam][i])
                push!(res_on_times[kdam], times[kdam][i])
            else
                push!(res_off[kdam], res[kdam][i])
                push!(res_off_times[kdam], times[kdam][i])
            end
        end
    end
    return res_on, res_off, res_on_times, res_off_times
end

function all_results_concat(res, df_lengths)
    all_res=[]
    for kdam in kdams
        push!(all_res, res[kdam])
    end
    tot_lengths = [sum(df_lengths[!,"$kdam"]) for kdam in kdams]
    groups = [fill(i, tot_lengths[i]) for i in eachindex(kdams)]
    kdam = [fill(kdams[i], tot_lengths[i]) for i in eachindex(kdams)]
    df_res = DataFrame(rtca=vcat(all_res...), kdam=vcat(kdam...), group=vcat(groups...))
    return df_res
end

function split_into_simulations(res, times_res, df_stops)
    sims = Dict(i=>Dict(kdam=>Float64[] for kdam in kdams) for i in 1:20)
    times = Dict(i=>Dict(kdam=>Float64[] for kdam in kdams) for i in 1:20)
    for kdam in eachindex(kdams)
        sims[1][kdams[kdam]] = res[kdams[kdam]][1:df_stops[1,"$(kdams[kdam])"]]
        times[1][kdams[kdam]] = times_res[kdams[kdam]][1:df_stops[1,"$(kdams[kdam])"]]
        for j in eachindex(kdams)
            if j < 20
                ind = j+1
            else 
                continue
            end
            sims[ind][kdams[kdam]] = res[kdams[kdam]][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
            times[ind][kdams[kdam]] = times_res[kdams[kdam]][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
        end
    end
    return sims, times
end




