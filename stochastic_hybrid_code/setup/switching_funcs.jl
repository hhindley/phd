# two versions of determining the switching rate
# if rtca and rtcb are both above a certain threshold
# or if rtca and rtcb are both increasing then we are in the on state, if they are not changing or are decreasing then we are in the off state 
function switch_times(res; threshold=nothing)
    df = res[!,[:time, :rtca, :rtcb]]
    if threshold !== nothing
        condition = (df.rtca .> threshold) .| (df.rtcb .> threshold)
    else
        diff_rtca = diff(df.rtca)
        diff_rtcb = diff(df.rtcb)
        condition = (diff_rtca .> 1e-10) .| (diff_rtcb .> 1e-10)
        pushfirst!(condition, !condition[1])
    end
    df = insertcols!(df, 1, :on => condition)
    condition_changes = diff(Int.(condition))
    start_indices = findall(x -> x == 1, condition_changes)
    stop_indices = findall(x -> x == -1, condition_changes)
    start_times = []
    for i in start_indices
        push!(start_times, df.time[i])
    end
    stop_times = []
    for i in stop_indices
        push!(stop_times, df.time[i])
    end
    if start_times[1] > stop_times[1] && start_times[end] > stop_times[end]
        # println("start_times[1] > stop_times[1] && start_times[end] > stop_times[end]")
        pushfirst!(start_times, df.time[1])
        push!(stop_times, df.time[end])
    elseif start_times[end] > stop_times[end]
        # println("start_times[end] > stop_times[end]")
        push!(stop_times, df.time[end])
    elseif start_times[1] > stop_times[1] 
        # println("start_times[1] > stop_times[1] ")
        pushfirst!(start_times, df.time[1])
    end
    return df, start_times, stop_times
end

function calc_switch_rate(df, start_times, stop_times)
    time_on = sum(stop_times.-start_times)
    on_res = time_on/length(start_times)
    time_off = df.time[end] - time_on
    off_res = time_off/length(stop_times)
    switch_rate_on = 1/on_res
    switch_rate_off = 1/off_res
    return switch_rate_on, switch_rate_off
end

function get_all_switch_rates(folder; threshold=nothing)
    switch_rates_on = []
    switch_rates_off = []
    for i in eachindex(dict_results[folder])
        # println(i)
        df, start_times, stop_times = switch_times(dict_results[folder][i], threshold=threshold)
        switch_rate_on, switch_rate_off = calc_switch_rate(df, start_times, stop_times)
        push!(switch_rates_on, switch_rate_on)
        push!(switch_rates_off, switch_rate_off)
    end
    return switch_rates_on, switch_rates_off
end

function calc_frac_times(switch_vals_on, switch_vals_off)
    frac_on = @. switch_vals_off/(switch_vals_on+switch_vals_off)
    frac_off = @. switch_vals_on/(switch_vals_on+switch_vals_off)
    return frac_on, frac_off
end

function get_av_conc_state(res, threshold, species; on=true, dynamic=false)
    df = res
    if typeof(species) == Symbol
        species = [species]
    end
    df_time = getproperty(df, :time)
    df_volume = getproperty(df, :volume)
    df_rtca = getproperty(df, :rtca)
    df_rtcb = getproperty(df, :rtcb)
    
    if species == [:rtca, :rtcb]
        all_species = [df_rtca, df_rtcb]
    else
        all_species = [getproperty(df,i) for i in species]
    end

    if dynamic == false
        condition = (df_rtca .> threshold) .| (df_rtcb .> threshold)
    else 
        condition = (df_rtca .> threshold_rtca) .& (df_rtcb .> threshold_rtcb)
    end

    condition_changes = diff(Int.(condition))

    condition_changes = diff(Int.(condition))
    start_indices = findall(x -> x == 1, condition_changes)
    stop_indices = findall(x -> x == -1, condition_changes)

    if start_indices[1] > stop_indices[1] && start_indices[end] > stop_indices[end]
        pushfirst!(start_indices, 1)
        push!(stop_indices, length(df_time))
    elseif start_indices[end] > stop_indices[end]
        push!(stop_indices, length(df_time))
    elseif start_indices[1] > stop_indices[1] 
        pushfirst!(start_indices, 1)
    end

    if on
        species_vals = [[@view all_species[i][start:stop] for (start, stop) in zip(start_indices, stop_indices)] for i in eachindex(all_species)]
        vol_vals = [@view df_volume[start:stop] for (start, stop) in zip(start_indices, stop_indices)]

        species_mean = [mean([mean(sf.*(species ./ volume)) for (species, volume) in zip(species_vals[i], vol_vals)]) for i in eachindex(species)]
    
    else
        new_starts = start_indices[2:end]
        new_stops = stop_indices[1:end-1]

        species_vals = [[@view all_species[i][stop:start] for (start, stop) in zip(new_starts, new_stops)] for i in eachindex(all_species)]
        vol_vals = [@view df_volume[stop:start] for (start, stop) in zip(new_starts, new_stops)]

        species_mean = [mean([mean(sf.*(species ./ volume)) for (species, volume) in zip(species_vals[i], vol_vals)]) for i in eachindex(species)]

    end

    return species_mean

end

function get_all_av_conc(folder, threshold, species)
    species_on = []
    species_off = []
    for i in eachindex(dict_results[folder])
        # println(i)
        if typeof(threshold) == Array
            for i in threshold
                means_on = get_av_conc_state(dict_results[folder][i], i, species, on=true)
                means_off = get_av_conc_state(dict_results[folder][i], i, species, on=false)
                push!(species_on, means_on)
                push!(species_off, means_off)
            end
        else
            means_on = get_av_conc_state(dict_results[folder][i], threshold, species, on=true)
            means_off = get_av_conc_state(dict_results[folder][i], threshold, species, on=false)
            push!(species_on, means_on)
            push!(species_off, means_off)
        end
    end
    return species_on, species_off
end

function get_unstab_threshold_array(folder)
    br = get_br_molec(rtc_model, ssvals_rtc_molec, params_rtc_molec, 1.5)
    bf = bf_point_df(br)
    df = create_br_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    unstab_rtca = df[!,:rtca][kdam1:kdam2]
    unstab_rtcb = df[!,:rtcb][kdam1:kdam2]
    unstab_kdam = df[!,:kdam][kdam1:kdam2]

    min_bs = ceil(minimum(unstab_kdam), digits=1)
    max_bs = floor(maximum(unstab_kdam), digits=1)

    min_kdam_ind = find_closest_index(dict_kdamvals[6][:kdam], min_bs)
    max_kdam_ind = find_closest_index(dict_kdamvals[6][:kdam], max_bs)

    thresholds_rtca=[]
    thresholds_rtcb=[]
    for i in dict_kdamvals[folder][:kdam][min_kdam_ind:max_kdam_ind]
        println(i)
        index = find_closest_index(unstab_kdam, i)
        push!(thresholds_rtca, unstab_rtca[index])
        push!(thresholds_rtcb, unstab_rtcb[index])
    end

    start_vals_rtca = fill(thresholds_rtca[1], length(dict_kdamvals[folder][:kdam][1:min_kdam_ind-1]))
    end_vals_rtca = fill(thresholds_rtca[end], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))

    start_vals_rtcb = fill(thresholds_rtcb[1], length(dict_kdamvals[folder][:kdam][1:min_kdam_ind-1]))
    end_vals_rtcb = fill(thresholds_rtcb[end], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))

    pushfirst!(thresholds_rtca, start_vals_rtca...)
    push!(thresholds_rtca, end_vals_rtca...)

    pushfirst!(thresholds_rtcb, start_vals_rtcb...)
    push!(thresholds_rtcb, end_vals_rtcb...)

    return thresholds_rtca, thresholds_rtcb
end