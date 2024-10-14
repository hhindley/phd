
function find_closest_index(array, value)
    differences = abs.(array .- value)
    closest_index = argmin(differences)
    return closest_index
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

    min_kdam_ind = find_closest_index(dict_kdamvals[folder][:kdam], min_bs)
    max_kdam_ind = find_closest_index(dict_kdamvals[folder][:kdam], max_bs)

    thresholds_rtca=[]
    thresholds_rtcb=[]
    for i in dict_kdamvals[folder][:kdam][min_kdam_ind:max_kdam_ind]
        println(i)
        index = find_closest_index(unstab_kdam, i)
        push!(thresholds_rtca, unstab_rtca[index])
        push!(thresholds_rtcb, unstab_rtcb[index])
    end

    start_vals_rtca = fill(thresholds_rtca[1], length(dict_kdamvals[folder][:kdam][1:min_kdam_ind-1]))
    # either fill the values after the unstable curve with the first threshold or the last (but the last is a lot higher than the first, and this means that sometimes RtcA or RtcB don't cross this threshold, so will go with the first)
    end_vals_rtca = fill(thresholds_rtca[1], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))
    # end_vals_rtca = fill(thresholds_rtca[end], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))

    start_vals_rtcb = fill(thresholds_rtcb[1], length(dict_kdamvals[folder][:kdam][1:min_kdam_ind-1]))
    end_vals_rtcb = fill(thresholds_rtcb[1], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))
    # end_vals_rtcb = fill(thresholds_rtcb[end], length(dict_kdamvals[folder][:kdam][max_kdam_ind+1:end]))

    pushfirst!(thresholds_rtca, start_vals_rtca...)
    push!(thresholds_rtca, end_vals_rtca...)

    pushfirst!(thresholds_rtcb, start_vals_rtcb...)
    push!(thresholds_rtcb, end_vals_rtcb...)

    return thresholds_rtca, thresholds_rtcb
end

function set_condition(df; threshold_rtca::Union{Int, Float64}=0, threshold_rtcb::Union{Int, Float64}=0) # just for one dataframe (one kdam val)
    df_rtca = df.rtca
    df_rtcb = df.rtcb

    condition = (df_rtca .> threshold_rtca) .| (df_rtcb .> threshold_rtcb)

    return condition
end

function start_stop_indices(df, condition)
    df_time = df.time

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

    return start_indices, stop_indices
end

function switch_times(df, start_indices, stop_indices)
    df_time = getproperty(df, :time)

    start_times = []
    for i in start_indices
        push!(start_times, df_time[i])
    end
    stop_times = []
    for i in stop_indices
        push!(stop_times, df_time[i])
    end

    return start_times, stop_times
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

function calc_frac_times(switch_vals_on, switch_vals_off)
    frac_on = @. switch_vals_off/(switch_vals_on+switch_vals_off)
    frac_off = @. switch_vals_on/(switch_vals_on+switch_vals_off)
    return frac_on, frac_off
end

function calc_av_state_conc(df, species, start_indices, stop_indices; on=true)
    df_volume = getproperty(df, :volume)
    if typeof(species) == Symbol
        species = [species]
    end
    all_species = [getproperty(df,i) for i in species]

    if on
        # calculate mean in batches 
        # species_vals = [[@view all_species[i][start:stop] for (start, stop) in zip(start_indices, stop_indices)] for i in eachindex(all_species)]
        # vol_vals = [@view df_volume[start:stop] for (start, stop) in zip(start_indices, stop_indices)]
        # species_mean = [mean([mean(sf.*(species ./ volume)) for (species, volume) in zip(species_vals[i], vol_vals)]) for i in eachindex(species)]
        
        # calculate mean over whole region
        species_vals_batch = [[@view all_species[i][start:stop] for (start, stop) in zip(start_indices, stop_indices)] for i in eachindex(all_species)]
        vol_vals_batch = [@view df_volume[start:stop] for (start, stop) in zip(start_indices, stop_indices)]

        species_vals = [vcat(species_vals_batch[i]...) for i in eachindex(all_species)]
        vol_vals = vcat(vol_vals_batch...)

        concs = [sf*(species_vals[i]./vol_vals) for i in eachindex(all_species)]

        species_mean = [mean(concs[i]) for i in eachindex(all_species)]

    else
        if length(start_indices) == 1 && length(stop_indices) == 1
            species_vals_batch1 = [@view all_species[i][1:start_indices[1]] for i in eachindex(all_species)]
            vol_vals_batch1 = @view df_volume[1:start_indices[1]]
            species_vals_batch2 = [@view all_species[i][stop_indices[1]:end] for i in eachindex(all_species)]
            vol_vals_batch2 = @view df_volume[stop_indices[1]:end]

            species_vals = [vcat(species_vals_batch1[i], species_vals_batch2[i]) for i in eachindex(all_species)]
            vol_vals = vcat(vol_vals_batch1, vol_vals_batch2)

            concs = [sf*(species_vals[i]./vol_vals) for i in eachindex(all_species)]

            species_mean = [mean(concs[i]) for i in eachindex(all_species)]
    
        else
            new_starts = start_indices[2:end]
            new_stops = stop_indices[1:end-1]

            species_vals_batch = [[@view all_species[i][stop:start] for (start, stop) in zip(new_starts, new_stops)] for i in eachindex(all_species)]
            vol_vals_batch = [@view df_volume[stop:start] for (start, stop) in zip(new_starts, new_stops)]

            species_vals = [vcat(species_vals_batch[i]...) for i in eachindex(all_species)]
            vol_vals = vcat(vol_vals_batch...)

            concs = [sf*(species_vals[i]./vol_vals) for i in eachindex(all_species)]

            species_mean = [mean(concs[i]) for i in eachindex(all_species)]
            
            # batch
            # species_mean = [mean([mean(sf.*(species ./ volume)) for (species, volume) in zip(species_vals[i], vol_vals)]) for i in eachindex(species)]
        end
    end

    return species_mean
end

# calculations for a folder at a time 
function folder_indices(folder, threshold_rtca::Union{Int, Vector{Any}}; threshold_rtcb::Union{Vector{Any}, Nothing}=nothing)
    kdamvals = dict_kdamvals[folder][:kdam]
    results = dict_results[folder]

    start_indices_f = Dict(key => Vector{Int}() for key in kdamvals)
    stop_indices_f = Dict(key => Vector{Int}() for key in kdamvals)
    
    if typeof(threshold_rtca) == Int
        threshold_rtca = fill(threshold_rtca, length(kdamvals))
    end

    threshold_rtcb = threshold_rtcb === nothing ? threshold_rtca : threshold_rtcb

    for (kdam, ind) in zip(kdamvals, eachindex(kdamvals))
        condition = set_condition(results[ind], threshold_rtca=threshold_rtca[ind], threshold_rtcb=threshold_rtcb[ind])
        start_indices_f[kdam], stop_indices_f[kdam] = start_stop_indices(results[ind], condition)
    end
    return start_indices_f, stop_indices_f
end

function folder_switchrates_fracs(folder, start_indices_f, stop_indices_f)
    kdamvals = dict_kdamvals[folder][:kdam]
    results = dict_results[folder]

    switch_rates_on_f = Dict(key => 0.0 for key in kdamvals)
    switch_rates_off_f = Dict(key => 0.0 for key in kdamvals)
    fracs_on_f = Dict(key => 0.0 for key in kdamvals)
    fracs_off_f = Dict(key => 0.0 for key in kdamvals)

    for (kdam, ind) in zip(kdamvals, eachindex(kdamvals))
        start_times, stop_times = switch_times(results[ind], start_indices_f[kdam], stop_indices_f[kdam])
        switch_rate_on, switch_rate_off = calc_switch_rate(results[ind], start_times, stop_times)
        switch_rates_on_f[kdam] = switch_rate_on
        switch_rates_off_f[kdam] = switch_rate_off

        frac_on, frac_off = calc_frac_times(switch_rate_on, switch_rate_off)
        fracs_on_f[kdam] = frac_on
        fracs_off_f[kdam] = frac_off
    end

    return switch_rates_on_f, switch_rates_off_f, fracs_on_f, fracs_off_f
end

function folder_concs(folder, species, start_indices_f, stop_indices_f)
    species_mean_on_f = Dict(key => Vector{Float64}() for key in dict_kdamvals[folder][:kdam])
    species_mean_off_f = Dict(key => Vector{Float64}() for key in dict_kdamvals[folder][:kdam])

    results = dict_results[folder]

    for (kdam, ind) in zip(dict_kdamvals[folder][:kdam], eachindex(dict_kdamvals[folder][:kdam]))
        species_mean_on_f[kdam] = calc_av_state_conc(results[ind], species, start_indices_f[kdam], stop_indices_f[kdam])
        species_mean_off_f[kdam] = calc_av_state_conc(results[ind], species, start_indices_f[kdam], stop_indices_f[kdam], on=false)
    end

    return species_mean_on_f, species_mean_off_f
end


# functions for doing all of the folders loaded in the current folders_dict
function all_indices(folders_dict, threshold_rtca::Union{Int, Vector{Any}}; threshold_rtcb::Union{Vector{Any}, Nothing}=nothing)
    all_start_indices = Dict(key => Dict() for key in eachindex(folders_dict))
    all_stop_indices = Dict(key => Dict() for key in eachindex(folders_dict))

    for folder in eachindex(folders_dict)
        start_indices, stop_indices = folder_indices(folder, threshold_rtca, threshold_rtcb=threshold_rtcb)
        all_start_indices[folder] = start_indices
        all_stop_indices[folder] = stop_indices
    end

    return all_start_indices, all_stop_indices
end

function all_switchrates_fracs(folders_dict, all_start_indices, all_stop_indices)
    all_switch_rates_on = Dict(key => Dict() for key in eachindex(folders_dict))
    all_switch_rates_off = Dict(key => Dict() for key in eachindex(folders_dict))
    all_fracs_on = Dict(key => Dict() for key in eachindex(folders_dict))
    all_fracs_off = Dict(key => Dict() for key in eachindex(folders_dict))

    for folder in eachindex(folders_dict)
        switch_rates_on_f, switch_rates_off_f, fracs_on_f, fracs_off_f = folder_switchrates_fracs(folder, all_start_indices[folder], all_stop_indices[folder])
        all_switch_rates_on[folder] = switch_rates_on_f
        all_switch_rates_off[folder] = switch_rates_off_f
        all_fracs_on[folder] = fracs_on_f
        all_fracs_off[folder] = fracs_off_f
    end

    return all_switch_rates_on, all_switch_rates_off, all_fracs_on, all_fracs_off
end

function all_concs(folders_dict, species, all_start_indices, all_stop_indices)
    all_species_mean_on = Dict(key => Dict() for key in eachindex(folders_dict))
    all_species_mean_off = Dict(key => Dict() for key in eachindex(folders_dict))

    for folder in eachindex(folders_dict)
        species_mean_on_f, species_mean_off_f = folder_concs(folder, species, all_start_indices[folder], all_stop_indices[folder])
        all_species_mean_on[folder] = species_mean_on_f
        all_species_mean_off[folder] = species_mean_off_f
    end

    return all_species_mean_on, all_species_mean_off
end

function calc_mean_std_vars(switch_rates, fracs, kdams)
    mean_switch_frac = Dict("on"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())), 
                        "off"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())))
    std_switch_frac = Dict("on"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())), 
                        "off"=>Dict("2"=>Dict("switch"=>Dict(), "frac"=>Dict()), "5"=>Dict("switch"=>Dict(), "frac"=>Dict()), "10"=>Dict("switch"=>Dict(), "frac"=>Dict()), "bs"=>Dict("switch"=>Dict(), "frac"=>Dict())))

    for onoff in ["on", "off"]
        for thresh in ["2", "5", "10", "bs"]
            for kdam in kdams
                rates_switch = [switch_rates[onoff][thresh][i][kdam] for i in eachindex(switch_rates[onoff][thresh])]
                mean_switch_frac[onoff][thresh]["switch"][kdam] = mean(rates_switch)
                std_switch_frac[onoff][thresh]["switch"][kdam] = std(rates_switch)
                
                fracs1 = [fracs[onoff][thresh][i][kdam] for i in eachindex(fracs[onoff][thresh])]
                mean_switch_frac[onoff][thresh]["frac"][kdam] = mean(fracs1)
                std_switch_frac[onoff][thresh]["frac"][kdam] = std(fracs1)
            end
        end
    end
    return mean_switch_frac, std_switch_frac
end