species_rtc = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]

# sort dataframes
# function df_sort(df)
#     df.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df.event]
#     df.event = [parse.(Float64, subarray) for subarray in df.event]
#     df_p = filter(row -> length(row[:event]) > 1, df)
#     df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)
#     df_r.event = first.(df_r.event)

#     props = df_p.event
#     react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

#     df_props = DataFrame([name => Float64[] for name in react_names])

#     for i in props
#         push!(df_props, [i[j] for j in eachindex(props[end])])
#     end
#     return df_props, df_r, df_p
# end

# sort multiple dataframes
# function df_sort_all(dfs)
#     props=[]; df_rs=[]; df_ps=[];
#     for i in ProgressBar(eachindex(dfs))
#         df_props, df_r, df_p = df_sort(dfs[i])
#         push!(props, df_props)
#         push!(df_rs, df_r)
#         push!(df_ps, df_p)
#     end
#     return props, df_rs, df_ps
# end

# look at stochastic reactions
# function create_react_df(df)
#     react_nums=[]
#     for i in 1:13
#         push!(react_nums, length(filter(row->row[:event] == i, df).event))
#     end

#     df_react = DataFrame(reaction=react_names, count=react_nums)
#     return df_react
# end

# function all_react_dfs(dfs, df_rs)
#     df_reacts=[]; tot_stoch=[];
#     for i in eachindex(dfs)
#         df_react = create_react_df(df_rs[i])
#         push!(df_reacts, df_react)
#         push!(tot_stoch, sum(df_react.count))
#     end
#     return df_reacts, tot_stoch
# end

# histograms
function create_hist(df_ps, x, n)
    df_ps[x].bins = ceil.(Int, (1:nrow(df_ps[x]))/n)
    df_grouped = combine(first, groupby(df_ps[x], :bins))
    nbins=length(df_grouped.bins)
    return df_grouped, nbins
end

function all_hists(df_ps, kdam_vals, n)
    groups=[]; bins=[];
    for i in eachindex(kdam_vals)
        df_grouped, nbins = create_hist(df_ps, i, n)
        push!(groups, df_grouped)
        push!(bins, nbins)
    end
    return groups, bins
end

# % expression
function calc_percentage_exp(df_ps, x, specie)
    exp = filter(x->x>1, df_ps[x][:,specie])
    return (length(exp)/length(df_ps[x][:,specie]))*100
end

function all_exp(df_ps, kdam_vals)
    exp_df = DataFrame(kdam=[], rm_a=[], rm_b=[], rm_r=[], rtca=[], rtcb=[], rtcr=[]);#, rh=[], rd=[], rt=[]);
    for x in eachindex(kdam_vals)
        push!(exp_df, [kdam_vals[x], calc_percentage_exp(df_ps, x, :rm_a), calc_percentage_exp(df_ps, x, :rm_b), calc_percentage_exp(df_ps, x, :rm_r), calc_percentage_exp(df_ps, x, :rtca), calc_percentage_exp(df_ps, x, :rtcb), calc_percentage_exp(df_ps, x, :rtcr)])#, calc_percentage_exp(x, :rh), calc_percentage_exp(x, :rd), calc_percentage_exp(x, :rt)])
    end
    return exp_df
end

# means vals
function mean_vals(df_ps, kdam_vals, n)
    df_av = DataFrame(kdam=[], rm_a=[], rm_b=[], rm_r=[], rtca=[], rtcb=[], rtcr=[], rh=[], rd=[], rt=[]);
    for x in eachindex(kdam_vals)
        ss = findfirst(x->x>n, df_ps[x].time)
        means = mean.(eachcol(df_ps[x][ss:end,3:11]))
        pushfirst!(means, kdam_vals[x])
        push!(df_av, means)
    end
    return df_av
end