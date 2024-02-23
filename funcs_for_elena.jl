function setup_ssvals_from_bfkit(rtc_mod, kdam_val, params2, initial, kdam)

    br2 = get_br(rtc_mod, initial, params2, kdam)
    df = create_br_df(br2)
    df_bf = bf_point_df(br2)
    first,middle,last=split_curves(df, df_bf)
    first=Float64.(first)
    last=Float64.(last)
    
    int_rma1 = QuadraticInterpolation(first.rm_a, first.kdam)
    int_rma2 = QuadraticInterpolation(last.rm_a, last.kdam)
    int_rtca1 = QuadraticInterpolation(first.rtca, first.kdam)
    int_rtca2 = QuadraticInterpolation(last.rtca, last.kdam)
    int_rmb1 = QuadraticInterpolation(first.rm_b, first.kdam)
    int_rmb2 = QuadraticInterpolation(last.rm_b, last.kdam)
    int_rtcb1 = QuadraticInterpolation(first.rtcb, first.kdam)
    int_rtcb2 = QuadraticInterpolation(last.rtcb, last.kdam)
    int_rmr1 = QuadraticInterpolation(first.rm_r, first.kdam)
    int_rmr2 = QuadraticInterpolation(last.rm_r, last.kdam)
    int_rtcr1 = QuadraticInterpolation(first.rtcr, first.kdam)
    int_rtcr2 = QuadraticInterpolation(last.rtcr, last.kdam)
    int_rh1 = QuadraticInterpolation(first.rh, first.kdam)
    int_rh2 = QuadraticInterpolation(last.rh, last.kdam)
    int_rd1 = QuadraticInterpolation(first.rd, first.kdam)
    int_rd2 = QuadraticInterpolation(last.rd, last.kdam)
    int_rt1 = QuadraticInterpolation(first.rt, first.kdam)
    int_rt2 = QuadraticInterpolation(last.rt, last.kdam)

    return DataFrame(species=species_rtc,ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val)],ss_val_off=[int_rma2(kdam_val),int_rtca2(kdam_val),int_rmb2(kdam_val),int_rtcb2(kdam_val),int_rmr2(kdam_val),int_rtcr2(kdam_val),int_rh2(kdam_val),int_rd2(kdam_val),int_rt2(kdam_val)]);
end

function get_all_ranges(func, branch_df, branch, n, l)
    all_ranges=[]
    for i in species_rtc#[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i]
        push!(all_ranges, func(branch_df, branch, i, n, l))
    end
    return @LArray [all_ranges[1], all_ranges[2],all_ranges[3],all_ranges[4],all_ranges[5],all_ranges[6],all_ranges[7],all_ranges[8],all_ranges[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt)  #[all_ranges[1], all_ranges[2],all_ranges[3],all_ranges[4],all_ranges[5],all_ranges[6],all_ranges[7],all_ranges[8],all_ranges[9],all_ranges[10]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i) 
end

function get_rh_init_switch_all_ranges(rtc_model, ranges, branch_ssval, specie, l, params, num_species, all_species)
    res=[];
    init_vals=[];
    unstable=[];
    params1 = deepcopy(params)
    for (range,i) in zip(ranges,range(1,9))
        initial = deepcopy(branch_ssval)
        # res=[]
        for j in range
            initial[i] = j
            # @show initial
            solu = sol(rtc_model, initial, tspan, params1)
            # if solu.retcode == ReturnCode.Unstable
            #     push!(unstable, initial)
            # end
            df = create_solu_df(solu, species_rtc)
            push!(res, get_ssval(df,specie))
            push!(init_vals, initial[i])
        end
        # push!(all_res,res)
    end

    res = Float64.(res)
    res = reshape(res, (l,num_species))
    # init_vals = Float64.init_vals
    init_vals = reshape(init_vals, (l,num_species))
    # return DataFrame(res,all_species), DataFrame(init_vals,all_species), unstable
    # @show size(res)
    # return DataFrame(res,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]), DataFrame(init_vals,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]);
    if num_species == 10
        return DataFrame(res,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtc_i]), DataFrame(init_vals,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtc_i]);
    else 
        return DataFrame(res,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]), DataFrame(init_vals,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]);
    end

end

function upper_or_lower(df, lower_branch, l, num_species)
    arr=[];
    # state1 = "$state"
    for col in eachcol(df)
        for i in col
        # for j in i 
            # if round(j;digits=3) == round(lower_branch[7];digits=3)
            if lower_branch-(0.001*lower_branch) < i < lower_branch+(0.001*lower_branch) 
                push!(arr, 0)
            else
                push!(arr, 1)
            end
            # end
        end
    end
    arr = reshape(arr, (l,num_species))
    if num_species == 10
        df = DataFrame(arr,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i])
    else
        df = DataFrame(arr,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    end 
    return df;
end




function get_switch_ind(binary_df,l) # for off to on l needs to be set to 0, otherwise l is l 
    switch_ind=[]
    for col in eachcol(binary_df)
        # @show length(findall(x->x==round(branches1.ss_val_on[7];digits=3),col))
        # if length(findall(x->x==round(branches1.ss_val_on[7];digits=3),col)) == 0
        if length(findall(x->x==1,col)) == l
            push!(switch_ind,NaN)
        else
            # push!(switch_ind,findall(x->x==round(branches1.ss_val_on[7];digits=3),col)[1])
            push!(switch_ind,findall(x->x==1,col)[1])

        end
    end
    return switch_ind
end

function get_switch_vals(switch_ind,init_vals)
    switch_vals=[]
    for (ind,val) in zip(switch_ind,eachcol(init_vals))
        if isnan(ind) == true
            push!(switch_vals,NaN)
        else
            push!(switch_vals,val[ind])
        end
    end
    return switch_vals
end

function plot_rtc_bf(df, kdam1, kdam2, specie, legendgroup, colour, name)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df[!,specie][1:kdam1], name=name, line=attr(width=6.5, color=colour), showlegend=true, legendgroup=legendgroup)#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df[!,specie][kdam1:kdam2], name="", mode="lines", line=attr(width=6.5,dash="dash", color=colour),showlegend=false, legendgroup=legendgroup)
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df[!,specie][kdam2:end], name="", line=attr(width=6.5, color=colour),showlegend=false, legendgroup=legendgroup)
    return rtcb1, rtcb2, rtcb3
end