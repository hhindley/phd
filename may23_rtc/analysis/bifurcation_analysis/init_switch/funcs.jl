function setup_ssvals(params)
    tspan = (0,1e9)
    # initial =  [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]
    initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

    solu = sol(rtc_model, initial, tspan, params)
    # plotly_plot_sol(solu, "", "", "")
    ss_init = ss_init_vals(solu)
    # soluss = sol(rtc_model, ss_init, tspan, params)
    # plotly_plot_sol(soluss, "", "", "")

    initial1 = [0., 0., 1., 0., 0., 0., 11.29, 0., 0.]
    solu1 = sol(rtc_model, initial1, tspan, params)
    ss_init_upper = ss_init_vals(solu1)

    #define ss for each branch, remember that rt and rd are going to be opposite to the rest of the species
    branches = DataFrame(species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"],ss_val_off=[ss_init[1],ss_init[2],ss_init[3],ss_init[4],ss_init[5],ss_init[6],ss_init[7],ss_init[8],ss_init[9]],ss_val_on=[ss_init_upper[1],ss_init_upper[2],ss_init_upper[3],ss_init_upper[4],ss_init_upper[5],ss_init_upper[6],ss_init_upper[7],ss_init_upper[8],ss_init_upper[9]],parameter=["ATP","kin","ω_ab","ω_r","λ","kdam",NaN,NaN,NaN], param_val=[params.atp,params.kin,params.ω_ab,params.ω_r,params.lam,params.kdam,NaN, NaN, NaN])
    return branches
end


function set_ss_range_zerotoNssval(branch_df, branch, specie, n, l)
    a = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]-(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),branch_df[branch_df.species .== specie,:][:,Symbol("$branch")],length=Int(l/2))
    b = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]+(n*branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),length=Int(l/2))
    return [(vcat(a,b)...)...]
end

function set_ss_range_zerotossval(branch_df, branch, specie, n, l)
    b = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]-branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),branch_df[branch_df.species .== specie,:][:,Symbol("$branch")],length=l)
    return [(b...)...]
end

function get_all_ranges(func, branch_df, branch, n, l)
    all_ranges=[]
    # species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"]
    for i in all_species
        push!(all_ranges, func(branch_df, branch, i, n, l))
    end
    return @LArray [all_ranges[1], all_ranges[2],all_ranges[3],all_ranges[4],all_ranges[5],all_ranges[6],all_ranges[7],all_ranges[8],all_ranges[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt) 
end


function get_ss_change_init(init, range, l, params)
    initial = deepcopy(init)
    res=[]
    for specie in all_species
        for i in range
            initial[7] = i
            @show initial
            solu = sol(rtc_model, initial, tspan, params)
            push!(res, get_ssval(solu,specie))
            # push!(res, [get_ssval(solu,specie) for specie in all_species])
        end
    end
    # @show res
    res = Float64.(res)
    res = (reshape(res, (l,9)))
    # @show typeof(res)
    df = DataFrame(res,all_species)
    return df
end


function get_rh_init_switch_all_ranges(ranges, branch_ssval, specie, l, params)
    res=[]
    init_vals=[]
    for (range,i) in zip(ranges,range(1,9))
        initial = deepcopy(branch_ssval)
        # res=[]
        for j in range
            initial[i] = j
            # @show initial
            solu = sol(rtc_model, initial, tspan, params)
            push!(res, get_ssval(solu,specie))
            push!(init_vals, initial[i])
        end
        # push!(all_res,res)
    end

    res = Float64.(res)
    res = (reshape(res, (l,9)))
    # init_vals = Float64.init_vals
    init_vals = reshape(init_vals, (l,9))
    return DataFrame(res,all_species), DataFrame(init_vals,all_species)
end




# function get_all_init_switch_result(ranges, init)
#     all=[]
#     for range in ranges
#         push!(all, get_ss_change_init(init, range))
#     end
#     return all
# end

function upper_or_lower(df, lower_branch, l)
    arr=[]
    # state1 = "$state"
    for col in eachcol(df)
        for i in col
            for j in i 
                # if round(j;digits=3) == round(lower_branch[7];digits=3)
                if lower_branch[7]-(0.001*lower_branch[7]) < j < lower_branch[7]+(0.001*lower_branch[7]) 
                    push!(arr, 0)
                else
                    push!(arr, 1)
                end
            end
        end
    end

    # if state1 == "lower"
    #     for col in eachcol(df)
    #         for i in col
    #             for j in i 
    #                 if round(j;digits=10) == round(col[1];digits=10)
    #                     push!(arr, 0)
    #                 else
    #                     push!(arr, 1)
    #                 end
    #             end
    #         end
    #     end
    # else
    #     for col in eachcol(df)
    #         for i in col
    #             for j in i 
    #                 if round(j;digits=10) == round(col[1];digits=10)
    #                     push!(arr, 1)
    #                 else
    #                     push!(arr, 0)
    #                 end
    #             end
    #         end
    #     end
    # end
    arr = reshape(arr, (l,9))
    df = DataFrame(arr,all_species)
    return df
end

# function all_upper_lower(all, state)
#     dfs=[]
#     for df in all
#         push!(dfs, upper_or_lower(df, state))
#     end
#     return dfs
# end

function set_shared_range(n,l)
    x1 = range(-100,0, length=Int(l/2))
    x2 = range(0,n*100, length=Int(l/2))
    return vcat(x1,x2)
end
# function set_shared_range(n)
#     # x1 = range(-100,0, length=200)
#     x2 = range(0,n*100, length=500)
#     return x2
# end


function double_init(branch,ranges,opposite_branch,init1,init2,params)
    res=[]
    sss=[]
    rtcbs=[]
    rtcrs=[]
    initial = @LArray [branch[1],branch[2],branch[3],branch[4],branch[5],branch[6],branch[7],branch[8],branch[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt)
    for k in ranges[Symbol("$init1")]
        initial[Symbol("$init1")] = k
        for j in ranges[Symbol("$init2")]
            initial[Symbol("$init2")] = j
            solu = sol(rtc_model, Array(initial), tspan, params)
            ss = get_ssval(solu,:rh)
            push!(sss,ss)
            push!(rtcbs, initial[4])
            push!(rtcrs,initial[6])
            if ss == opposite_branch[7]
                push!(res, (initial[4],initial[6]))
            end
            @show initial    
        end
    end
    return sss, rtcbs, rtcrs
end

function get_binary(lower_branch, sss)
    binary=[]
    for i in sss
        if lower_branch[7]-(0.001*lower_branch[7]) < i < lower_branch[7]+(0.001*lower_branch[7]) 
            push!(binary, 0)
        else
            push!(binary, 1)
        end
    end
    return binary
end

function get_xy_percentage(shared_range)
    rtcbs1=[]
    rtcrs1=[]
    for i in shared_range
        for j in shared_range
            push!(rtcbs1,i)
            push!(rtcrs1,j)
        end
    end
    return rtcbs1, rtcrs1
end

function plot_binary(sss, shared_range, lower_branch,init1, init2, title)
    binary = get_binary(lower_branch, sss)
    rtcbs1, rtcrs1 = get_xy_percentage(shared_range)
    return plot(contour(x=rtcbs1, y=rtcrs1, z=binary, contours_start=0, contours_end=1, contours_size=1, colorscale=[[0,"aqua"],[1,"blue"]]),
    Layout(xaxis_title="$init1 - % from ss_val", yaxis_title="$init2 - % from ss_val", title=title))
end


function triple_init(branch,ranges,opposite_branch,init1,init2,init3,params)
    res=[]
    sss=[]
    init1s=[]
    init2s=[]
    init3s=[]
    initial = @LArray [branch.ss_val[1],branch.ss_val[2],branch.ss_val[3],branch.ss_val[4],branch.ss_val[5],branch.ss_val[6],branch.ss_val[7],branch.ss_val[8],branch.ss_val[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt)
    for i in ranges[Symbol("$init3")]
        initial[Symbol("$init3")] = i
        for k in ranges[Symbol("$init1")]
            initial[Symbol("$init1")] = k
            for j in ranges[Symbol("$init2")]
                initial[Symbol("$init2")] = j
                solu = sol(rtc_model, Array(initial), tspan, params)
                ss = get_ssval(solu,:rh)
                push!(sss,ss)
                push!(init1s, initial[Symbol("$init1")])
                push!(init2s,initial[Symbol("$init2")])
                push!(init3s,initial[Symbol("$init3")])
                if ss == opposite_branch.ss_val[7]
                    push!(res, (initial[4],initial[6]))
                end
                @show initial    
            end
        end
    end
    return sss, init1s, init2s, init3s
end





function round_ssvals(all)
    rounded=[]
    for col in eachcol(all)
        new_rh=[]
        for i in col
            push!(new_rh,round(i;digits=3))
        end
        push!(rounded,new_rh)
    end
    return DataFrame(rounded,all_species)
end

function get_switch_ind(binary_df)
    switch_ind=[]
    for col in eachcol(binary_df)
        # @show length(findall(x->x==round(branches1.ss_val_on[7];digits=3),col))
        # if length(findall(x->x==round(branches1.ss_val_on[7];digits=3),col)) == 0
        if length(findall(x->x==1,col)) == 0
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

function get_percentages(switch_vals,branch) #used to check the original plots
    #work out percentage different between switch_vals and ss_val_lower to see if it corresponds to plot already made
    perc=[]
    for i in range(1,9)
        if isnan(switch_vals[i]) == true
            push!(perc, NaN)
        else
            push!(perc,100*((switch_vals[i]-branch.ss_val_off[i])/branch.ss_val_off[i]))
        end
    end
    return perc
end

function get_diffs(switch_vals,branch)
    #plot how much more than ss_val is needed in uM 
    diffs=[]
    for i in range(1,9)
        if isnan(switch_vals[i]) == true
            push!(diffs,NaN)
        else
            push!(diffs,(switch_vals[i]-branch.ss_val_off[i]))
        end
    end
    return diffs
end

function full_find_differences_or_percs(all,func,init_vals,branches1)
    binary_df = upper_or_lower(all,branches1.ss_val_off,l)

    switch_ind = get_switch_ind(binary_df)

    switch_vals = get_switch_vals(switch_ind,init_vals)

    diffs = func(switch_vals,branches1)
    return diffs

end 



function setup_ssvals_from_bfkit(kdam_val)
    params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014)
    initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]
    br2 = get_br(rtc_mod, params2, initial, 3.)
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

    return DataFrame(species=all_species,ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val)],ss_val_off=[int_rma2(kdam_val),int_rtca2(kdam_val),int_rmb2(kdam_val),int_rtcb2(kdam_val),int_rmr2(kdam_val),int_rtcr2(kdam_val),int_rh2(kdam_val),int_rd2(kdam_val),int_rt2(kdam_val)])
end