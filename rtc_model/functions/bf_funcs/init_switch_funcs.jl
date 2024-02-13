

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


function set_ss_range_Nssval(branch_df, branch, specie, n, l)
    # a = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]-(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),branch_df[branch_df.species .== specie,:][:,Symbol("$branch")],length=Int(l/2))
    b = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]+(n*branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),length=l)
    # return [(vcat(a,b)...)...]
    return [(b...)...]
end

function set_ss_range_zerotossval(branch_df, branch, specie, n, l)
    b = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]-branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),length=l)
    return [(b...)...]
end

function set_ss_range_full(branch_df, branch, specie, n, l)
    a = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]-(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),branch_df[branch_df.species .== specie,:][:,Symbol("$branch")],length=Int(l/2))
    b = range((branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]),(branch_df[branch_df.species .== specie,:][:,Symbol("$branch")]+(n*branch_df[branch_df.species .== specie,:][:,Symbol("$branch")])),length=Int(l/2))
    # return [(vcat(a,b)...)...]
    return [(vcat(a,b)...)...]
end

function get_all_ranges(func, branch_df, branch, n, l)
    all_ranges=[]
    # species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"]
    for i in species_rtc#[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i]
        push!(all_ranges, func(branch_df, branch, i, n, l))
    end
    return @LArray [all_ranges[1], all_ranges[2],all_ranges[3],all_ranges[4],all_ranges[5],all_ranges[6],all_ranges[7],all_ranges[8],all_ranges[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt)  #[all_ranges[1], all_ranges[2],all_ranges[3],all_ranges[4],all_ranges[5],all_ranges[6],all_ranges[7],all_ranges[8],all_ranges[9],all_ranges[10]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i) 
end


function get_ss_change_init(init, range, l, params)
    initial = deepcopy(init)
    res=[]
    for specie in all_species
        for i in range
            initial[7] = i
            # @show initial
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
            push!(res, get_ssval(solu,specie,all_species))
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




# function get_all_init_switch_result(ranges, init)
#     all=[]
#     for range in ranges
#         push!(all, get_ss_change_init(init, range))
#     end
#     return all
# end

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
function set_shared_range_0ton(n,l)
    # x1 = range(-100,0, length=200)
    x2 = range(0,n*100, length=l)
    return x2
end


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
            # if ss == opposite_branch[7]
            #     push!(res, (initial[4],initial[6]))
            # end
            # @show initial    
        end
    end
    return sss, rtcbs, rtcrs
end

function get_binary(lower_branch, sss)
    binary=[]
    for i in sss
        if lower_branch-(0.1*lower_branch) < i < lower_branch+(0.1*lower_branch) 
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
    initial = @LArray [branch[1],branch[2],branch[3],branch[4],branch[5],branch[6],branch[7],branch[8],branch[9]] (:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt)
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
                # if ss == opposite_branch.ss_val[7]
                #     push!(res, (initial[4],initial[6]))
                # end
                # @show initial    
            end
        end
    end
    return sss, init1s, init2s, init3s
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

function get_percentages(switch_vals,branch,branch_label) #used to check the original plots
    #work out percentage different between switch_vals and ss_val_lower to see if it corresponds to plot already made
    perc=[];
    if branch_label == "off"
    for i in range(1,9)
        if isnan(switch_vals[i]) == true
            push!(perc, NaN)
        else
            push!(perc,100*((switch_vals[i]-branch[i])/branch[i]))
        end
    end
    else
        for i in range(1,9)
            if isnan(switch_vals[i]) == true
                push!(perc, NaN)
            else
                push!(perc,100*((branch[i]-switch_vals[i])/branch[i]))
            end
        end
    end
    return perc;
end

function get_diffs(switch_vals,branch,branch_label)
    #plot how much more than ss_val is needed in uM 
    diffs=[]
    # if branch == branches1.ss_val_off
    for i in range(1,9)
        if isnan(switch_vals[i]) == true
            push!(diffs,NaN)
        else
            push!(diffs,(switch_vals[i]-branch[i]))
        end
    end
    # else
    #     for i in range(1,9)
    #         if isnan(switch_vals[i]) == true
    #             push!(diffs,NaN)
    #         else
    #             push!(diffs,(branch[i]-switch_vals[i]))
    #         end
    #     end
    # end
    return diffs
end

function get_multiples(switch_vals,branch,branch_label)
    #plot how much more than ss_val is needed in uM 
    diffs=[]
    if branch_label == "off"
        for i in range(1,9)
            if isnan(switch_vals[i]) == true
                push!(diffs,NaN)
            else
                push!(diffs,(switch_vals[i]/branch[i]))
            end
        end
    else
        for i in range(1,9)
            if isnan(switch_vals[i]) == true
                push!(diffs,NaN)
            else
                push!(diffs,(branch[i]/switch_vals[i]))
            end
        end
    end 
    return diffs
end

function full_find_differences_or_percs(all,func,init_vals,lower_branch,l,branch,l1,branch_label)
    binary_df = upper_or_lower(all,lower_branch,l,9)

    switch_ind = get_switch_ind(binary_df,l1)

    switch_vals = get_switch_vals(switch_ind,init_vals)

    diffs = func(switch_vals,branch,branch_label)
    return diffs
    # return switch_vals
end 

function split_curves(df, df_bf)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    middle=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    last=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    
    for i in df.kdam[1:kdam1]
        push!(first.kdam, i)
    end
    for i in df.kdam[kdam1:kdam2]
        push!(middle.kdam, i)
    end
    for i in df.kdam[kdam2:end]
        push!(last.kdam, i)
    end
    for (col,col1) in zip(eachcol(df)[1:9],eachcol(first)[2:end])
        for i in col[1:kdam1]
            push!(col1, i)
        end
    end
    for (col,col1) in zip(eachcol(df)[1:9],eachcol(middle)[2:end])
        for i in col[kdam1:kdam2]
            push!(col1, i)
        end
    end
    for (col,col1) in zip(eachcol(df)[1:9],eachcol(last)[2:end])
        for i in col[kdam2:end]
            push!(col1, i)
        end
    end
    return first, middle, last
end

function setup_ssvals_from_bfkit(rtc_mod, kdam_val, params2, initial, kdam)
    # params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    # θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    # krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    # kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    # kdam =  0.01, lam = 0.014)
    br2 = get_br(rtc_mod, params2, initial, kdam)
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

function setup_ssvals_from_bfkit_inhib(rtc_mod, kdam_val, params2)
    # params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    # θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    # krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    # kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    # kdam =  0.01, lam = 0.014)
    initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0., 0.]
    br2 = get_br(rtc_mod, params2, initial, 3.)
    df = create_br_df_inhib(br2)
    df_bf = bf_point_df_inhib(br2)
    first,middle,last=split_curves_inhib(df, df_bf)
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
    int_rtcbi1 = QuadraticInterpolation(first.rtcb_i, first.kdam)
    int_rtcbi2 = QuadraticInterpolation(last.rtcb_i, last.kdam)

    return DataFrame(species=[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i],ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val),int_rtcbi1(kdam_val)],ss_val_off=[int_rma2(kdam_val),int_rtca2(kdam_val),int_rmb2(kdam_val),int_rtcb2(kdam_val),int_rmr2(kdam_val),int_rtcr2(kdam_val),int_rh2(kdam_val),int_rd2(kdam_val),int_rt2(kdam_val),int_rtcbi2(kdam_val)]);
end


function create_resdf(all_res,kdam_range)
    df = DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    for (res,kdam) in zip(all_res,kdam_range)
        push!(df.rm_a,res[1])
        push!(df.rtca,res[2])
        push!(df.rm_b,res[3])
        push!(df.rtcb,res[4])
        push!(df.rm_r,res[5])
        push!(df.rtcr,res[6])
        push!(df.rh,res[7])
        push!(df.rd,res[8])
        push!(df.rt,res[9])
        push!(df.kdam,kdam)
    end
    return df
end




function get_plot_vals(binary, l, n)
    binary_matrix = reshape(binary, (l,l))
    # ind = findall(x->x==1,binary_matrix)
    # tups=[]
    # for i in ind
    #     push!(tups, Tuple(i))
    # end
    # first_tuples_dict = OrderedDict{Int, Tuple{Int, Int}}()

    # for tuple in tups
    #     second_index = tuple[2]
    #     if !haskey(first_tuples_dict, second_index)
    #         # If this is the first tuple with the current second index, store it
    #         first_tuples_dict[second_index] = tuple
    #     end
    # end
    last_tuples_dict = OrderedDict{Int, Tuple{Int, Int}}()
    ind1 = findall(x->x==0,binary_matrix)
    tups1=[]
    for i in ind1
        push!(tups1, Tuple(i))
    end
    for tuple in tups1
        second_index = tuple[2]
        last_tuples_dict[second_index] = tuple  # Update the dictionary with the current tuple
    end

    xs,ys = get_xy_percentage(set_shared_range_0ton(n,l))
    new_arr = reshape(tuple.(xs,ys), (l,l))
    vals = Tuple(values(last_tuples_dict))

    plot_vals=[]
    for i in vals
        push!(plot_vals, new_arr[i[2],i[1]])
    end

    x=[];y=[];
    for i in plot_vals
        push!(x, i[2])
        push!(y, i[1])
    end
    return x, y
end


function double_param_vary(param_range1, param1, param_range2, param2, params1)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            # @show params[:k_inhib], params[:inhib]
            br = get_br(rtc_mod, params, initial, 3.)
            if length(br.specialpoint) == 2
                push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end    
    end
    return df
end

function double_param_vary_inhib(param_range1, param1, param_range2, param2, params1)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            @show params
            br = get_br(rtc_inhib_mod, params, initial_i, 3.)
            if length(br.specialpoint) == 2
                push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end    
    end
    return df
end

function plot_bistable_region(param_range1, param1, param_range2, param2, params1, colour, title)
    df = double_param_vary_inhib(param_range1, param1, param_range2, param2, params1)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
    end
    if length(max_) == length(param_range1)
        p = plot(param_range1, max_; fillrange=(min_), fillalpha = 0.35, fillcolor=colour, linecolor=:white, title=title, legend=false, xlims=(minimum(param_range1), maximum(param_range1)), ylims=(minimum(param_range2), maximum(param_range2)))
        # p = plot!(plotatpt, plotlamt, c=:black)
    else
        Nothing
    end

    return (p)
end


function bistable_region(param_range1, param1, param_range2, param2, params1)
    df = double_param_vary(param_range1, param1, param_range2, param2, params1)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
    end
    return max_, min_
end

function bistable_region_inhib(param_range1, param1, param_range2, param2, params1)
    df = double_param_vary_inhib(param_range1, param1, param_range2, param2, params1)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
    end
    return max_, min_
end

function plot_bs_region_same_plot(param_range1, param_range2, results, title, range1, param1, param2)
    colours =palette(:tab10)
    p = plot()
    for (i,j) in zip(range(1,length(results)), range(1,5))
        if length(results[i][1]) == length(param_range1)
            p = plot!(param_range1, results[i][1]; fillrange=(results[i][2]), fillalpha = 0.45, fillcolor=colours[j], 
            linecolor=colours[j], title=title, label="$(@sprintf "%g" (range1[j]))", xlims=(minimum(param_range1), maximum(param_range1)), 
            ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2")
            # display(p)
        else
            Nothing
        end
    end
    return p 
    # return p = (plot!(plot_var1, plot_var2, c=:black, label=""))
end

function get_bs_region_results_wr(param_range1, param1, param_range2, param2, wab)
    results_wr=[]
    for i in wr_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = wab, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wr, (max_,min_))
    end
    return results_wr
end

function get_bs_region_results_wab(param_range1, param1, param_range2, param2, wr)
    results_wab=[]
    for i in wab_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = i, ω_r = wr, 
        kdam =  0.01, lam = 0.014) 	
        # @show (params1[:ω_ab], params1[:atp])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wab, (max_,min_))
    end
    return results_wab
end

function get_bs_region_results_wab_inhib(param_range1, param1, param_range2, param2, wr)
    results_wab=[]
    for i in wab_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = i, ω_r = wr, 
        kdam =  0.01, lam = 0.014, k_inhib=0.1, inhib=10)	
        # @show (params1[:ω_ab], params1[:atp])
        max_,min_ = bistable_region_inhib(param_range1, param1, param_range2, param2, params1)
        push!(results_wab, (max_,min_))
    end
    return results_wab
end