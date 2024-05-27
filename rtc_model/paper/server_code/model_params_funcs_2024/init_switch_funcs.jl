

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

function get_all_ranges(func, branch_df, branch, n, l)
    all_ranges=[]
    # species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"]
    for i in all_species#[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i]
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
        return DataFrame(res,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i]), DataFrame(init_vals,[:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtcb_i]);
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
            if lower_branch-(0.1*lower_branch) < i < lower_branch+(0.1*lower_branch) 
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

function full_find_differences_or_percs(all,func,init_vals,lower_branch,l,branch,l1,branch_label)
    binary_df = upper_or_lower(all,lower_branch,l,9)

    switch_ind = get_switch_ind(binary_df,l1)

    switch_vals = get_switch_vals(switch_ind,init_vals)

    diffs = func(switch_vals,branch,branch_label)
    return diffs
    # return switch_vals
end 



function setup_ssvals_from_bfkit(rtc_mod, kdam_val, params2, initial)
    # params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    # θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    # krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    # kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    # kdam =  0.01, lam = 0.014)
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

    return DataFrame(species=all_species,ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val)],ss_val_off=[int_rma2(kdam_val),int_rtca2(kdam_val),int_rmb2(kdam_val),int_rtcb2(kdam_val),int_rmr2(kdam_val),int_rtcr2(kdam_val),int_rh2(kdam_val),int_rd2(kdam_val),int_rt2(kdam_val)]);
end

function numerical_bistability_analysis(model, params, init, specie, all_species, kdam_range) # used to be 'checking_bistability' - used when bifurcationkit not working 
    param_init = deepcopy(params)
    new_params = deepcopy(params)
    first_params = deepcopy(params)
    first_params[:kdam]=kdam_range[1]
    solu = sol(model, init, tspan, first_params)
    ss = get_ssval(solu, specie, all_species)
    init_first = ss_init_vals(solu, all_species)
    res =[]
    for i in range(2, length(kdam_range))
        param_init[:kdam]=kdam_range[i-1]
        solu_init = sol(model, init_first, tspan, param_init)
        init_ss = ss_init_vals(solu_init, all_species)
        new_params[:kdam] = kdam_range[i]
        solu_new = sol(model, init_ss, tspan, new_params)
        push!(res, get_ssval(solu_new, specie, all_species))
    end
    pushfirst!(res, ss)
    return res
end

function ss_init_vals(sol, species)
    if length(species) == 9
        rm_a = get_ssval(sol, :rm_a, species)
        rm_b = get_ssval(sol, :rm_b, species)
        rm_r = get_ssval(sol, :rm_r, species)
        rtca = get_ssval(sol, :rtca, species)
        rtcb = get_ssval(sol, :rtcb, species)
        rtcr = get_ssval(sol, :rtcr, species)
        rt = get_ssval(sol, :rt, species)
        rd = get_ssval(sol, :rd, species)
        if species[7] == :rh
            rh = get_ssval(sol, :rh, species)
            return [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt]
        else
            trna = get_ssval(sol, :trna, species)
            return [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt]
        end
    else
        rm_a = get_ssval(sol, :rm_a, species)
        rm_b = get_ssval(sol, :rm_b, species)
        rm_r = get_ssval(sol, :rm_r, species)
        rtca = get_ssval(sol, :rtca, species)
        rtcb = get_ssval(sol, :rtcb, species)
        rtcr = get_ssval(sol, :rtcr, species)
        rt = get_ssval(sol, :rt, species)
        rd = get_ssval(sol, :rd, species)
        rtc_i = get_ssval(sol, :rtc_i, species)
        if species[7] == :rh
            rh = get_ssval(sol, :rh, species)
            return [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i]
        else
            trna = get_ssval(sol, :trna, species)
            return [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i]
        end
    end
end

function setup_ssvals_from_bfkit_trna(rtc_mod, kdam_val, params2, initial)
    br2 = get_br(rtc_mod, params2, initial, 400.)
    df = create_br_df(br2)
    bf = bf_point_df(br2)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    df = Float64.(df)

    
    int_rma1 = QuadraticInterpolation(df.rm_a[1:kdam1], df.kdam[1:kdam1])
    int_rtca1 = QuadraticInterpolation(df.rtca[1:kdam1], df.kdam[1:kdam1])
    int_rmb1 = QuadraticInterpolation(df.rm_b[1:kdam1], df.kdam[1:kdam1])
    int_rtcb1 = QuadraticInterpolation(df.rtcb[1:kdam1], df.kdam[1:kdam1])
    int_rmr1 = QuadraticInterpolation(df.rm_r[1:kdam1], df.kdam[1:kdam1])
    int_rtcr1 = QuadraticInterpolation(df.rtcr[1:kdam1], df.kdam[1:kdam1])
    int_rh1 = QuadraticInterpolation(df.rh[1:kdam1], df.kdam[1:kdam1])
    int_rd1 = QuadraticInterpolation(df.rd[1:kdam1], df.kdam[1:kdam1])
    int_rt1 = QuadraticInterpolation(df.rt[1:kdam1], df.kdam[1:kdam1])



    return DataFrame(species=all_species,ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val)]);
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

using BifurcationKit
const BK = BifurcationKit


# sup norm
norminf(x) = norm(x, Inf)

# parameter values
# params = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
# θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
# krep = 137., ktag = 9780., km_a = 20., km_b = 16., g_max = 2.0923, kdeg = 0.001, 
# kdam =  0.01,
# ω_ab = 2., ω_r = 0.0089, atp = 3000., kin = 0.022222, lam = 0.04)

params_bf = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014) 		
#ω_ab = 1, ω_r = 0.0001


initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

initial1 = [0., 0., 1., 0., 0., 0., 11.29, 0., 0.]


function get_br(model, params, initial, kdam_max)
    # Bifurcation Problem
    if length(params) == 25 
        prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001, 
        dsmax = 0.15, # 0.001, 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 4, max_bisection_steps = 10, #3,2,10
        # number of eigenvalues
        nev = 2,  #4 
        # maximum number of continuation steps
        max_steps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution
    elseif length(params) == 28
        prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001, 
        dsmax = 0.05, # 0.001, 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 4, max_bisection_steps = 10, #3,2,10
        # number of eigenvalues
        nev = 2,  #4 
        # maximum number of continuation steps
        max_steps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution

    else
        prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        # continuation options
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001, dsmax = 0.01,
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 8, max_bisection_steps = 25,
        # number of eigenvalues
        nev = 2, 
        # maximum number of continuation steps
        max_steps = 1000,) # only using parameters that make a difference to solution

    end
    # continuation of equilibria
    br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    # br = continuation(prob, PALC(θ=0.75), opts_br; plot = false, bothside=true, normC = norminf)
    return br
end


function sol(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4())
    return solu
end

all_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]



function get_ssval(sol, species, all_species)
    df = create_solu_df(sol, all_species)
    species = df[end, species]
    return species
end



function create_br_df(br)
    if length(br.sol[1][1]) == 9
        df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],kdam=[]);
    elseif length(br.sol[1][1]) == 10
        df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[],kdam=[]);
    elseif length(br.sol[1][1]) == 11
        df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtca_i=[],rtcb_i=[],kdam=[]);
    else
        Nothing
    end
    # for i in eachcol(df)
    #     println(i)
    # end
    for i in range(1,length(br.sol))
        for (s,d) in zip(range(1,length(br.sol[1][1])), eachcol(df))
            push!(d, br.sol[i][1][s])
        end
        push!(df.kdam, br.sol[i][2])
    end
    
    return df;
end

function bf_point_df(br2)
    df_bf = DataFrame(rm_a=Float64[],rtca=Float64[],rm_b=Float64[],rtcb=Float64[],rm_r=Float64[],rtcr=Float64[],rh=Float64[],rd=Float64[],rt=Float64[],kdam=Float64[]);
    # df_bf.rm_a;
    for i in br2.specialpoint
        if i.type == :bp
            push!(df_bf.rm_a, i.x[1])
            push!(df_bf.rtca, i.x[2])
            push!(df_bf.rm_b, i.x[3])
            push!(df_bf.rtcb, i.x[4])
            push!(df_bf.rm_r, i.x[5])
            push!(df_bf.rtcr, i.x[6])
            push!(df_bf.rh, i.x[7])
            push!(df_bf.rd, i.x[8])
            push!(df_bf.rt, i.x[9])
            push!(df_bf.kdam, i.param)
            # [push!(col,i.x[num] for (num, col) in zip(range(1,9), eachcol(df_bf)))]
        
        end
    end
    return df_bf;
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

function create_solu_df(sol, all_species)
    df = DataFrame(sol)
    names = [:time;all_species]
    rename!(df, names)
    return df
end