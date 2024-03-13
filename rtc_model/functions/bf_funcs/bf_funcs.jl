using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

function get_br(model, init, params, kdam_max)
    prob = ODEProblem(model, init, tspan, params; jac=true)
    odefun = prob.f
    F = (u,p) -> odefun(u,p,0)
    J = (u,p) -> odefun.jac(u,p,0)
    id_kdam = indexof(kdam, parameters(model))
    par_tm = prob.p
    # Bifurcation Problem
    if nameof(model) == :rtc_model || nameof(model) == :rtc_inhib_model
        # print("original rtc model")
        prob = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); J=J,
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001, a=0.1,
        dsmax = 0.05, # 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 4, max_bisection_steps = 20, #3,2,10
        # number of eigenvalues
        nev = 2, 
        # maximum number of continuation steps
        max_steps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution
        # continuation of equilibria
        br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    else
        # print("tRNA rtc model")
        prob = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); J=J,
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001,# a=0.1,
        dsmax = 0.15, dsmin = 0.0001,# 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 20, #3,2,10
        # number of eigenvalues
        # nev =100, #tolParamBisectionEvent=1e-30, 
        # maximum number of continuation steps
        max_steps = 50000,)# dsmin_bisection=1e-30)#, tol_bisection_eigenvalue=1e-10)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution
        # continuation of equilibria

        br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    end
    return br
end



function numerical_bistability_analysis(model, params, init, specie, all_species, kdam_range, kdam) # used to be 'checking_bistability' - used when bifurcationkit not working 
    param_init = deepcopy(params)
    new_params = deepcopy(params)
    first_params = deepcopy(params)
    first_params[kdam]=kdam_range[1]
    solu = sol(model, init, tspan, first_params)
    df_sol = create_solu_df(solu, all_species)
    if specie == :lam
        ss_vals = get_all_ssvals(solu, all_species)
        ssvals_dict = Dict([i => j for (i,j) in zip(all_species, ss_vals)])
        ss = calc_lam(first_params, ssvals_dict)
    else
        ss = get_ssval(df_sol, specie)
    end

    init_first = ss_init_vals(df_sol, all_species)
    res =[]
    for i in ProgressBar(range(2, length(kdam_range)))
        param_init[kdam]=kdam_range[i-1]
        solu_init = sol(model, init_first, tspan, param_init)
        df_sol_init = create_solu_df(solu_init, all_species)
        init_ss = ss_init_vals(df_sol_init, all_species)
        new_params[kdam] = kdam_range[i]
        solu_new = sol(model, init_ss, tspan, new_params)
        df_sol_new = create_solu_df(solu_new, all_species)
        if specie == :lam
            ss_vals_new = get_all_ssvals(solu_new, all_species)
            ssvals_dict_new = Dict([i => j for (i,j) in zip(all_species, ss_vals_new)])
            push!(res, calc_lam(new_params,ssvals_dict_new))
        else
            push!(res, get_ssval(df_sol_new, specie))
        end
    end
    pushfirst!(res, ss)
    return res
end

# function numerical_bistability_analysis(model, params, init1, specie, all_species, kdam_range, kdam) # used to be 'checking_bistability' - used when bifurcationkit not working 
#     param_init = deepcopy(params)
#     new_params = deepcopy(params)
#     first_params = deepcopy(params)
#     first_params[kdam]=kdam_range[1]
#     init = steady_states(model, init1, first_params)
#     ss_vals = Dict(zip(all_species, init))
#     res =[]
#     for i in ProgressBar(range(2, length(kdam_range)))
#         param_init[kdam]=kdam_range[i-1]
#         init_ss = steady_states(model, init, param_init)
#         new_params[kdam] = kdam_range[i]
#         final_ss = Dict(zip(all_species, steady_states(model, init_ss, new_params)))
#         push!(res, final_ss[specie])
#     end
#     pushfirst!(res, ss_vals[specie])
#     return res
# end

function full_numerical_bistab(model, params, init1, specie, all_species, kdam_range, kdam_range_rev, kdam)
    res = numerical_bistability_analysis(model, params, init1, specie, all_species, kdam_range, kdam)
    res2 = numerical_bistability_analysis(model, params, init1, specie, all_species, kdam_range_rev, kdam)
    return res, res2
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

    for i in range(1,length(br.sol))
        for (s,d) in zip(range(1,length(br.sol[1][1])), eachcol(df))
            push!(d, br.sol[i][1][s])
        end
        push!(df.kdam, br.sol[i][2])
    end
    
    return df;
end



function bf_point_df_inhib(br2)
    df_bf = DataFrame(rm_a=Float64[],rtca=Float64[],rm_b=Float64[],rtcb=Float64[],rm_r=Float64[],rtcr=Float64[],rh=Float64[],rd=Float64[],rt=Float64[],rtcb_i=Float64[],kdam=Float64[]);
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
            push!(df_bf.rtcb_i, i.x[10])
            push!(df_bf.kdam, i.param)
            # [push!(col,i.x[num] for (num, col) in zip(range(1,9), eachcol(df_bf)))]
        
        end
    end
    return df_bf;
end




function create_br_df_inhib(br)
    df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[],kdam=[]);
    for i in range(1,length(br.sol))
        for (s,d) in zip(range(1,10), eachcol(df))
            push!(d, br.sol[i][1][s])
        end
        push!(df.kdam, br.sol[i][2])
    end
    
    return df;
end


function different_levels_inhibition(rtc_inhib_mod, init_inhib, params_br_inhib, k_inhib1_val, kdam_max)
    params = deepcopy(params_br_inhib)
    params[k_inhib1] = k_inhib1_val
    br = get_br(rtc_inhib_mod, init_inhib, params, kdam_max)
    bf = bf_point_df_inhib(br)
    df = create_br_df_inhib(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    return bf, df, kdam1, kdam2
end

function plot_rtc_bf(df, kdam1, kdam2, specie, legendgroup, colour, name)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df[!,specie][1:kdam1], name=name, line=attr(width=6.5, color=colour), showlegend=true, legendgroup=legendgroup)#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df[!,specie][kdam1:kdam2], name="", mode="lines", line=attr(width=6.5,dash="dash", color=colour),showlegend=false, legendgroup=legendgroup)
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df[!,specie][kdam2:end], name="", line=attr(width=6.5, color=colour),showlegend=false, legendgroup=legendgroup)
    return rtcb1, rtcb2, rtcb3
end

function plot_rtc_bf_init(df, kdam1, kdam2, specie, legendgroup)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df[!,specie][1:kdam1], name="", line=attr(width=6.5, color="#117733ff"), showlegend=true, legendgroup=legendgroup)#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df[!,specie][kdam1:kdam2], name="", mode="lines", line=attr(width=6.5,dash="dash", color=:black),showlegend=false, legendgroup=legendgroup)
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df[!,specie][kdam2:end], name="", line=attr(width=6.5, color="#882255ff"),showlegend=false, legendgroup=legendgroup)
    return rtcb1, rtcb2, rtcb3
end
# functions for the double param vary to produce fig3 banana plot 

function double_param_vary(rtc_model, ssvals_rtc, param_range1, param1, param_range2, param2, params_bf, kdam_max)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    # params = deepcopy(params_bf)
    for i in ProgressBar(param_range1)
        # params_bf = merge(params_bf, (param1=>i,))
        params_bf[param1] = i
        for j in param_range2
            # params_bf = merge(params_bf, (param2=>j,))
            params_bf[param2] = j
            # @show params_bf[:ω_ab], params_bf[:ω_r], params_bf[:atp], params_bf[:lam]
            br = get_br(rtc_model, ssvals_rtc, params_bf, kdam_max)
            if length(br.specialpoint) == 2
                push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end    
    end
    return df
end

function bistable_region(rtc_model, ssvals_rtc, param_range1, param1, param_range2, param2, params_bf, kdam_max)
    df = double_param_vary(rtc_model, ssvals_rtc, param_range1, param1, param_range2, param2, params_bf, kdam_max)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    param_vals = []
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(param_vals, bsp[bsp[:,1] .== i, :][:,1][1])
        end
    end
    return max_, min_, param_vals
end
function get_bs_region_results(rtc_model, ssvals_rtc, param_range1, param1, param_range2, param2, param_range3, param3, params_bf, kdam_max)
    results_wr=[]
    xvals=[]
    params_new = deepcopy(params_bf)
    for i in (param_range3)
        # params_new = merge(params_new, (param3=>i,))
        params_new[param3] = i
        # @show (params_new[param3])
        max_,min_,param_vals = bistable_region(rtc_model, ssvals_rtc, param_range1, param1, param_range2, param2, params_new, kdam_max)
        push!(results_wr, (max_,min_))
        push!(xvals, param_vals)
    end


    return results_wr, xvals
end

function plot_bs_region_same_plot(xvals, param_range1, param_range2, results, title, range1, param1, param2)
    colours = ["#f4e5ffff","#e6c5ffff","#d7a1ffff","#c06affff","#a730ffff"]
    p = Plots.plot()
    for (i,j) in zip(range(1,length(results)), range(1,5))
        if length(results[i][1]) == length(param_range1)
            p = Plots.plot!(xvals[i], results[i][1]; fillrange=(results[i][2]), fillalpha = 0.45, fillcolor=colours[j], 
            linecolor=colours[j], title=title, label="$(@sprintf "%g" (range1[j]))", xlims=(minimum(param_range1), maximum(param_range1)), 
            ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1 (μM)", ylabel="$param2 (min\$^{-1}\$)", tickfontsize=16,
            guidefontsize=20,guidefont="sans-serif",grid=false)
            # display(p)
        else
            Nothing
        end
    end
    return p 
end

function area_under_curve(xvals,yvals)
    x=xvals
    y=yvals
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end
function creating_rtc_inhib_plot(rtc_model, ssvals_rtc, params_rtc, rtc_inhib_mod, ssvals_inhib, params_inhib, specie, kdam_max, colours, k_inhib_vals)
    br = get_br(rtc_model, ssvals_rtc, params_rtc, kdam_max)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    rtcb_01, rtcb_02, rtcb_03 = plot_rtc_bf(df0, kdam01, kdam02, specie, "1", "969696ff", "orig")

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[1], kdam_max)
    rtcb1, rtcb2, rtcb3 = plot_rtc_bf(df, kdam1, kdam2, specie, "2", colours[1], "least inhib")

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[2], kdam_max)
    rtcb1a, rtcb2a, rtcb3a = plot_rtc_bf(dfa, kdam1a, kdam2a, specie, "3", colours[2], "mid inhib")

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[3], kdam_max)
    rtcb1b, rtcb2b, rtcb3b = plot_rtc_bf(dfb, kdam1b, kdam2b, specie, "4", colours[3], "most inhib")
    return [rtcb_01, rtcb_02, rtcb_03, rtcb1, rtcb2, rtcb3, rtcb1a, rtcb2a, rtcb3a, rtcb1b, rtcb2b, rtcb3b]
end


function area_under_curve_rh(df,kdam1)
    df=Float64.(df)
    x = df.kdam[1:2:kdam1]
    y = df.rh[1:2:kdam1]
    # Interpolations.deduplicate_knots!(x, move_knots=false)
    # Interpolations.deduplicate_knots!(y, move_knots=false)
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end

function all_area_under_curve_rh(rtc_inhib_mod, params_bf_inhib, ssvals_inhib, kdam_max, k_inhib_vals, rtc_model, ssvals_rtc, params_rtc)
    br = get_br(rtc_model, ssvals_rtc, params_rtc, kdam_max)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_bf_inhib, k_inhib_vals[1], kdam_max)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_bf_inhib, k_inhib_vals[2], kdam_max)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_bf_inhib, k_inhib_vals[3], kdam_max)

    res0 = area_under_curve_rh(df0,kdam01)
    res = area_under_curve_rh(df,kdam1)
    res1 = area_under_curve_rh(dfa,kdam1a)
    res2 = area_under_curve_rh(dfb,kdam1b)

    return [res0,res1,res2,res]
end


function bf_size(rtc_inhib_mod, ssvals_inhib, params_inhib, kdam_max, k_inhib_vals, rtc_model, ssvals_rtc, params_rtc)

    br = get_br(rtc_model, ssvals_rtc, params_rtc, kdam_max)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[1], kdam_max)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[2], kdam_max)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[3], kdam_max)

    s0 = bf0.kdam[1]-bf0.kdam[2]
    s1 = bf.kdam[1]-bf.kdam[2]
    s2 = bfa.kdam[1]-bfa.kdam[2]
    s3 = bfb.kdam[1]-bfb.kdam[2]

    rtcb_sizes = [s0, s2, s3, s1]
    percentage_of_original_size = [100*(rtcb_sizes[i]/rtcb_sizes[1]) for i in range(2,4)]
    # percentage_decrease_rtcb = [100*((rtcb_sizes[i]-rtcb_sizes[1])/rtcb_sizes[1]) for i in range(2,4)]
    return percentage_of_original_size
    # return percentage_decrease_rtcb
end


function protein_decrease(rtc_inhib_mod, ssvals_inhib, params_inhib, specie, kdam_max, k_inhib_vals, rtc_model, ssvals_rtc, params_rtc)
    br = get_br(rtc_model, ssvals_rtc, params_rtc, kdam_max)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[1], kdam_max)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[2], kdam_max)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, ssvals_inhib, params_inhib, k_inhib_vals[3], kdam_max)

    df=Float64.(df)
    df0=Float64.(df0)
    dfa=Float64.(dfa)
    dfb=Float64.(dfb)

    int_orig = QuadraticInterpolation(df0[!,specie][1:kdam01], df0.kdam[1:kdam01])
    inhib_int_rtcba = QuadraticInterpolation(dfa[!,specie][1:kdam1a], dfa.kdam[1:kdam1a])
    inhib_int_rtcbb = QuadraticInterpolation(dfb[!,specie][1:kdam1b], dfb.kdam[1:kdam1b])
    inhib_int_rtcb1 = QuadraticInterpolation(df[!,specie][1:kdam1], df.kdam[1:kdam1])

    orig_rtcb = [int_orig(i) for i in range(0,df0.kdam[kdam01],length=1000)]
    inhib_rtcba = [inhib_int_rtcba(i) for i in range(0,dfa.kdam[kdam1a],length=1000)]
    inhib_rtcbb = [inhib_int_rtcbb(i) for i in range(0,dfb.kdam[kdam1b],length=1000)]
    inhib_rtcb1 = [inhib_int_rtcb1(i) for i in range(0,df.kdam[kdam1],length=1000)]

    perc_deca = [(i/o)*100 for (i,o) in zip(inhib_rtcba, orig_rtcb)]
    perc_decb = [(i/o)*100 for (i,o) in zip(inhib_rtcbb, orig_rtcb)]
    perc_dec1 = [(i/o)*100 for (i,o) in zip(inhib_rtcb1, orig_rtcb)]

    return [perc_deca, perc_decb, perc_dec1]
end
# function dashed_lines_species(df, df_bf, colors, type)
#     kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
#     kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
#     first=[]
#     middle=[]
#     last=[]
#     names=["RtcBA mRNA","RtcA","RtcB mRNA","RtcB","RtcR mRNA","RtcR"]
#     if type == ""
#         for (col,i) in zip(eachcol(df)[1:6],range(1,9))
#             push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=(names[i]*" $type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
#             push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#             push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#         end
#     else
#         for (col,i) in zip(eachcol(df)[1:6],range(1,9))
#             push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=("$type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
#             push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#             push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#         end
#     end
#     return first, middle, last
# end
# function dashed_lines_ribosomes(df, df_bf, colors, type)
#     kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
#     kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
#     first=[]
#     middle=[]
#     last=[]
#     names=["Healthy ribosomes","Damaged ribosomes","Tagged ribosomes"]
#     if type == ""
#         for (col,i) in zip(eachcol(df)[7:9],range(1,9))
#             push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=(names[i]*" $type"), yaxis="y2", line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
#             push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", yaxis="y2", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#             push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", yaxis="y2", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#         end
#     else
#         for (col,i) in zip(eachcol(df)[7:9],range(1,9))
#             push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=("$type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
#             push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#             push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
#         end
#     end
#     return first, middle, last
# end


# function bf_scatter(df_bf, color)
#     bf_rma = scatter(x=df_bf.kdam, y=df_bf.rm_a, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcBA mRNA")
#     bf_rtca = scatter(x=df_bf.kdam, y=df_bf.rtca, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcA")
#     bf_rmb = scatter(x=df_bf.kdam, y=df_bf.rm_b, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcB mRNA")
#     bf_rtcb = scatter(x=df_bf.kdam, y=df_bf.rtcb, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcB")
#     bf_rmr = scatter(x=df_bf.kdam, y=df_bf.rm_r, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcR mRNA")
#     bf_rtcr = scatter(x=df_bf.kdam, y=df_bf.rtcr, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcR")
#     bf_rh = scatter(x=df_bf.kdam, y=df_bf.rh, mode="markers", yaxis="y2", name="", line=attr(color=color),showlegend=false, legendgroup="Healthy ribosomes")
#     bf_rd = scatter(x=df_bf.kdam, y=df_bf.rd, mode="markers", yaxis="y2", name="", line=attr(color=color),showlegend=false, legendgroup="Damaged ribosomes")
#     bf_rt = scatter(x=df_bf.kdam, y=df_bf.rt, mode="markers", yaxis="y2", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="Tagged ribosomes")
#     return bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr, bf_rh, bf_rd, bf_rt
# end


# function full_lines(df,width, colors)
#     rma_p = scatter(x=df.kdam, y=df.rm_a, name="RtcBA mRNA", line=attr(width=width,color=colors[1]))
#     rtca_p = scatter(x=df.kdam, y=df.rtca, name="RtcA", line=attr(width=width,color=colors[2]))
#     rmb_p = scatter(x=df.kdam, y=df.rm_b, name="rm_b", line=attr(width=width,color=colors[3]))
#     rtcb_p = scatter(x=df.kdam, y=df.rtcb, name="RtcB", line=attr(width=width,color=colors[4]))
#     rmr_p = scatter(x=df.kdam, y=df.rm_r, name="rm_r", line=attr(width=width,color=colors[5]))
#     rtcr_p = scatter(x=df.kdam, y=df.rtcr, name="RtcR", line=attr(width=width,color=colors[6]))
#     rh_p = scatter(x=df.kdam, y=df.rh, name="Rh", yaxis="y2", line=attr(width=width,color=colors[7]))
#     rd_p = scatter(x=df.kdam, y=df.rd, name="Rd", yaxis="y2", line=attr(width=width,color=colors[8]))
#     rt_p = scatter(x=df.kdam, y=df.rt, name="Rt", yaxis="y2", line=attr(width=width,color=colors[9]))
#     return rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, rh_p, rd_p, rt_p
# end


# function split_curves(df, df_bf)
#     kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
#     kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
#     first=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
#     middle=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
#     last=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    
#     for i in df.kdam[1:kdam1]
#         push!(first.kdam, i)
#     end
#     for i in df.kdam[kdam1:kdam2]
#         push!(middle.kdam, i)
#     end
#     for i in df.kdam[kdam2:end]
#         push!(last.kdam, i)
#     end
#     for (col,col1) in zip(eachcol(df)[1:9],eachcol(first)[2:end])
#         for i in col[1:kdam1]
#             push!(col1, i)
#         end
#     end
#     for (col,col1) in zip(eachcol(df)[1:9],eachcol(middle)[2:end])
#         for i in col[kdam1:kdam2]
#             push!(col1, i)
#         end
#     end
#     for (col,col1) in zip(eachcol(df)[1:9],eachcol(last)[2:end])
#         for i in col[kdam2:end]
#             push!(col1, i)
#         end
#     end
#     return first, middle, last
# end


# function split_curves_inhib(df, df_bf)
#     kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
#     kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
#     first=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
#     middle=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
#     last=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
    
#     for i in df.kdam[1:kdam1]
#         push!(first.kdam, i)
#     end
#     for i in df.kdam[kdam1:kdam2]
#         push!(middle.kdam, i)
#     end
#     for i in df.kdam[kdam2:end]
#         push!(last.kdam, i)
#     end
#     for (col,col1) in zip(eachcol(df)[1:10],eachcol(first)[2:end])
#         for i in col[1:kdam1]
#             push!(col1, i)
#         end
#     end
#     for (col,col1) in zip(eachcol(df)[1:10],eachcol(middle)[2:end])
#         for i in col[kdam1:kdam2]
#             push!(col1, i)
#         end
#     end
#     for (col,col1) in zip(eachcol(df)[1:10],eachcol(last)[2:end])
#         for i in col[kdam2:end]
#             push!(col1, i)
#         end
#     end
#     return first, middle, last
# end



# function plot_all_curves_bistable(br, colors2, colors_r, title)
#     df = create_br_df(br)
#     bf = bf_point_df(br)
#     bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt = bf_scatter(bf, "darkblue")
#     first_r, middle_r, last_r = dashed_lines_ribosomes(df, bf, colors_r, "")
#     first, middle, last = dashed_lines_species(df, bf, colors2, "")
#     return plot([first[1],middle[1],last[1], bfp_rma,
#     first[2],middle[2],last[2], bfp_rtca,
#     first[3],middle[3],last[3], bfp_rmb,
#     first[4],middle[4],last[4], bfp_rtcb,
#     first[6],middle[6],last[6], bfp_rtcr,
#     first_r[1],middle_r[1],last_r[1], bfp_rh,
#     first_r[2],middle_r[2],last_r[2], bfp_rd,
#     first_r[3],middle_r[3],last_r[3], bfp_rt],
#     Layout(yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
#     yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)", title="$title"))
# end

# function get_all_curves_for_bistab_plotting(br, colors2, colors_r)
#     df = create_br_df(br)
#     bf = bf_point_df(br)
#     bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt = bf_scatter(bf, "darkblue")
#     first_r, middle_r, last_r = dashed_lines_ribosomes(df, bf, colors_r, "")
#     first, middle, last = dashed_lines_species(df, bf, colors2, "")
#     return bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt, first_r, middle_r, last_r, first, middle, last
# end

# function plot_species_separately_ribosomes(br, colors_r, type)
#     df = create_br_df(br)
#     bf = bf_point_df(br)
#     first, middle, last = dashed_lines_ribosomes(df, bf, colors_r, "$type")
#     return first, middle, last
# end

# function plot_species_separately(br, colors_r, type)
#     df = create_br_df(br)
#     bf = bf_point_df(br)
#     first, middle, last = dashed_lines_species(df, bf, colors_r, "$type")
#     return first, middle, last
# end





