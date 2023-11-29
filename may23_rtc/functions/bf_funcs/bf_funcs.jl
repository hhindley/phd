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
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = kdam_max, ds = 0.001, 
        dsmax = 0.15, # 0.001, 0.15
        # options to detect bifurcations
        detectBifurcation = 3, nInversion = 4, maxBisectionSteps = 10, #3,2,10
        # number of eigenvalues
        nev = 2,  #4 
        # maximum number of continuation steps
        maxSteps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution
    elseif length(params) == 28
        prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = kdam_max, ds = 0.001, 
        dsmax = 0.05, # 0.001, 0.15
        # options to detect bifurcations
        detectBifurcation = 3, nInversion = 4, maxBisectionSteps = 10, #3,2,10
        # number of eigenvalues
        nev = 2,  #4 
        # maximum number of continuation steps
        maxSteps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
        # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
        # only using parameters that make a difference to solution

    else
        prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        # continuation options
        opts_br = ContinuationPar(pMin = 0., pMax = kdam_max, ds = 0.001, dsmax = 0.01,
        # options to detect bifurcations
        detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25,
        # number of eigenvalues
        nev = 2, 
        # maximum number of continuation steps
        maxSteps = 1000,) # only using parameters that make a difference to solution

    end
    # continuation of equilibria
    br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    # br = continuation(prob, PALC(θ=0.75), opts_br; plot = false, bothside=true, normC = norminf)
    return br
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


function split_curves_inhib(df, df_bf)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
    middle=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
    last=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],rtcb_i=[])
    
    for i in df.kdam[1:kdam1]
        push!(first.kdam, i)
    end
    for i in df.kdam[kdam1:kdam2]
        push!(middle.kdam, i)
    end
    for i in df.kdam[kdam2:end]
        push!(last.kdam, i)
    end
    for (col,col1) in zip(eachcol(df)[1:10],eachcol(first)[2:end])
        for i in col[1:kdam1]
            push!(col1, i)
        end
    end
    for (col,col1) in zip(eachcol(df)[1:10],eachcol(middle)[2:end])
        for i in col[kdam1:kdam2]
            push!(col1, i)
        end
    end
    for (col,col1) in zip(eachcol(df)[1:10],eachcol(last)[2:end])
        for i in col[kdam2:end]
            push!(col1, i)
        end
    end
    return first, middle, last
end

function dashed_lines_species(df, df_bf, colors, type)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=[]
    middle=[]
    last=[]
    names=["RtcBA mRNA","RtcA","RtcB mRNA","RtcB","RtcR mRNA","RtcR"]
    if type == ""
        for (col,i) in zip(eachcol(df)[1:6],range(1,9))
            push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=(names[i]*" $type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
            push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
            push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
        end
    else
        for (col,i) in zip(eachcol(df)[1:6],range(1,9))
            push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=("$type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
            push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
            push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
        end
    end
    return first, middle, last
end
function dashed_lines_ribosomes(df, df_bf, colors, type)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=[]
    middle=[]
    last=[]
    names=["Healthy ribosomes","Damaged ribosomes","Tagged ribosomes"]
    if type == ""
        for (col,i) in zip(eachcol(df)[7:9],range(1,9))
            push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=(names[i]*" $type"), yaxis="y2", line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
            push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", yaxis="y2", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
            push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", yaxis="y2", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
        end
    else
        for (col,i) in zip(eachcol(df)[7:9],range(1,9))
            push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=("$type"), line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"*"$type"))
            push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
            push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"*"$type"))
        end
    end
    return first, middle, last
end

function bf_scatter(df_bf, color)
    bf_rma = scatter(x=df_bf.kdam, y=df_bf.rm_a, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcBA mRNA")
    bf_rtca = scatter(x=df_bf.kdam, y=df_bf.rtca, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcA")
    bf_rmb = scatter(x=df_bf.kdam, y=df_bf.rm_b, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcB mRNA")
    bf_rtcb = scatter(x=df_bf.kdam, y=df_bf.rtcb, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcB")
    bf_rmr = scatter(x=df_bf.kdam, y=df_bf.rm_r, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcR mRNA")
    bf_rtcr = scatter(x=df_bf.kdam, y=df_bf.rtcr, mode="markers", name="", line=attr(color=color),showlegend=false, legendgroup="RtcR")
    bf_rh = scatter(x=df_bf.kdam, y=df_bf.rh, mode="markers", yaxis="y2", name="", line=attr(color=color),showlegend=false, legendgroup="Healthy ribosomes")
    bf_rd = scatter(x=df_bf.kdam, y=df_bf.rd, mode="markers", yaxis="y2", name="", line=attr(color=color),showlegend=false, legendgroup="Damaged ribosomes")
    bf_rt = scatter(x=df_bf.kdam, y=df_bf.rt, mode="markers", yaxis="y2", name="Bifurcation point", line=attr(color=color),showlegend=false, legendgroup="Tagged ribosomes")
    return bf_rma, bf_rtca, bf_rmb, bf_rtcb, bf_rmr, bf_rtcr, bf_rh, bf_rd, bf_rt
end


function full_lines(df,width, colors)
    rma_p = scatter(x=df.kdam, y=df.rm_a, name="RtcBA mRNA", line=attr(width=width,color=colors[1]))
    rtca_p = scatter(x=df.kdam, y=df.rtca, name="RtcA", line=attr(width=width,color=colors[2]))
    rmb_p = scatter(x=df.kdam, y=df.rm_b, name="rm_b", line=attr(width=width,color=colors[3]))
    rtcb_p = scatter(x=df.kdam, y=df.rtcb, name="RtcB", line=attr(width=width,color=colors[4]))
    rmr_p = scatter(x=df.kdam, y=df.rm_r, name="rm_r", line=attr(width=width,color=colors[5]))
    rtcr_p = scatter(x=df.kdam, y=df.rtcr, name="RtcR", line=attr(width=width,color=colors[6]))
    rh_p = scatter(x=df.kdam, y=df.rh, name="Rh", yaxis="y2", line=attr(width=width,color=colors[7]))
    rd_p = scatter(x=df.kdam, y=df.rd, name="Rd", yaxis="y2", line=attr(width=width,color=colors[8]))
    rt_p = scatter(x=df.kdam, y=df.rt, name="Rt", yaxis="y2", line=attr(width=width,color=colors[9]))
    return rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, rh_p, rd_p, rt_p
end


function different_levels_inhibition(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)
    params_for_ssval_setup_inhib = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib)
    br = get_br(rtc_inhib_mod, params_for_ssval_setup_inhib, initial_i, 3.)
    bf = bf_point_df_inhib(br)
    df = create_br_df_inhib(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    return bf, df, kdam1, kdam2
end

function plot_rtc_bf(bf, df, kdam1, kdam2, specie)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df[!,specie][1:kdam1], name="$specie", line=attr(width=5, color="#005356ff"), showlegend=true)#, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df[!,specie][kdam1:kdam2], name="", line=attr(width=5,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df[!,specie][kdam2:end], name="", line=attr(width=5, color="#f04e53ff"),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf[!,specie], mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    return rtcb1, rtcb2, rtcb3, bf_rtcb
end


function plot_all_curves_bistable(br, colors2, colors_r, title)
    df = create_br_df(br)
    bf = bf_point_df(br)
    bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt = bf_scatter(bf, "darkblue")
    first_r, middle_r, last_r = dashed_lines_ribosomes(df, bf, colors_r, "")
    first, middle, last = dashed_lines_species(df, bf, colors2, "")
    return plot([first[1],middle[1],last[1], bfp_rma,
    first[2],middle[2],last[2], bfp_rtca,
    first[3],middle[3],last[3], bfp_rmb,
    first[4],middle[4],last[4], bfp_rtcb,
    first[6],middle[6],last[6], bfp_rtcr,
    first_r[1],middle_r[1],last_r[1], bfp_rh,
    first_r[2],middle_r[2],last_r[2], bfp_rd,
    first_r[3],middle_r[3],last_r[3], bfp_rt],
    Layout(yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
    yaxis_title="Proteins and mRNAs (μM)", yaxis2_title="Ribosomal species (μM)", title="$title"))
end

function get_all_curves_for_bistab_plotting(br, colors2, colors_r)
    df = create_br_df(br)
    bf = bf_point_df(br)
    bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt = bf_scatter(bf, "darkblue")
    first_r, middle_r, last_r = dashed_lines_ribosomes(df, bf, colors_r, "")
    first, middle, last = dashed_lines_species(df, bf, colors2, "")
    return bfp_rma, bfp_rtca, bfp_rmb, bfp_rtcb, bfp_rmr, bfp_rtcr, bfp_rh, bfp_rd, bfp_rt, first_r, middle_r, last_r, first, middle, last
end

function plot_species_separately_ribosomes(br, colors_r, type)
    df = create_br_df(br)
    bf = bf_point_df(br)
    first, middle, last = dashed_lines_ribosomes(df, bf, colors_r, "$type")
    return first, middle, last
end

function plot_species_separately(br, colors_r, type)
    df = create_br_df(br)
    bf = bf_point_df(br)
    first, middle, last = dashed_lines_species(df, bf, colors_r, "$type")
    return first, middle, last
end