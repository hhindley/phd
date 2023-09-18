using BifurcationKit
const BK = BifurcationKit


# sup norm
norminf(x) = norm(x, Inf)

# vector field

function rtc_mod!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)

    # translation
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    
    dz
    
end


# parameter values
params = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., km_a = 20., km_b = 16., g_max = 2.0923, kdeg = 0.001, 
kdam =  0.01,
ω_ab = 2., ω_r = 0.0089, atp = 3000., kin = 0.022222, lam = 0.04)

params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014) 		


# initial condition
# using DifferentialEquations, DataFrames, StaticArrays
# include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
# tspan=(0,1e9)
# params1 = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# solu = sol(rtc_model, initial, tspan, params1)
# ss_init = ss_init_vals(solu)
initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]
initss = [0.0008735009426379191,
1.7342351146181998e-5,
0.0008735009426379191,
1.4366947763258614e-5,
0.044729189117567736,
9.403498079216664e-5,
0.048110392167741795,
0.2609572290416992,
2.77192276995412]
initial1 = [0., 0., 1., 0., 0., 0., 11.29, 0., 0.]

rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)

function get_br(model, params, initial, kdam_max)
    # Bifurcation Problem
    prob = BifurcationProblem(model, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    # # continuation options
    # opts_br = ContinuationPar(pMin = 0., pMax = 1.,
    # # parameters to have a smooth result
    # ds = 0.001, dsmax = 0.05,)
    opts_br = ContinuationPar(pMin = 0., pMax = kdam_max, ds = 0.001, dsmax = 0.01, 
    # options to detect bifurcations
    detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, 
    # number of eigenvalues
    nev = 2, 
    # maximum number of continuation steps
    maxSteps = 1000,)
    # continuation of equilibria
    # br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)

    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    return br
end



function rtc_mod_t!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = z


    # dilution by growth and degradation 
    dil(species) = lam(t)*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp(t)/(Km_init+atp(t))
    tscr_el_a = ω_ab*atp(t)/(θtscr+atp(t))
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp(t)/(θtscr+atp(t))
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp(t)/(θtscr+atp(t))

    # translation
    tlr_el = g_max*atp(t)/(θtlr+atp(t))
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp(t)*rtca)/(atp(t)+(km_a*rd)) 
    rtcb1 = (atp(t)*rtcb)/(atp(t)+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin(t)*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb)
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    # dz
    
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
    df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],kdam=[]);
    # for i in eachcol(df)
    #     println(i)
    # end
    for i in range(1,length(br.sol))
        for (s,d) in zip(range(1,9), eachcol(df))
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
    # for i in eachcol(df)
    #     println(i)
    # end
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

function dashed_lines_species(df, df_bf, colors)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=[]
    middle=[]
    last=[]
    names=["RtcBA mRNA","RtcA","RtcB mRNA","RtcB","RtcR mRNA","RtcR"]
    for (col,i) in zip(eachcol(df)[1:6],range(1,9))
        push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=names[i], line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"))
        push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"))
        push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"))
    end
    return first, middle, last
end
function dashed_lines_ribosomes(df, df_bf, colors)
    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=[]
    middle=[]
    last=[]
    names=["Healthy ribosomes","Damaged ribosomes","Tagged ribosomes"]
    for (col,i) in zip(eachcol(df)[7:9],range(1,9))
        push!(first,scatter(x=df.kdam[1:kdam1], y=col[1:kdam1], name=names[i], yaxis="y2", line=attr(width=3, color=colors[i]), legendgroup="$(names[i])"))
        push!(middle,scatter(x=df.kdam[kdam1:kdam2], y=col[kdam1:kdam2], name="", yaxis="y2", line=attr(width=3,dash="dash", color=colors[i]),showlegend=false, legendgroup="$(names[i])"))
        push!(last,scatter(x=df.kdam[kdam2:end], y=col[kdam2:end], name="", yaxis="y2", line=attr(width=3, color=colors[i]),showlegend=false, legendgroup="$(names[i])"))
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
    bf_rt = scatter(x=df_bf.kdam, y=df_bf.rt, mode="markers", yaxis="y2", name="Bifurcation point", line=attr(color=color),showlegend=true, legendgroup="Tagged ribosomes")
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

function plot_rtcb_bf(bf, df, kdam1, kdam2)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtcb[1:kdam1], name="RtcB", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtcb[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtcb[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf.rtcb, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    return rtcb1, rtcb2, rtcb3, bf_rtcb
end

function plot_rtca_bf(bf, df, kdam1, kdam2)
    rtcb1 = scatter(x=df.kdam[1:kdam1], y=df.rtca[1:kdam1], name="RtcA", line=attr(width=3, color=:green), showlegend=false, legendgroup="1")#, fill="tozeroy")
    rtcb2 = scatter(x=df.kdam[kdam1:kdam2], y=df.rtca[kdam1:kdam2], name="", line=attr(width=3,dash="dash", color=:black),showlegend=false, legendgroup="1")
    rtcb3 = scatter(x=df.kdam[kdam2:end], y=df.rtca[kdam2:end], name="", line=attr(width=3, color=:red),showlegend=false, legendgroup="1")
    bf_rtcb = scatter(x=bf.kdam, y=bf.rtca, mode="markers", name="Bifurcation point", line=attr(color=:black),showlegend=false, legendgroup="1")
    return rtcb1, rtcb2, rtcb3, bf_rtcb
end