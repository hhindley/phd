using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra, DataFrames
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

				


# initial condition
# using DifferentialEquations, DataFrames, StaticArrays
# include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
# tspan=(0,1e9)
# params1 = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# solu = sol(rtc_model, initial, tspan, params1)
# ss_init = ss_init_vals(solu)
initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)

# Bifurcation Problem
prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.0), (@lens _.kdam);
recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)

# # continuation options
# opts_br = ContinuationPar(pMin = 0., pMax = 1.,
# # parameters to have a smooth result
# ds = 0.001, dsmax = 0.05,)

opts_br = ContinuationPar(pMin = 0., pMax = 1.0,  dsmax = 0.01, ds = 0.001,
	# options to detect bifurcations
	detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25,
	# number of eigenvalues
	nev = 2,
	# maximum number of continuation steps
	maxSteps = 1000,)


# continuation of equilibria
br = continuation(prob, PALC(), opts_br;
plot = false, bothside=true, normC = norminf)




plot(br, plotfold=true, plotspecialpoints=true, markersize=3, legend=:best)

plot(br, vars = (:param, :rh))
plot(br, vars = (:param, :rt))
plot(br, vars = (:param, :rd))
plot(br, vars = (:param, :rm_r))


# cont2 = continuation(br, 1, (@lens _.lam), ContinuationPar(opts_br, pMax = 0.04, dsmax=0.05);
# normC = norminf, detectCodim2Bifurcation=2, updateMinAugEveryStep=1, bothside=true)

wab_range = 10 .^range(-5,stop=0,length=10)
wr_range = 10 .^(range(-7,stop=0,length=10))

atp_range = range(500,stop=5000,length=10)
kin_range = range(0,stop=0.2,length=10)
lam_range = range(0.001,stop=0.04,length=10)

wr_range = (range(0.0001,stop=0.001,length=10)) # use this range with the above 3 ranges and get ~3000 combos that give bistability
wab_range = range(0.01, stop=1., length=10)

function run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)
    bistab1_df=DataFrame(wab = Float64[], wr = Float64[], kin = Float64[], lam = Float64[], bs_type = Symbol[])
    bistab1 = []
    for i in wab_range
        params = merge(params, (ω_ab=i,))
        for j in wr_range
            params = merge(params, (ω_r=j,))
            # for k in atp_range 
            #     params = merge(params, (atp=k,))
                for l in kin_range
                    params = merge(params, (kin=l,))

                    for m in lam_range
                        params = merge(params, (lam=m,))
                        # println("lam = $m")
                        # println("kin = $l, lam = $m")
                        # println("atp = $k, kin = $l, lam = $m")
                        # println("wr = $j, atp = $k, kin = $l, lam = $m")
                        # println("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m")

                        prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
                        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
                        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
                        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
                        for b in range(1,length(br.specialpoint))
                            if br.specialpoint[b].type == :endpoint
                                Nothing
                            else
                                # push!(bistab, ("lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab1, ("wab = $i, wr = $j, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab1, ("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))

                                push!(bistab1_df, (i, j, l, m, br.specialpoint[b].type))
                            end
                        end                     
                    end
                end
            # end
        end
    end
    return bistab1_df
end
bistab1_df = @time run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)


bistab1_df = bistab1_df[1:2:end,:]

CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/param_sweep_bs.csv", bistab1_df)

bistab # big 7000 element vector for bistab values using these ranges - atp_range = range(500,stop=5000,length=20), kin_range = range(0,stop=0.2,length=20), lam_range = range(0.001,stop=0.04,length=20), wr_range = (range(0.00001,stop=0.01,length=20))

bistab_5 = bistab[1:2:end]

bistab_5[1][1]


param_combos_bs = bistab_5[1000]
bistab_5[1000:1200]

[println(i) for i in bistab[1:2:end]]
lam_bi = [0.005578947368421052, 0.006, 0.0063236842105263156, 0.006394736842105263, 0.0074444444444444445, 0.007902105263157894, 0.00796842105263158, 0.009480526315789474, 0.009542105263157895, 0.01, 0.010105105105105105]

function bifu1()
    # params = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    #     θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    #     krep = 137., ktag = 9780., atp = 500., km_a = 20., km_b = 16., g_max = 2.0923, 
    #     kdeg = 0.001, kin = 0.010526315789473684, ω_ab = 4., ω_r = 2e-7, 
    #     kdam =  0.01, lam = 0.01)

    params = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 4000., km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.054, ω_ab = 4., ω_r = 1e-6, 
    kdam =  0.01, lam = 0.01)

    initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

    rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)

    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1., ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)

    scene = plot(br, plotfold=true, plotspecialpoints=true, markersize=3, legend=:best)
end
# plot(br, vars = (:param, :rh))
bifu1()



# params_wab = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
#     θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
#     krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
#     kdeg = 0.001, kin = 0.010526315789473684, ω_ab = 0.00022342342342342343, ω_r = 1e-5, 
#     kdam =  0.01, lam = 0.009210526315789473)


# wab_range = range(1, stop=4, length=100)

# df=DataFrame(wab = Float64[], kdam1 = Float64[], kdam2 = Float64[])

# for i in wab_range
#     # @show parameter
#     params = merge(params_wab, (ω_ab=i,))
#     println("wab = $i")
#     prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
#     recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
#     opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
#     br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
#     for b in range(1,length(br.specialpoint))
#         if br.specialpoint[b].type == :endpoint
#             Nothing
#             # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
#         else
#             push!(df, (i, br.specialpoint[2].param, br.specialpoint[3].param))
#         end
#     end
# end
# p = scatter(df.wab, df.kdam1, xlabel="wab", ylabel="kdam")
# p = scatter!(df.wab, df.kdam2)



params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014)


atp_range = range(500,stop=5000,length=100)

df_atp=DataFrame(atp = Float64[], kdam1 = Float64[], kdam2 = Float64[])

for i in atp_range
    # @show parameter
    params = merge(params_atp, (atp=i,))
    println("atp = $i")
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            Nothing
            # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
        else
            push!(df_atp, (i, br.specialpoint[2].param, br.specialpoint[3].param))
        end
    end
end
df_atp = df_atp[1:2:end, :]
p = scatter(df_atp.atp, df_atp.kdam1, xlabel="ATP", ylabel="kdam")
p = scatter!(df_atp.atp, df_atp.kdam2)



kin_range = range(0,stop=0.2,length=500)

df_kin=DataFrame(kin = Float64[], kdam1 = Float64[], kdam2 = Float64[])

for i in kin_range
    # @show parameter
    params = merge(params_atp, (kin=i,))
    println("kin = $i")
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            Nothing
            # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
        else
            push!(df_kin, (i, br.specialpoint[2].param, br.specialpoint[3].param))
        end
    end
end
df_kin = df_kin[1:2:end, :]
p = scatter(df_kin.kin, df_kin.kdam1, xlabel="kin", ylabel="kdam")
p = scatter!(df_kin.kin, df_kin.kdam2)



lam_range = range(0.001,stop=0.04,length=500)

df_lam=DataFrame(lam = Float64[], kdam1 = Float64[], kdam2 = Float64[])

for i in lam_range
    # @show parameter
    params = merge(params_atp, (lam=i,))
    println("lam = $i")
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            Nothing
            # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
        else
            push!(df_lam, (i, br.specialpoint[2].param, br.specialpoint[3].param))
        end
    end
end
df_lam = df_lam[1:2:end, :]
p = scatter(df_lam.lam, df_lam.kdam1, xlabel="lam", ylabel="kdam")
p = scatter!(df_lam.lam, df_lam.kdam2)




wab_range = range(0.01, stop=4, length=100)

df_wab=DataFrame(wab = Float64[], kdam1 = Float64[], kdam2 = Float64[])

for i in wab_range
    params = merge(params_atp, (ω_ab=i,))
    println("wab = $i")
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            Nothing
        else
            push!(df_wab, (i, br.specialpoint[2].param, br.specialpoint[3].param))
        end
    end
end
df_wab = df_wab[1:2:end, :]
p = scatter(df_wab.wab, df_wab.kdam1, xlabel="wab", ylabel="kdam")
p = scatter!(df_wab.wab, df_wab.kdam2)



wr_range = (range(0.00001,stop=0.001,length=200))
df_wr=DataFrame(wr = Float64[], kdam1 = Float64[], kdam2 = Float64[])

for i in wr_range
    # @show parameter
    params = merge(params_atp, (ω_r=i,))
    println("wr = $i")
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            Nothing
            # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
        else
            push!(df_wr, (i, br.specialpoint[2].param, br.specialpoint[3].param))
        end
    end
end
df_wr = df_wr[1:2:end, :]
p = scatter(df_wr.wr, df_wr.kdam1, xlabel="wr", ylabel="kdam")
p = scatter!(df_wr.wr, df_wr.kdam2)










atp_range = range(500,stop=5000,length=100)
kin_range = range(0,stop=0.2,length=100)
lam_range = range(0.001,stop=0.04,length=100)
wr_range = (range(0.00001,stop=0.001,length=100))
wab_range = range(0.01, stop=4, length=100)

# atp vs wab 
df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in atp_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (atp=i,))
    # println("atp = $i")
    for j in wab_range
        params_atp = merge(params_atp, (ω_ab=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df[df.bs .== :bp, :].atp, df[df.bs .== :bp, :].wab, xlabel="ATP", ylabel="wab", label="bistable region")
scatter!(df[df.bs .== :endpoint, :].atp, df[df.bs .== :endpoint, :].wab, label="no bistability")

plot(df[df.bs .== :bp, :].atp, df[df.bs .== :bp, :].wab)
scatter!(df[df.bs .== :bp, :].atp, df[df.bs .== :bp, :].wab)

function plot_bistable_region(df, param_range, param1, param2 )
    bsp = df[df.bs .== :bp, :]
    max_wab = []
    min_wab = []
    for i in param_range
        push!(max_wab, maximum(bsp[bsp[:,1] .== i, :][:,2]))
        push!(min_wab, minimum(bsp[bsp[:,1] .== i, :][:,2]))
    end

    return display(plot(param_range, [min_wab, max_wab], fillrange = (max_wab.-0), fillalpha = 0.35, linecolor=:lightblue, legend=false, xlabel="$param1", ylabel="$param2"))
end

plot_bistable_region(df, atp_range, "atp", "wab")



# atp vs wr 
df1 = DataFrame(atp = Float64[], wr = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in atp_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (atp=i,))
    # println("atp = $i")
    for j in wr_range
        params_atp = merge(params_atp, (ω_r=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df1, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df1, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df1[df1.bs .== :bp, :].atp, df1[df1.bs .== :bp, :].wr, xlabel="ATP", ylabel="wr", label="bistable region")
scatter!(df1[df1.bs .== :endpoint, :].atp, df1[df1.bs .== :endpoint, :].wr, label="no bistability")
plot_bistable_region(df1, atp_range, "atp", "wr")


# atp vs kin
df2 = DataFrame(atp = Float64[], kin = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in atp_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (atp=i,))
    # println("atp = $i")
    for j in kin_range
        params_atp = merge(params_atp, (kin=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df2, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df2, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df2[df2.bs .== :bp, :].atp, df2[df2.bs .== :bp, :].kin, xlabel="ATP", ylabel="kin", label="bistable region")
scatter!(df2[df2.bs .== :endpoint, :].atp, df2[df2.bs .== :endpoint, :].kin, label="no bistability")
plot_bistable_region(df2, atp_range, "atp", "kin")


# atp vs lam
df3 = DataFrame(atp = Float64[], lam = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in atp_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (atp=i,))
    # println("atp = $i")
    for j in lam_range
        params_atp = merge(params_atp, (lam=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df3, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df3, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df3[df3.bs .== :bp, :].atp, df3[df3.bs .== :bp, :].lam, xlabel="ATP", ylabel="lam", label="bistable region")
scatter!(df3[df3.bs .== :endpoint, :].atp, df3[df3.bs .== :endpoint, :].lam, label="no bistability")
plot_bistable_region(df3, atp_range, "atp", "lam")


# wr vs wab
df4 = DataFrame(wr = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in wr_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (ω_r=i,))
    # println("atp = $i")
    for j in wab_range
        params_atp = merge(params_atp, (ω_ab=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df4, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df4, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df4[df4.bs .== :bp, :].wr, df4[df4.bs .== :bp, :].wab, xlabel="wr", ylabel="wab", label="bistable region")
scatter!(df4[df4.bs .== :endpoint, :].wr, df4[df4.bs .== :endpoint, :].wab, label="no bistability")
plot_bistable_region(df4, wr_range, "wr", "wab")


# wr vs kin
df5 = DataFrame(wr = Float64[], kin = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in wr_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (ω_r=i,))
    # println("atp = $i")
    for j in kin_range
        params_atp = merge(params_atp, (kin=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df5, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df5, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df5[df5.bs .== :bp, :].wr, df5[df5.bs .== :bp, :].kin, xlabel="wr", ylabel="kin", label="bistable region")
scatter!(df5[df5.bs .== :endpoint, :].wr, df5[df5.bs .== :endpoint, :].kin, label="no bistability")
plot_bistable_region(df5, wr_range, "wr", "kin")


# wr vs lam 
df6 = DataFrame(wr = Float64[], lam = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in wr_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (ω_r=i,))
    # println("atp = $i")
    for j in lam_range
        params_atp = merge(params_atp, (lam=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df6, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df6, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df6[df6.bs .== :bp, :].wr, df6[df6.bs .== :bp, :].lam, xlabel="wr", ylabel="lam", label="bistable region")
scatter!(df6[df6.bs .== :endpoint, :].wr, df6[df6.bs .== :endpoint, :].lam, label="no bistability")
plot_bistable_region(df6, wr_range, "wr", "lam")


# wab vs kin
df7 = DataFrame(wab = Float64[], kin = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in wab_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (ω_ab=i,))
    # println("atp = $i")
    for j in kin_range
        params_atp = merge(params_atp, (kin=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df7, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df7, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df7[df7.bs .== :bp, :].wab, df7[df7.bs .== :bp, :].kin, xlabel="wab", ylabel="kin", label="bistable region")
scatter!(df7[df7.bs .== :endpoint, :].wab, df7[df7.bs .== :endpoint, :].kin, label="no bistability")
plot_bistable_region(df7, wab_range, "wab", "kin")


# wab vs lam 
df8 = DataFrame(wab = Float64[], lam = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in wab_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (ω_ab=i,))
    # println("atp = $i")
    for j in lam_range
        params_atp = merge(params_atp, (lam=j,))
        prob = BifurcationProblem(rtc_mod, intial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df8, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df8, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df8[df8.bs .== :bp, :].wab, df8[df8.bs .== :bp, :].lam, xlabel="wab", ylabel="lam", label="bistable region")
scatter!(df8[df8.bs .== :endpoint, :].wab, df8[df8.bs .== :endpoint, :].lam, label="no bistability")
plot_bistable_region(df8, wab_range, "wab", "lam")


# kin vs lam 
df9 = DataFrame(kin = Float64[], lam = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
for i in kin_range
    params_atp = deepcopy(params1)
    params_atp = merge(params_atp, (kin=i,))
    # println("atp = $i")
    for j in lam_range
        params_atp = merge(params_atp, (lam=j,))
        prob = BifurcationProblem(rtc_mod, initial, setproperties(params_atp; kdam=0.), (@lens _.kdam);
        recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
        br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
        if length(br.specialpoint) == 2
            push!(df9, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
        else
            push!(df9, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
        end
    end    
end
scatter(df9[df9.bs .== :bp, :].kin, df9[df9.bs .== :bp, :].lam, xlabel="wab", ylabel="lam", label="bistable region")
scatter!(df9[df9.bs .== :endpoint, :].kin, df9[df9.bs .== :endpoint, :].lam, label="no bistability")
plot_bistable_region(df9, kin_range, "kin", "lam")










kin_range = range(0,stop=0.2,length=10)
atp_range = range(500,stop=5000,length=10)
lam_range = range(0.001,stop=0.04,length=10)
wr_range = (range(0.00001,stop=0.001,length=10))
wab_range = range(0.01, stop=4, length=10)

p1 = plot()
for i in atp_range
    params = merge(params_atp, (atp=i,))
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    display(plot!(p1, br, legend=false))
end


p1 = plot()
for i in kin_range
    params = merge(params_atp, (kin=i,))
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    display(plot!(p1, br, legend=false))
end


p1 = plot()
for i in lam_range
    params = merge(params_atp, (lam=i,))
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    display(plot!(p1, br, legend=false))
end

p1 = plot()
for i in wab_range
    params = merge(params_atp, (ω_ab=i,))
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    display(plot!(p1, br, legend=false))
end

p1 = plot()
for i in wr_range
    params = merge(params_atp, (ω_r=i,))
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    display(plot!(p1, br, legend=false))
end