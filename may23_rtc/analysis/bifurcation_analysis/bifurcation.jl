using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames
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

rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)

function get_br(params, initial)
    # Bifurcation Problem
    prob = BifurcationProblem(rtc_mod, initial, setproperties(params; kdam=0.), (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    # # continuation options
    # opts_br = ContinuationPar(pMin = 0., pMax = 1.,
    # # parameters to have a smooth result
    # ds = 0.001, dsmax = 0.05,)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, ds = 0.001, dsmax = 0.01, 
    # options to detect bifurcations
    detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, 
    # number of eigenvalues
    nev = 2, 
    # maximum number of continuation steps
    maxSteps = 1000,)
    # continuation of equilibria
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    return br
end

br = get_br(params1, initial)


plot(br, vars=(:param, :rh), putspecialptlegend=false, label="rfw")

p_rh = plot(br, vars = (:param, :rh), label="Healthy ribosomes", c=:palevioletred, putspecialptlegend = false)
p_rh1 = plot!(twinx(), br, vars = (:param, :rtca), c=:grey60, label="RtcBA")

savefig(p_rm_a, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_rma.svg")
savefig(p_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_rh.svg")
savefig(p_rh1, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_both.svg")



# wab_range = 10 .^range(-5,stop=0,length=10)
# wr_range = 10 .^(range(-7,stop=0,length=10))

# atp_range = range(500,stop=5000,length=10)
# kin_range = range(0,stop=0.2,length=10)
# lam_range = range(0.001,stop=0.04,length=10)

# wr_range = (range(0.0001,stop=0.001,length=10)) # use this range with the above 3 ranges and get ~3000 combos that give bistability
# wab_range = range(0.01, stop=1., length=10)

# function run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)
#     bistab1_df=DataFrame(wab = Float64[], wr = Float64[], kin = Float64[], lam = Float64[], bs_type = Symbol[])
#     bistab1 = []
#     params = deepcopy(params1)
#     for i in wab_range
#         params = merge(params, (ω_ab=i,))
#         for j in wr_range
#             params = merge(params, (ω_r=j,))
#             # for k in atp_range 
#             #     params = merge(params, (atp=k,))
#                 for l in kin_range
#                     params = merge(params, (kin=l,))

#                     for m in lam_range
#                         params = merge(params, (lam=m,))
#                         # println("lam = $m")
#                         # println("kin = $l, lam = $m")
#                         # println("atp = $k, kin = $l, lam = $m")
#                         # println("wr = $j, atp = $k, kin = $l, lam = $m")
#                         # println("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m")

#                         br = get_br(params, initial)
#                         for b in range(1,length(br.specialpoint))
#                             if br.specialpoint[b].type == :endpoint
#                                 Nothing
#                             else
#                                 # push!(bistab, ("lam = $m", br.specialpoint[b].type))
#                                 # push!(bistab, ("kin = $l, lam = $m", br.specialpoint[b].type))
#                                 # push!(bistab, ("atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
#                                 # push!(bistab, ("wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
#                                 # push!(bistab1, ("wab = $i, wr = $j, kin = $l, lam = $m", br.specialpoint[b].type))
#                                 # push!(bistab1, ("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))

#                                 push!(bistab1_df, (i, j, l, m, br.specialpoint[b].type))
#                             end
#                         end                     
#                     end
#                 end
#             # end
#         end
#     end
#     return bistab1_df
# end
# bistab1_df = @time run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)


# bistab1_df = bistab1_df[1:2:end,:]

# CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/param_sweep_bs.csv", bistab1_df)


function plot_bf_against_kdam(param_range, param)
    df=DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[])
    paramscopy = deepcopy(params1)
    for i in param_range
        params = merge(paramscopy, (param=>i,))
        # println("$param = $i")
        br = get_br(params, initial)
        for b in range(1,length(br.specialpoint))
            if br.specialpoint[b].type == :endpoint
                Nothing
                # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
            else
                push!(df, (i, br.specialpoint[2].param, br.specialpoint[3].param))
            end
        end
    end
    df = df[1:2:end, :]
    p = plot(df.param, df.kdam1, xlabel="$param", ylabel="kdam", ylims=(0,1), label="bf1")
    p = plot!(df.param, df.kdam2, label="bf2")

end

atp_range = range(500,stop=5000,length=500)
kin_range = range(0,stop=0.2,length=500)
lam_range = range(0.001,stop=0.04,length=500)
wab_range = range(0.01, stop=4, length=500)
wr_range = (range(0.00001,stop=0.001,length=200))

p_atp = plot_bf_against_kdam(atp_range, :atp)
p_kin = plot_bf_against_kdam(kin_range, :kin)
p_lam = plot_bf_against_kdam(lam_range, :lam)
p_wab = plot_bf_against_kdam(wab_range, :ω_ab)
p_wr = plot_bf_against_kdam(wr_range, :ω_r)


l = @layout [a b; c d; e]
all_bf = plot(p_atp, p_kin, p_lam, p_wab, p_wr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Bifurcation points for kdam vs. parameter")
savefig(all_bf, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_bf.svg")






atp_range = range(500,stop=5000,length=100)
kin_range = range(0.001,stop=0.2,length=100)
lam_range = range(0.001,stop=0.04,length=100)
wr_range = (range(0.00001,stop=0.001,length=100))
wab_range = range(0.01, stop=4, length=100)


function double_param_vary(param_range1, param1, param_range2, param2)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            br = get_br(params, initial)
            if length(br.specialpoint) == 2
                push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end    
    end
    return df
end
function plot_bistable_region(param_range1, param1, param_range2, param2)
    df = double_param_vary(param_range1, param1, param_range2, param2)
    bsp = df[df.bs .== :bp, :]
    max_wab = []
    min_wab = []
    for i in param_range1
        push!(max_wab, maximum(bsp[bsp[:,1] .== i, :][:,2]))
        push!(min_wab, minimum(bsp[bsp[:,1] .== i, :][:,2]))
    end
    p = plot(param_range1, fill(maximum(param_range2),length(param_range1)); fillrange=(minimum(param_range2)), fillcolor=:lightblue2, linecolor=:lightblue, label="", xlabel="$param1", ylabel="$param2")
    p = plot!(param_range1, max_wab; fillrange=(minimum(param_range2)), fillcolor=:cyan4, linecolor=:lightblue, label="bistable region", legend=false)
    p = plot!(param_range1, min_wab; fillrange=(minimum(param_range2)), fillcolor=:lightblue2, linecolor=:lightblue, label="")
    return p
end

# atp vs wab
atp_wab = plot_bistable_region(atp_range, :atp, wab_range, :ω_ab)
# atp vs wr 
atp_wr = plot_bistable_region(atp_range, :atp, wr_range, :ω_r)
# atp vs kin
atp_kin = plot_bistable_region(atp_range, :atp, kin_range, :kin)
# atp vs lam
atp_lam = plot_bistable_region(atp_range, :atp, lam_range, :lam)
# wr vs wab
wr_wab = plot_bistable_region(wr_range, :ω_r, wab_range, :ω_ab)
# wr vs kin
wr_kin = plot_bistable_region(wr_range, :ω_r, kin_range, :kin)
# wr vs lam 
wr_lam = plot_bistable_region(wr_range, :ω_r, lam_range, :lam)
# wab vs kin
wab_kin = plot_bistable_region(wab_range, :ω_ab, kin_range, :kin)
# wab vs lam 
wab_lam = plot_bistable_region(wab_range, :ω_ab, lam_range, :lam)
# kin vs lam - get instability here 
kin_lam = plot_bistable_region(kin_range, :kin, lam_range, :lam)

l = @layout [a b c; d e f; g h i j]
all_ranges = plot(kin_lam, wab_lam, wab_kin, wr_lam, wr_wab, wr_kin, atp_kin, atp_lam, atp_wab, atp_wr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Parameter spaces where bistability is present (dark area = bs)")

savefig(all_ranges, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_ranges.svg")
savefig(atp_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wab_atp_bigger_range.svg")

function plot_triple_param_vary(param_range1, param1, param_range2, param2, param_range3, param3)
    copyparams=deepcopy(params1)
    p = plot();
    for i in param_range3
        params = merge(copyparams, (param3=>i,))
        df = double_param_vary(param_range1, param1, param_range2, param2, params)
        bsp = df[df.bs .== :bp, :]
        max_wab = []
        min_wab = []
        for i in param_range1
            push!(max_wab, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_wab, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
        plot!(p, param_range1, max_wab, linecolor=:lightblue, label="$param3 = $i", xlabel="$param1", ylabel="$param2")
        plot!(p, param_range1, min_wab, linecolor=:lightblue, label="$param3 = $i")
    end
    return p 
end

plot_triple_param_vary(atp_range, :atp, kin_range, :kin, wab_range, :ω_ab)



kin_range = range(0,stop=0.2,length=10)
atp_range = range(500,stop=5000,length=10)
lam_range = range(0.001,stop=0.04,length=10)
wr_range = (range(0.00001,stop=0.001,length=10))
wab_range = range(0.01, stop=4, length=10)

function plot_different_bp_lines(param_range, params1, param, specie)
    p = plot();
    copyparams = deepcopy(params1)
    for i in param_range
        params = merge(copyparams, (param=>i,))
        br = get_br(params, initial)
        plot!(p, br, vars = (:param, specie), legend=false)#, label="$param = $i", putspecialptlegend=false)
    end
    return p
end



p1_rma = plot_different_bp_lines(atp_range, params1, :atp, :rm_a)
p2_rma = plot_different_bp_lines(kin_range, params1, :kin, :rm_a)
p3_rma = plot_different_bp_lines(lam_range, params1, :lam, :rm_a)
p4_rma = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_a)
p5_rma = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_a)

l = @layout [a b; c d; e]
all_rma = plot(p1_rma, p2_rma, p3_rma, p4_rma, p5_rma, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rma, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rma.svg")

p1_rmr = plot_different_bp_lines(atp_range, params1, :atp, :rm_r)
p2_rmr = plot_different_bp_lines(kin_range, params1, :kin, :rm_r)
p3_rmr = plot_different_bp_lines(lam_range, params1, :lam, :rm_r)
p4_rmr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rm_r)
p5_rmr = plot_different_bp_lines(wr_range, params1, :ω_r, :rm_r)

all_rmr = plot(p1_rmr, p2_rmr, p3_rmr, p4_rmr, p5_rmr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rmr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rmr.svg")

p1_rtca = plot_different_bp_lines(atp_range, params1, :atp, :rtca)
p2_rtca = plot_different_bp_lines(kin_range, params1, :kin, :rtca)
p3_rtca = plot_different_bp_lines(lam_range, params1, :lam, :rtca)
p4_rtca = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtca)
p5_rtca = plot_different_bp_lines(wr_range, params1, :ω_r, :rtca)

all_rtca = plot(p1_rtca, p2_rtca, p3_rtca, p4_rtca, p5_rtca, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rtca, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtca.svg")

p1_rtcr = plot_different_bp_lines(atp_range, params1, :atp, :rtcr)
p2_rtcr = plot_different_bp_lines(kin_range, params1, :kin, :rtcr)
p3_rtcr = plot_different_bp_lines(lam_range, params1, :lam, :rtcr)
p4_rtcr = plot_different_bp_lines(wab_range, params1, :ω_ab, :rtcr)
p5_rtcr = plot_different_bp_lines(wr_range, params1, :ω_r, :rtcr)

all_rtcr = plot(p1_rtcr, p2_rtcr, p3_rtcr, p4_rtcr, p5_rtcr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rtcr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rtcr.svg")

p1_rh = plot_different_bp_lines(atp_range, params1, :atp, :rh)
p2_rh = plot_different_bp_lines(kin_range, params1, :kin, :rh)
p3_rh = plot_different_bp_lines(lam_range, params1, :lam, :rh)
p4_rh = plot_different_bp_lines(wab_range, params1, :ω_ab, :rh)
p5_rh = plot_different_bp_lines(wr_range, params1, :ω_r, :rh)

all_rh = plot(p1_rh, p2_rh, p3_rh, p4_rh, p5_rh, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rh.svg")

p1_rt = plot_different_bp_lines(atp_range, params1, :atp, :rt)
p2_rt = plot_different_bp_lines(kin_range, params1, :kin, :rt)
p3_rt = plot_different_bp_lines(lam_range, params1, :lam, :rt)
p4_rt = plot_different_bp_lines(wab_range, params1, :ω_ab, :rt)
p5_rt = plot_different_bp_lines(wr_range, params1, :ω_r, :rt)

all_rt = plot(p1_rt, p2_rt, p3_rt, p4_rt, p5_rt, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rt, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rt.svg")

p1_rd = plot_different_bp_lines(atp_range, params1, :atp, :rd)
p2_rd = plot_different_bp_lines(kin_range, params1, :kin, :rd)
p3_rd = plot_different_bp_lines(lam_range, params1, :lam, :rd)
p4_rd = plot_different_bp_lines(wab_range, params1, :ω_ab, :rd)
p5_rd = plot_different_bp_lines(wr_range, params1, :ω_r, :rd)

all_rd = plot(p1_rd, p2_rd, p3_rd, p4_rd, p5_rd, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Vary param bistability plot")
savefig(all_rd, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_rd.svg")


l = @layout [a b; c d; e f h]
atp_plots = plot(p1_rma, p1_rtca, p1_rmr, p1_rtcr, p1_rh, p1_rt, p1_rd, layout=l, size=(1200,800), plot_title="ATP = $([round.(i; digits=4) for i in atp_range])", titlefontsize=6, left_margin=4Plots.mm)
kin_plots = plot(p2_rma, p2_rtca, p2_rmr, p2_rtcr, p2_rh, p2_rt, p2_rd, layout=l, size=(1200,800), plot_title="kin = $([round.(i; digits=4) for i in kin_range])", titlefontsize=6, left_margin=4Plots.mm)
lam_plots = plot(p3_rma, p3_rtca, p3_rmr, p3_rtcr, p3_rh, p3_rt, p3_rd, layout=l, size=(1200,800), plot_title="λ = $([round.(i; digits=4) for i in lam_range])", titlefontsize=6, left_margin=4Plots.mm)
wab_plots = plot(p4_rma, p4_rtca, p4_rmr, p4_rtcr, p4_rh, p4_rt, p4_rd, layout=l, size=(1200,800), plot_title="ω_ab = $([round.(i; digits=4) for i in wab_range])", titlefontsize=6, left_margin=4Plots.mm)
wr_plots = plot(p5_rma, p5_rtca, p5_rmr, p5_rtcr, p5_rh, p5_rt, p5_rd, layout=l, size=(1200,800), plot_title="ω_r = $([round.(i; digits=4) for i in wr_range])", titlefontsize=6, left_margin=4Plots.mm)

savefig(atp_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/atp_plots.svg")
savefig(kin_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/kin_plots.svg")
savefig(lam_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/lam_plots.svg")
savefig(wab_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wab_plots.svg")
savefig(wr_plots, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wr_plots.svg")



using PlotlyJS
p1 = PlotlyJS.scatter(
    bsp,
    x=:atp, y=:wab, z=:kdam1,
    type="scatter3d", mode="markers", name="bf1"
)
p2 = PlotlyJS.scatter(
    bsp,
    x=:atp, y=:wab, z=:kdam2,
    type="scatter3d", mode="markers"
)

layout = Layout(
    scene=attr(
        zaxis=attr(
            range=[0,1],
            title="kdam"),
        xaxis=attr(
            title="ATP"),    
        yaxis=attr(
            title="ω_ab"),)
)
plot([p1,p2], layout)




function plot_3d_bf_against_kdam(param_range1, param1, param_range2, param2)
    df = double_param_vary(atp_range, :atp, kin_range, :kin)
    bsp = df[df.bs .== :bp, :]
    p1 = PlotlyJS.scatter(
        bsp,
        x=:atp, y=:wab, z=:kdam1,
        type="scatter3d", mode="markers", label="bf1"
    )
    p2 = PlotlyJS.scatter(
        bsp,
        x=:atp, y=:wab, z=:kdam2,
        type="scatter3d", mode="markers"
    )
    
    layout = Layout(
        scene=attr(
            zaxis=attr(
                range=[0,1],
                title="kdam"),
            xaxis=attr(
                title="ATP"),    
            yaxis=attr(
                title="ω_ab"),)
    )
    plot([p1,p2], layout)

end

