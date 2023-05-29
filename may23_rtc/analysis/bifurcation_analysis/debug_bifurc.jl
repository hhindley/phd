using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
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
    krep = 137., ktag = 9780., atp = 4000., km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.054, ω_ab = 4., ω_r = 2e-7, 
    kdam =  0.01, lam = 0.033)


# initial condition
initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)

wab_range = range(1,2)
wr_range = range(1,2)
atp_range = range(1,2)
kin_range = range(1,2)
lam_range = range(1,2)
kdam_range = range(1,2)

bistab = []
for i in lam_range
    params = merge(params, (lam=i,))
    prob = BifurcationProblem(rtc_mod, initial, params, (@lens _.kdam);
    recordFromSolution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
    opts_br = ContinuationPar(pMin = 0., pMax = 1.0, dsmax = 0.01, detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25, nev = 2, maxSteps = 1000,)
    br = continuation(prob, PALC(), opts_br; plot = false, bothside=true, normC = norminf)
    for b in range(1,length(br.specialpoint))
        if br.specialpoint[b].type == :endpoint
            push!(bistab, ("lam = $i", br.specialpoint[i].type))
        else
            Nothing
        end
    end  
end
bistab