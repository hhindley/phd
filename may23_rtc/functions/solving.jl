
using DifferentialEquations
# module Solving

all_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
all_species_OD = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD]
all_species_atp = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp]

# export solcb, sol_with_t, sol, get_all_curves, get_all_curves_df, check_get_ssval, get_curve, get_ssval

function solcb(model, init, tspan, params, callback)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params, callback=callback)
    solu = solve(prob, Rodas4())
    return solu
end

function sol_with_t(model, init, params, tspan, t)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4(), saveat=t)
    return solu
end

function sol(model, init, tspan, params)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4())
    # solu = solve(prob, Rodas5(), isoutofdomain=(y,p,t)->any(x->x<0,y), abstol=1e-15, reltol=1e-12);
    # @show solu.retcode
    # solu = solve(prob, alg_hints=[:auto], isoutofdomain=(y,p,t)->any(x->x<0,y))#, abstol=1e-15, reltol=1e-12);

    return solu
end

function get_all_curves(sol, species) 
    dict_res = OrderedDict(name => [] for name in species)
    for (i,j) in zip(values(dict_res), species)
        push!(i, get_curve(sol, j))
    end
    return ((dict_res))
end

function get_all_curves_df(sol, species) 
    df_res = DataFrame(name => [] for name in species)
    for (i,j) in zip(values(df_res), species)
        push!(i, get_curve(sol, j))
    end
    return df_res
end

function check_get_ssval(sol, species, n=3)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    end

    ss_vals = Dict(:rm_a=>[], :rtca=>[], :rm_b=>[], :rtcb=>[], :rm_r=>[], :rtcr=>[], :rh=>[], :rd=>[], :rt=>[])
    for i in collect(0:2)
        alpha = df[!,:rt][end-i]/kr 
        fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
        ra = fa*df[!,:rtcr][end-i]
        
        Vinit = ra*Vmax_init*atp/(Km_init+atp)
        tscr_el_a = ω_ab*atp/(θtscr+atp)
        tscr_a = Vinit*tscr_el_a
        tscr_el_b = ω_ab*atp/(θtscr+atp)
        tscr_b = Vinit*tscr_el_b
        tscr_r = ω_r*atp/(θtscr+atp)

        tlr_el = g_max*atp/(θtlr+atp)

        rtca1 = (atp*df[!,:rtca][end-i])/(atp+(km_a*df[!,:rd][end-i])) 
        rtcb1 = (atp*df[!,:rtcb][end-i])/(atp+(km_b*df[!,:rt][end-i])) 

        Vrep = krep*rtcb1*df[!,:rt][end-i]
        Vdam = kdam*df[!,:rh][end-i]
        Vinflux = kin*tlr_el 
        Vtag = ktag*rtca1*df[!,:rd][end-i]

        push!(ss_vals[:rm_a], tscr_a - lam(lam[end-i])*df[!,:rm_a][end-i] - d*df[!,:rm_a][end-i])
        push!(ss_vals[:rtca], (1/na)*df[!,:rh][end-i]*df[!,:rm_a][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtca][end-i])  
        push!(ss_vals[:rm_b], tscr_b - lam(lam[end-i])*df[!,:rm_b][end-i] - d*df[!,:rm_b][end-i])
        push!(ss_vals[:rtcb], (1/nb)*df[!,:rh][end-i]*df[!,:rm_b][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtcb][end-i])
        push!(ss_vals[:rm_r], tscr_r - lam(lam[end-i])*df[!,:rm_r][end-i] - d*df[!,:rm_r][end-i])
        push!(ss_vals[:rtcr], (1/nr)*df[!,:rh][end-i]*df[!,:rm_r][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtcr][end-i])
        push!(ss_vals[:rh], Vrep - Vdam + Vinflux - lam(lam[end-i])*df[!,:rh][end-i])
        push!(ss_vals[:rd], Vdam - Vtag - kdeg*df[!,:rd][end-i] - lam(lam[end-i])*df[!,:rd][end-i])
        push!(ss_vals[:rt], Vtag - Vrep - lam(lam[end-i])*df[!,:rt][end-i])
    end
    return std(ss_vals[species])
    # return (std(df[length(sol.t)-n:length(sol.t), species])/mean(df[length(sol.t)-n:length(sol.t), species]))*100
end

function get_curve(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp])
    end
    species = df[:, species]
    return species
end

function get_ssval(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp])
    end
    species = df[end, species]
    return species
end

# end

function ss_init_vals(sol)
    rm_a = get_ssval(sol, :rm_a)
    rm_b = get_ssval(sol, :rm_b)
    rm_r = get_ssval(sol, :rm_r)
    rtca = get_ssval(sol, :rtca)
    rtcb = get_ssval(sol, :rtcb)
    rtcr = get_ssval(sol, :rtcr)
    rh = get_ssval(sol, :rh)
    rt = get_ssval(sol, :rt)
    rd = get_ssval(sol, :rd)

    return [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt]
end