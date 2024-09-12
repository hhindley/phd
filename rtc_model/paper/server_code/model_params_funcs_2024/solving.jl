function steady_states(model, init, params)
    prob = SteadyStateProblem(model, init, params; jac=true)
    return solve(prob, DynamicSS())#; abstol=1e-15, reltol=1e-12)
end 

function sol(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params; jac=true)
    if nameof(model) == :combined_model || nameof(model) == :reduced_combined_model
        # solu = solve(prob, Rosenbrock23())#, abstol=1e-10, reltol=1e-10)
        solu = solve(prob, QNDF(), abstol=1e-9, reltol=1e-6);
    # elseif nameof(model) == :rtc_trna_inhib_model
        # solu = solve(prob, TRBDF2())
    elseif nameof(model) == :growth_model 
        solu = solve(prob, Rodas4());
    else
        solu = solve(prob, Rodas4())
    end
    # solu = solve(prob, Rosenbrock23(), isoutofdomain=(y,p,t)->any(x->x<0,y), abstol=1e-18, reltol=1e-4);
    # @show solu.retcode
    # solu = solve(prob, alg_hints=[:auto], isoutofdomain=(y,p,t)->any(x->x<0,y))#, abstol=1e-15, reltol=1e-12);

    return solu
end

function calc_lam(params, ssvals_dict, gm_or_comb)
    if gm_or_comb == :comb    
        gamma = @. params[gmax] * ssvals_dict[:a]/(params[Kgamma] + ssvals_dict[:a])
        ttrate = @. (ssvals_dict[:c_q] + ssvals_dict[:c_rh] + ssvals_dict[:c_t] + ssvals_dict[:c_m] + ssvals_dict[:c_R] + ssvals_dict[:c_A] + ssvals_dict[:c_B])*gamma
    else
        gamma = @. params[gmax] * ssvals_dict[:a]/(params[Kgamma] + ssvals_dict[:a])
        ttrate = @. (ssvals_dict[:cq] + ssvals_dict[:cr] + ssvals_dict[:ct] + ssvals_dict[:cm])*gamma
    end
    return @. ttrate/params[M]
end

function calc_rmf(params, ssvals_dict, gm_or_comb)
    if gm_or_comb == :comb
        # return params[nrh]*(ssvals_dict[:rh] + ssvals_dict[:rt] + ssvals_dict[:rd] + ssvals_dict[:c_rh] + ssvals_dict[:c_t] + ssvals_dict[:c_m] + ssvals_dict[:c_q] + ssvals_dict[:c_A] + ssvals_dict[:c_B] + ssvals_dict[:c_R] + ssvals_dict[:z_rh] + ssvals_dict[:z_t] + ssvals_dict[:z_m] + ssvals_dict[:z_q] + ssvals_dict[:z_A] + ssvals_dict[:z_R] + ssvals_dict[:z_B])/params[M]
        return params[nrh]*(ssvals_dict[:rh] + ssvals_dict[:rt] + ssvals_dict[:rd] + ssvals_dict[:c_rh] + ssvals_dict[:c_t] + ssvals_dict[:c_m] + ssvals_dict[:c_q] + ssvals_dict[:c_A] + ssvals_dict[:c_B] + ssvals_dict[:c_R])/params[M]
    else
        return params[nr]*(ssvals_dict[:r] + ssvals_dict[:cr] + ssvals_dict[:ct] + ssvals_dict[:cm] + ssvals_dict[:cq] + ssvals_dict[:zmr] + ssvals_dict[:zmt] + ssvals_dict[:zmm] + ssvals_dict[:zmq])/params[M]
    end
end

function plot_solu(df)
    return plot([scatter(x=df.time, y=col, name="$(names(df)[i])") for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log",  yaxis_tickformat=".2e"))
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

function create_solu_df(sol, all_species)
    df = DataFrame(sol)
    names = [:time;all_species]
    rename!(df, names)
    return df
end

function get_curve(df, species)
    species = df[:, species]
    return species
end

function get_ssval(df, specie)
    # df = create_solu_df(sol, all_species)
    species = df[end, specie]
    return species
end
function get_all_ssvals(sol, all_species)
    df = create_solu_df(sol, all_species)
    ssvals=[]
    for i in all_species
        push!(ssvals, get_ssval(df, i))
    end
    return ssvals
end

function ss_init_vals(df, all_species)
    return [get_ssval(df, i) for i in all_species]
end


function all_var_df(sol, all_species)
    df = DataFrame(sol)
    names = [:time;all_species]
    rename!(df, names)
    alpha = @. df.rt/kr
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = @. fa*df.rtcr
    Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_ab = Vinit*tscr_el_a
    tscr_r = ω_r*atp/(θtscr+atp)
    tlr_el = g_max*atp/(θtlr+atp)
    tlr_r = @. (1/nr)*df.rh*df.rm_r*tlr_el
    tlr_a = @. (1/na)*df.rh*df.rm_a*tlr_el
    tlr_b = @. (1/nb)*df.rh*df.rm_b*tlr_el
    rtca1 = @. (atp*df.rtca)/(atp+(km_a*df.rd)) 
    rtcb1 = @. (atp*df.rtcb)/(atp+(km_b*df.rt)) 
    Vrep = @. krep*rtcb1*df.rt
    Vdam = @. kdam*df.rh
    Vinflux = kin*tlr_el
    Vtag = @. ktag*rtca1*df.rd
    all_df = DataFrame(:alpha=>alpha,:fa=>fa,:ra=>ra,:Vinit=>Vinit,:tscr_ab=>tscr_ab, :tscr_r=>tscr_r, :tlr_r=>tlr_r, :tlr_b=>tlr_b, :tlr_a=>tlr_a, :rtca1=>rtca1, :rtcb1=>rtcb1, :Vrep=>Vrep, :Vdam=>Vdam, :Vinflux=>Vinflux, :Vtag=>Vtag)
    return all_df
end

function var_param(model, kdam, params_rtc1, kdam_range, init1)
    # new_params = deepcopy(params_rtc1)
    ssvals=[]
    @show params_rtc1[c], params_rtc1[L]
    for i in kdam_range
        params_rtc1[kdam] = i
        solu_rtc = sol(model, init1, tspan, params_rtc1)
        df = create_solu_df(solu_rtc, species_rtc)
        ss = [i[end] for i in eachcol(df[:,2:end])]
        # ss = steady_states(test, init_rtc, new_params)
        push!(ssvals, ss)
    end
    df_ssvals = DataFrame(vcat(transpose(ssvals)...), :auto)
    rename!(df_ssvals, species_rtc)
    return df_ssvals
end