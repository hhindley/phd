using DifferentialEquations, Plots, DataFrames

include("rtcR_models.jl");
include("params.jl");
include("functions.jl");

# prob1 = ODEProblem(comparison_rtc!, init, tspan, comp_params);
# sol1 = solve(prob1, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
# plot(sol1)

a_range = range(0, 100, length=101);

function rtcr_change(a_range)
    ss_p = [];
    comp_params = [k1, k2, k3, k4, katp, a_0, b_0];
    for i in a_range
        comp_params[6] = i
        prob = ODEProblem(comparison_rtc!, init, tspan, comp_params);
        sol = solve(prob, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
        b, c, d, p = separate_species!(species2, sol);
        push!(ss_p, last(p))
    end
    return ss_p
end

p1 = plot(a_range, rtcr_change(a_range), title="Best model RtcR", xlabel="RtcR", ylabel="\\sigma*", legend=false);

function rt_change(a_range)
    ss_p = [];
    comp_params = [k1, k2, k3, k4, katp, a_0, b_0];
    init = [b_0, c_0, d_0, p_0];
    for i in a_range
        comp_params[7] = i
        init[1] = i
        prob = ODEProblem(comparison_rtc!, init, tspan, comp_params);
        sol = solve(prob, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
        b, c, d, p = separate_species!(species2, sol);
        push!(ss_p, last(p))
    end
    return ss_p
end

p2 = plot(a_range, rt_change(a_range), title="Best model Rt", xlabel="Rt", ylabel="\\sigma*", legend=false);


function rtcr_full(a_range, change)
    ss_p = [];
    init_full = [a_0, b_0, c_0, d_0, e_0, p_0];
    for i in a_range
        init_full[change] = i
        prob = ODEProblem(hill!, init_full, tspan, params_full);
        sol = solve(prob, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
        b, c, d, p = separate_species2!(species_full2, sol);
        push!(ss_p, last(p))
    end
    return ss_p
end

p3 = plot(a_range, rtcr_full(a_range, 1), title="Full hill model RtcR", xlabel="RtcR", ylabel="\\sigma*", legend=false);

p4 = plot(a_range, rtcr_full(a_range, 2), title="Full hill model Rt", xlabel="Rt", ylabel="\\sigma*", legend=false);

plot(p1, p2, p3, p4, layout=(2,2), size=(600,600))



# prob1 = ODEProblem(comparison_rtc!, init, tspan, comp_params);
# sol1 = solve(prob1, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
# pl1 = plot(sol1);

# prob = ODEProblem(hill!, init_full, tspan, comp_params);
# sol = solve(prob, alg=Rodas4(), abstol=1e-6, reltol=1e-6);
# pl2 = plot(sol);

# plot(pl1, pl2, layout=(2,1))