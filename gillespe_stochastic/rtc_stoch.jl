using DifferentialEquations, JumpProcesses, Random, Catalyst, Distributions, DataFrames, GLMakie

include("/Users/s2257179/phd/general_funcs/solving.jl")
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))

Vref = 1
v_sf(u) = Vref/u[vidx(:V)]

sf1(u) = (1e6/(6.022e23*(u[vidx(:V)]*1e-15)))
atp(u) = p[pidx(:atp)]/sf1(u)

alpha(u) = u[vidx(:rt)]/p[pidx(:kr)] * v_sf(u) # unitless
fa(u) = (1+alpha(u))^6/(p[pidx(:L)]*((1+p[pidx(:c)]*alpha(u))^6)+(1+alpha(u))^6)  # unitless 
ra(u) = fa(u)*u[vidx(:rtcr)] # uM 

# transcription
Voc(u) = (p[pidx(:Vmax_init)]*atp(u)/((p[pidx(:Km_init)]/v_sf(u))+atp(u)))
sig_o(u) = ra(u)*Voc(u)/p[pidx(:k_diss)] # uM

tscr_ab(u) = (sig_o(u)*p[pidx(:ω_ab)]*atp(u)/(p[pidx(:θtscr)]/v_sf(u)+atp(u)) )# uM min-1
tscr_r(u) = (p[pidx(:ω_r)]/v_sf(u)*atp(u)/((p[pidx(:θtscr)]/v_sf(u))+atp(u))) # uM min-1

tlr_el(u) = (p[pidx(:g_max)]*atp(u)/((p[pidx(:θtlr)]/v_sf(u))+atp(u)))
tlr_a(u) = (1/p[pidx(:na)])*p[pidx(:kc)]*u[vidx(:rh)]*u[vidx(:rm_a)]*tlr_el(u)
tlr_b(u) = (1/p[pidx(:nb)])*p[pidx(:kc)]*u[vidx(:rh)]*u[vidx(:rm_b)]*tlr_el(u)
tlr_r(u) = (1/p[pidx(:nr)])*p[pidx(:kc)]*u[vidx(:rh)]*u[vidx(:rm_r)]*tlr_el(u)

# # ribosomes
Vrep(u) = u[vidx(:rtcb)]*u[vidx(:rt)]*p[pidx(:krep)]/(u[vidx(:rt)]+p[pidx(:km_b)]/v_sf(u)) # uM min-1 
Vdam(u) = p[pidx(:kdam)]*u[vidx(:rh)] # uM min-1
Vinflux(u) = p[pidx(:kin)] * p[pidx(:g_max)]*atp(u)/(p[pidx(:θtlr)]/v_sf(u)+atp(u)) # uM min-1 
Vtag(u) = u[vidx(:rtca)]*u[vidx(:rd)]*p[pidx(:ktag)]/(u[vidx(:rd)]+p[pidx(:km_a)]/v_sf(u)) # uM min-1 

rate1(u,p,t) = tscr_ab(u)
rate2(u,p,t) = tscr_r(u)
rate3(u,p,t) = tlr_a(u) * v_sf(u)
rate4(u,p,t) = tlr_b(u) * v_sf(u)
rate5(u,p,t) = tlr_r(u) * v_sf(u)
rate6(u,p,t) = Vinflux(u) / v_sf(u)
rate7(u,p,t) = Vdam(u)
rate8(u,p,t) = Vtag(u)
rate9(u,p,t) = Vrep(u)
rate10(u,p,t) = p[pidx(:kdeg)] * u[vidx(:rd)]
rate11(u,p,t) = p[pidx(:d)] * u[vidx(:rm_a)]
rate12(u,p,t) = p[pidx(:d)] * u[vidx(:rm_b)]
rate13(u,p,t) = p[pidx(:d)] * u[vidx(:rm_r)]

function affect1!(integrator)
    integrator.u[vidx(:rm_a)] += 1       # * → rm_a  
    integrator.u[vidx(:rm_b)] += 1       # * → rm_b 
    nothing
end

function affect2!(integrator)
    integrator.u[vidx(:rm_r)] += 1       # * → rm_r 
    nothing
end

function affect3!(integrator)
    integrator.u[vidx(:rtca)] += 1       # rm_a + rh → rtca
    nothing
end

function affect4!(integrator)
    integrator.u[vidx(:rtcb)] += 1       # rm_b + rh → rtcb
    nothing
end

function affect5!(integrator)
    integrator.u[vidx(:rtcr)] += 1       # rm_r + rh → rtcr
    nothing
end

function affect6!(integrator)
    integrator.u[vidx(:rh)] += 1       # * → rh
    nothing
end

function affect7!(integrator)
    integrator.u[vidx(:rh)] -= 1       # rh → rd
    integrator.u[vidx(:rd)] += 1       # rh → rd
    nothing
end

function affect8!(integrator)
    integrator.u[vidx(:rd)] -= 1       # rd + rtca → rt
    integrator.u[vidx(:rt)] += 1       # rd + rtca → rt
    nothing
end

function affect9!(integrator)
    integrator.u[vidx(:rt)] -= 1       # rt + rtcb → rh
    integrator.u[vidx(:rh)] += 1       # rt + rtcb → rh
    nothing
end

function affect10!(integrator)
    integrator.u[vidx(:rd)] -= 1       # rd → *
    nothing
end

function affect11!(integrator)
    integrator.u[vidx(:rm_a)] -= 1       # rm_a → *
    nothing
end

function affect12!(integrator)
    integrator.u[vidx(:rm_b)] -= 1       # rm_a → *
    nothing
end

function affect13!(integrator)
    integrator.u[vidx(:rm_r)] -= 1       # rm_a → *
    nothing
end

react1 = VariableRateJump(rate1, affect1!)
react2 = VariableRateJump(rate2, affect2!)
react3 = VariableRateJump(rate3, affect3!)
react4 = VariableRateJump(rate4, affect4!)
react5 = VariableRateJump(rate5, affect5!)
react6 = VariableRateJump(rate6, affect6!)
react7 = VariableRateJump(rate7, affect7!)
react8 = VariableRateJump(rate8, affect8!)
react9 = VariableRateJump(rate9, affect9!)
react10 = VariableRateJump(rate10, affect10!)
react11 = VariableRateJump(rate11, affect11!)
react12 = VariableRateJump(rate12, affect12!)
react13 = VariableRateJump(rate13, affect13!)

jumpset = JumpSet(react1, react2, react3, react4, react5, react6, react7, react8, react9, react10, react11, react12, react13)

function get_X0(indV, ssvals)
    X0 = zeros(10)

    X0[vidx(:rm_a)] = ssvals[1]
    X0[vidx(:rtca)] = ssvals[2]
    X0[vidx(:rm_b)] = ssvals[3]
    X0[vidx(:rtcb)] = ssvals[4]
    X0[vidx(:rm_r)] = ssvals[5]
    X0[vidx(:rtcr)] = ssvals[6]
    X0[vidx(:rh)] = ssvals[7]
    X0[vidx(:rd)] = ssvals[8]
    X0[vidx(:rt)] = ssvals[9]
    X0[vidx(:V)] = 1 

    return X0
end

p = collect(get_par(indP))
u0 = collect(get_X0(indV, init_molec))
tspan = (0.0,500.0)

function volume_ode(du, u, p, t)
    du[vidx(:V)] = u[vidx(:V)] * p[pidx(:lam)]
    nothing
end
div_times = [i*(log(2)/p[pidx(:lam)]) for i in 1:floor(tspan[2]/(log(2)/p[pidx(:lam)]))]
condition_div(u,t,integrator) = t in div_times

function affect_division!(integrator)
    pDiv=rand(Beta(200.,200.))
    for i in 1:length(integrator.u)-2
        n = integrator.u[i]
        integrator.u[i] = n > 0 ? float(rand(Binomial(round(Int, n), pDiv))) : 0
    end
    integrator.u[vidx(:V)] = 1 #integrator.u[vidx(:V)]/2
    nothing
end
division_cb = DiscreteCallback(condition_div, affect_division!) # is called at every time step during integration, if t is greater than division time then make division = true and terminate integration

p[pidx(:kdam)] = 0.0

prob_vol = ODEProblem(volume_ode, u0, tspan, p)
jump_prob = JumpProblem(prob_vol, Direct(), jumpset)

solu = solve(jump_prob, Tsit5(), callback=division_cb, tstops=div_times, abstol=1e-15, reltol=1e-12)

df = DataFrame(solu)
lines(df.timestamp, sf*(df.value7./df.value10))


df = create_solu_df(solu, [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :V])

lines(df.time, df.rtca)