using DifferentialEquations
using JumpProcesses
using Random, Catalyst

β = 0.1 / 1000.0
ν = 0.01;
lam = 0.014
p = (β=β, ν=ν, lam=lam)
rate1(u, p, t) = p.β * u[1] * u[2]  # β*S*I
function affect1!(integrator)
    integrator.u[1] -= 1         # S -> S - 1
    integrator.u[2] += 1         # I -> I + 1
    nothing
end
jump = ConstantRateJump(rate1, affect1!)

rate2(u, p, t) = p[2] * u[2]         # ν*I
function affect2!(integrator)
    integrator.u[2] -= 1        # I -> I - 1
    integrator.u[3] += 1        # R -> R + 1
    nothing
end
jump2 = ConstantRateJump(rate2, affect2!)


u0 = [999, 10, 0]
tspan = (0.0,250.0)

prob = DiscreteProblem(u0, tspan, p)

jump_prob = JumpProblem(prob, Direct(), jump, jump2)

solu = solve(jump_prob, SSAStepper())

using DataFrames, GLMakie
df = DataFrame(solu)

lines(df.timestamp, df.value1)
lines!(df.timestamp, df.value2)
lines!(df.timestamp, df.value3)




function f(du, u, p, t)
    du[4] = u[4] * lam
    nothing
end
div_times = [i*(log(2)/lam) for i in 1:floor(250/(log(2)/lam))]
condition_div(u,t,integrator) = t in div_times

function affect_division!(integrator)
    integrator.u[1] /= 2
    integrator.u[2] /= 2
    integrator.u[3] /= 2
    integrator.u[4] /= 2
    nothing
end
division_cb = DiscreteCallback(condition_div, affect_division!) # is called at every time step during integration, if t is greater than division time then make division = true and terminate integration

u₀ = [990.0, 10.0, 0.0, 1.0]
prob = ODEProblem(f, u₀, tspan, p)

jump_prob = JumpProblem(prob, Direct(), jump, jump2)
sol = solve(jump_prob, Tsit5(), callback=division_cb, tstops=div_times)

df1 = DataFrame(sol)
p1=lines(df1.timestamp, df1.value4)
DataInspector(p1)
lines(df1.timestamp, df1.value1)
lines!(df1.timestamp, df1.value2)
lines!(df1.timestamp, df1.value3)
lines!(df1.timestamp, df1.value4)


