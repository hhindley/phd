using Catalyst, DataFrames, DifferentialEquations, GLMakie
include("/Users/s2257179/phd/general_funcs/solving.jl")

sir_model = @reaction_network begin
    β, S + I --> 2I
    ν, I --> R
    lam, V --> V+(V*lam)
end

p = (:β => 0.1 / 1000, :ν => 0.01, :lam => 0.014)
u₀ = [:S => 990, :I => 10, :R => 0, :V => 1]
tspan = (0.0, 250.0)
prob = DiscreteProblem(sir_model, u₀, tspan, p)

jump_prob = JumpProblem(sir_model, prob, Direct())

solu = solve(jump_prob, SSAStepper())

df = create_solu_df(solu, [:S, :I, :R, :V])

lines(df.time, df.V)