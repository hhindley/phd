using DifferentialEquations, GLMakie, DataFrames, Catalyst

include("/Users/s2257179/phd/general_funcs/solving.jl")

# here we order S = 1, I = 2, and R = 3
# substrate stoichiometry:
substoich = [[1 => 1, 2 => 1],    # 1*S + 1*I
    [2 => 1]]                     # 1*I
# net change by each jump type
netstoich = [[1 => -1, 2 => 1],   # S -> S-1, I -> I+1
    [2 => -1, 3 => 1]]            # I -> I-1, R -> R+1
# rate constants for each jump
p = (0.1 / 1000, 0.01)

# p[1] is rate for S+I --> 2I, p[2] for I --> R
pidxs = [1, 2]

maj = MassActionJump(substoich, netstoich; param_idxs = pidxs)

u₀ = [999, 1, 0]       #[S(0),I(0),R(0)]
tspan = (0.0, 250.0)
dprob = DiscreteProblem(u₀, tspan, p)

# use the Direct method to simulate
jprob = JumpProblem(dprob, maj)

# solve as a pure jump process, i.e. using SSAStepper
solu = solve(jprob)

df = create_solu_df(solu, [:S, :I, :R])

lines(df.time, df.S)
lines!(df.time, df.I)
lines!(df.time, df.R)


