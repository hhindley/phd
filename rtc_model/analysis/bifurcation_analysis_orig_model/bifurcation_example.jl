using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

# vector field
function TMvf!(dz, z, p, t)
	@unpack J, α, E0, τ, τD, τF, U0 = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	dz[1] = (-E + SS1) / τ
	dz[2] =	(1.0 - x) / τD - u * x * E
	dz[3] = (U0 - u) / τF +  U0 * (1.0 - u) * E
	dz
end

# out of place method
TMvf(z, p) = TMvf!(similar(z), z, p, 0)

# parameter values
par_tm = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007)

# initial condition
z0 = [0.238616, 0.982747, 0.367876]

# Bifurcation Problem
prob = BifurcationProblem(TMvf, z0, par_tm, (@lens _.E0);
recordFromSolution = (x, p) -> (E = x[1], x = x[2], u = x[3]),)

# continuation options
opts_br = ContinuationPar(pMin = -10.0, pMax = -0.9,
# parameters to have a smooth result
ds = 0.04, dsmax = 0.05,)

# continuation of equilibria
br = continuation(prob, PALC(tangent=Bordered()), opts_br;
plot = true, normC = norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft)

