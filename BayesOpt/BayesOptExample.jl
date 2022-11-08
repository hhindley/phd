using BayesOpt
using Plots

config = ConfigParameters()         # initiates parameters 
# print(config)

config.n_iterations = 50
config.n_iter_relearn = 5
config.n_init_samples = 2       # how to change parameter values from original
# print(config.n_iterations)

# set_kernel!(config, "kMaternARD5")  # calls set_kernel of the C API

# config.sc_type = SC_MAP
# f(x) =  sum(x .^ 2) # defines the function 
function testfunc(Xin)
    total = 5
    for value in Xin
        total=total+(value-0.33)*(value-0.33)
    end
    return total
end
# x = range(-2,2, length=101)
# plot(x,@. sum(x^2))
n = 5
lb = zeros((n,)); ub = ones((n,));

function bo()
    @time optimizer, optimum = bayes_optimization(testfunc, lb, ub, config)
end

bo()
# print(optimizer)
# print(optimum)
# scatter!(optimizer)
# plot!([optimum])

# print(lowerbound)
# plot!(lowerbound)
# plot!(upperbound)

# scatter!([0,1,2,3,4])







config = ConfigParameters()         # calls initialize_parameters_to_default of the C API
set_kernel!(config, "kMaternARD5")  # calls set_kernel of the C API
config.sc_type = SC_MAP
f(x) = sum(x .^ 2)
lowerbound = [-2.]; upperbound = [2.]
optimizer, optimum = bayes_optimization(f, lowerbound, upperbound, config)

x = -20:20
plot(x.^2)

vline!(optimizer)

x=0:0.01:1
total = 5
res = []
for i in x
    total = total + (i-0.33) * (i-0.33)
    push!(res, total)
end

print(res)
plot(x, res)