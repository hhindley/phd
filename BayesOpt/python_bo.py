import sys
import bayesopt
from bayesoptmodule import BayesOptContinuous
import numpy as np

from time import clock

# Python3 compat
if hasattr(__builtins__, 'raw_input'):
    input = raw_input


# Function for testing.
def testfunc(Xin):
    total = 5.0
    print(Xin)
    for value in Xin:
        total = total + (value - 0.33) * (value - 0.33)

    return total

def new(x):
    return (x**2)

# Class for OO testing.
class BayesOptTest(BayesOptContinuous):
    def evaluateSample(self,Xin):
        return testfunc(Xin)
        


# Let's define the parameters
# For different options: see parameters.h and cpp
# If a parameter is not define, it will be automatically set
# to a default value.
params = {}
params['n_iterations'] = 50
params['n_iter_relearn'] = 5
params['n_init_samples'] = 2

print("Callback implementation")

n = 5                     # n dimensions
lb = np.zeros((n,))
ub = np.ones((n,))

start = clock()
mvalue, x_out, error = bayesopt.optimize(testfunc, n, lb, ub, params) # standard 

print("Result", mvalue, "at", x_out)
print("Running time:", clock() - start, "seconds")

# input('Press INTRO to continue')

# print("OO implementation") # object oriented implementation
# bo_test = BayesOptTest(n)
# bo_test.parameters = params
# bo_test.lower_bound = lb
# bo_test.upper_bound = ub

# start = clock()
# mvalue, x_out, error = bo_test.optimize()

# print("Result", mvalue, "at", x_out)
# print("Running time:", clock() - start, "seconds")
# input('Press INTRO to continue')

# print("Callback discrete implementation") # discrete 
# x_set = np.random.rand(100, n)
# start = clock()

# mvalue, x_out, error = bayesopt.optimize_discrete(testfunc, x_set, params)

# print("Result", mvalue, "at", x_out)
# print("Running time:", clock() - start, "seconds")

# value = np.array([testfunc(i) for i in x_set])
# print("Optimum", value.min(), "at", x_set[value.argmin()])
