import numpy as np
import matplotlib as plt
from scipy.integrate import solve_ivp

from init_rp import params, rates, init
from ode_rhs import model



tspan = np.arange(0, 1000)

sol = solve_ivp(lambda t, y: model(t, y, params, rates),
                tspan, init, method='Radau')

plt.plot(sol.t, sol.y[0])
