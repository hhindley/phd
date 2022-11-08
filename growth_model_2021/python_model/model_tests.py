from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.integrate import solve_ivp

def model(y, t, params):

    # LHS
    E = y[0]
    S = y[1]
    ES = y[2]
    P = y[3]
    
    # params
    k1 = params[0]
    k2 = params[1]
    k3 = params[2]
    
    # rates
    v1 = k1*S*E
    v2 = k2*ES
    v3 = k3*ES

    # RHS (dydt)
    dydt = np.zeros(4)
    # first reaction
    dydt[0] = dydt[0] - v1
    dydt[1] = dydt[1] - v1
    dydt[2] = dydt[2] + v1
    #second reaction
    dydt[0] = dydt[0] + v2
    dydt[1] = dydt[1] + v2
    dydt[2] = dydt[2] - v2
    #third reaction
    dydt[0] = dydt[0] + v3
    dydt[2] = dydt[2] - v3
    dydt[3] = dydt[3] + v3
    
    return dydt