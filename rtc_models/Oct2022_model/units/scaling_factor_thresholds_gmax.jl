θ = 4.38; max = 1260; thr = 7;

# to go from molecules/cell to μM/min
SF = 1e6/(6.022e23*1e-15)

a_SF = 22000

θtscr = θ*SF*a_SF
θtlr = thr*SF*a_SF

# max rate of translation 
g_max = max*SF