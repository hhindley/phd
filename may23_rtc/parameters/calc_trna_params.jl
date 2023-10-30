# to go from molecules/cell to μM/min
SF = 1e6/(6.022e23*1e-15)

#number of ribosomes/cell at different growth rates 
r1 = 6800 # gr = 0.6/h
r2 = 13500 # gr = 1/h
r3 = 26300 # gr = 1.5/h
r4 = 45100 # gr = 2/h
r5 = 72000 # gr = 2.5/h

# concentration of ribosomes at those above growth rates 
rc1 = r1*SF
rc2 = r2*SF
rc3 = r3*SF
rc4 = r4*SF
rc5 = r5*SF

# relationship between number of ribosomes and tRNAs at different growth rates 
grs=[0.4,0.7,1.07,1.6,2.5]
trna_per_ribo=[12.8,9.5,7.8,7.6,6.8]
plot(scatter(x=grs,y=trna_per_ribo), Layout(xaxis_title="Growth rate (1/h)", yaxis_title="Number of tRNA molecules per ribosome"))

# number of tRNAs for each growth as above 
trnas1 = rc1*trna_per_ribo[1]
trnas2 = rc2*trna_per_ribo[2]
trnas3 = rc3*trna_per_ribo[3]
trnas4 = rc4*trna_per_ribo[4]
trnas5 = rc5*trna_per_ribo[5]

# kin for normal model - taken from the ode for rh assuming no damage or repair
g_max = 2.0923; atp = 3578.9473684210525; θtlr = 255.73; rh = 120; λ = 0.014;
λ = 0.0475;
kin = λ*rh/(g_max*atp/(θtlr+atp))

# kin for tRNA model - taken from the ode for tRNA assuming no damage or repair 
trna = 135.5; thr_t = 5;
kin_trna = λ*trna/((g_max*atp/(θtlr+atp)) * trna/(thr_t+trna)) 

