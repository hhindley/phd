SF = 1e6/(6.022e23*1e-15)

# ribosomes
ribo_cell =  [6800,13500,26300,45100,72000]
ribo_conc = []
for i in ribo_cell
    push!(ribo_conc, i*SF)
end
print(ribo_conc)

#checking gr_c values from real data 
grs = [0.6,1,1.5,2,2.5]/60
print(grs)
cs = []
for (i,j) in zip(grs, ribo_conc)
    push!(cs, i/j)
end
print(cs)
growth_rate_constant = mean(cs)

# growth rate constant
lam = 2.77/60 # taken from colD data in google colab notebook - need more WT data to base this off I think
lambda = grs[end]
rh =  ribo_conc[end] # taken from ribosome table at the highest growth rate 
gr_c = lam/rh

#initial value of rh 
growth_rate = 2.5/60 # /60 to get from per hour to per minute
rh_0 = growth_rate/growth_rate_constant

# influx of healthy ribosomes
kin = lambda*rh/g_max