using CSV, DataFrames, PlotlyJS

include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


v = 80/1000000
g = 0.0045*200
g1 = (0.0045*1000000) #*(200*10^-6)
k = (138*10^6)#*1000000
k1 = 33.4 

ONPG_mgml = 1.1 # or 4
ONPG_mw = 301.251
ONPG = (ONPG_mgml/ONPG_mw)#*10^3

SF = 1e6/(6.022e23*1e-15)

WT_mRNA_conc_mu = []
for (i,j) in zip(csv[!,2], WT1)
    push!(WT_mRNA_conc_mu, ((i*0.5/20)*SF))
end
WT_mRNA_conc_mu


WT_mRNA_conc = []
for (i,j) in zip(csv[!,2], WT1)
    push!(WT_mRNA_conc, ((i*v*j/(1000*g1*k1*ONPG))/4)/20)
end
WT_mRNA_conc
WT_mRNA_conc.*1000000


WT_mRNA_conc = []
for (i,j) in zip(csv[!,2], WT1)
    push!(WT_mRNA_conc, (i*v*j/(1000*g*k*ONPG))/20)
end
WT_mRNA_conc

WT_mRNA_conc_std = []
for (i,j) in zip(csv[!,2], WT1_std)
    push!(WT_mRNA_conc_std, (i*v*j/(1000*g*k*ONPG))/20)
end
WT_mRNA_conc_std