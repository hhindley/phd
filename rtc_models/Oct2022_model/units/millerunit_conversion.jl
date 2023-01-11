using CSV, DataFrames, PlotlyJS

include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


v = 80
g = 0.0045*200
k = (1.38e9)*1000000

ONPG_mgml = 4
ONPG_mw = 301.251
ONPG = (ONPG_mgml/ONPG_mw)*10^6

SF = 1e6/(6.022e23*1e-15)

WT_mRNA_conc_mu = []
for (i,j) in zip(csv[!,2], WT1)
    push!(WT_mRNA_conc_mu, ((i*0.5/20)*SF))
end
WT_mRNA_conc_mu



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