using CSV, DataFrames, PlotlyJS

include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


v = 80/1000000
# g = (0.0045/1000000)*0.0002
g1 = (0.0045/1000000) #*(200*10^-6)
k = (138*10^6)#*1000000

ONPG_mgml = 4 # or 4
ONPG_mw = 301.251
ONPG = (ONPG_mgml/ONPG_mw)#*10^3


WT_mRNA_conc_M = []
for (i,j) in zip(WT1, WT2)
    push!(WT_mRNA_conc_M, (i*v*j/(1000*g1*k*(ONPG)))/80)
end
WT_mRNA_conc_M

mRNA_conc_mM = @. WT_mRNA_conc_M*1000

mRNA_conc_uM = @. WT_mRNA_conc_M*1000000

conc_std = []
for (i,j) in zip(WT1_std, WT2_std)
    push!(conc_std, (i*v*j/(1000*g1*k*(ONPG)))/80)
end
conc_std = @. conc_std*1000000





p = Subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)

p1 = plot(scatter(x=dfc[:,1], y=WT1, showlegend=false), Layout(yaxis_title="Miller Units"))
p2 = plot(scatter(x=dfc[:,1], y=WT2, showlegend=false), Layout(yaxis_title="OD600"))

p3 = plot(scatter(x=dfc[:,1], y=mRNA_conc_uM, showlegend=false), Layout(xaxis_title="time (hours)", yaxis_title="[LacZ mRNA] (μM)"))

p = [p1; p2; p3]




WT1
WT1_std

(WT1[1]-WT1_std[1])/WT1[1]
WT1_std[1]/WT1[1]
mRNA_conc_uM[1]*0.214

stds=[]
for (i,j) in zip(WT1, WT1_std)
    print(i,j)
    push!(stds, j/i)
end 
stds

new_stds=[]
for (i,j) in zip(mRNA_conc_uM, stds)
    push!(new_stds, i*j)
end 
new_stds

conc_df = DataFrame(mRNA_conc_uM=mRNA_conc_uM, conc_std=new_stds)

CSV.write("/home/holliehindley/phd/data/mRNA_conc_uM.csv", string.(conc_df))
