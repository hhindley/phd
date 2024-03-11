using CSV, DataFrames, PlotlyJS

include("$PATHParam_inf/inf_setup.jl")


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


hpx_mRNA_conc_M = []
for (i,j) in zip(hpx, hpx_rtcoff)
    push!(hpx_mRNA_conc_M, (i*v*j/(1000*g1*k*(ONPG)))/80)
end
hpx_mRNA_conc_M

hpx_mRNA_conc_uM = @. hpx_mRNA_conc_M*1000000

((hpx[end]*v*hpx_rtcoff[end]/(1000*g1*k*(ONPG)))/80)*1000000
((WT1[end]*v*WT2[end]/(1000*g1*k*(ONPG)))/80)*1000000

# conc_std = []
# for (i,j) in zip(WT1_std, WT2_std)
#     push!(conc_std, (i*v*j/(1000*g1*k*(ONPG)))/80)
# end
# conc_std = @. conc_std*1000000



WT1
hpx

p1 = plot(scatter(x=dfc[:,1], y=WT1, showlegend=false), Layout(yaxis_title="Miller Units"));#, title="WT"));
p2 = plot(scatter(x=dfc[:,1], y=WT2, showlegend=false), Layout(yaxis_title="OD600"));#, title="WT"));

p3 = plot(scatter(x=dfc[:,1], y=mRNA_conc_uM, showlegend=false), Layout(yaxis_range=[0,190], xaxis_title="time (hours)", yaxis_title="WT [LacZ mRNA] (μM)"));#, title="WT"));

p4 = plot(scatter(x=dfc[:,1], y=hpx, showlegend=false), Layout(yaxis_title="Miller Units"));#, title="Hpx-"));
p5 = plot(scatter(x=dfc[:,1], y=hpx_rtcoff, showlegend=false), Layout(yaxis_title="OD600"));#, title="Hpx-"));

p6 = plot(scatter(x=dfc[:,1], y=hpx_mRNA_conc_uM, showlegend=false), Layout(yaxis_range=[0,190],xaxis_title="time (hours)", yaxis_title="Hpx- [LacZ mRNA] (μM)"));#, title="Hpx-"));
[p1;p4]
p = [p1; p2; p3]
p_hpx = [p4;p5;p6]
p_both = [p3;p6]


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

# CSV.write("$PATHdata/mRNA_conc_uM.csv", string.(conc_df))



paper_conv = []
for i in WT1
    push!(paper_conv, ((i*0.5)/(6.022e23*1e-15))/80)
end
paper_conc_mM = @. paper_conv*1000

paper_conc = @. paper_conv*1000000

paper_stds=[]

for (i,j) in zip(paper_conc, stds)
    push!(paper_stds, i*j)
end
paper_stds


p7 = plot(scatter(x=dfc[:,1], y=paper_conc, showlegend=false), Layout(xaxis_title="time (hours)", yaxis_title="[LacZ mRNA] (μM) (paper)"))#, title="Hpx-"))

[p3;p6;p7]

paper_df = DataFrame(mRNA_conc_uM=paper_conc, conc_std=paper_stds)
CSV.write("$PATHdata/paper_conv.csv", string.(paper_df))
