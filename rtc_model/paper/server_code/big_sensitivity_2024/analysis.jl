using StatsPlots, Distributions, Statistics, LinearAlgebra, DataFrames, JLD2, KernelDensity, Combinatorics, Measures

pb = load("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/params_nokin.jld2")
pn = load("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/params_notbs_nokin.jld2")

df1 = DataFrame(pb)
df2 = DataFrame(pn)

df_b = DataFrame(pb)
df_n = DataFrame(pn)
rename!(df_b, "Km_init"=>:Km)
rename!(df_b, "Vmax_init"=>:Vmax)
rename!(df_n, :Km_init=>:Km)
rename!(df_n, :Vmax_init=>:Vmax)

df_b = transform(df_b, :Km=>(x->"bistable")=>:type)
df_n = transform(df_n, :Km=>(x->"not_bistable")=>:type)

df = vcat(df_b, df_n)

df = stack(df, Not(:type))

df_b = stack(df_b, Not(:type))
df_n = stack(df_n, Not(:type))

theme(:default, 
    background_color = "white", 
    fontfamily="helvetica", 
    grid=false, 
    guidefontsize=14, 
    legendfontsize=14, 
    tickfontsize=14, 
    legend=false)

p = @df df_b violin(:variable, :value, group=:variable, color=:red, yaxis=:log, side=:left, alpha=0.6, label="bistable", legend=:false, yticks=10 .^(-2.0:1:5))
@df df_n violin!(:variable, :value, group=:variable, color=:blue, side=:right, yaxis=:log, alpha=0.6, label="not bistable", legend=:false)

savefig(p, "/home/hollie_hindley/Documents/paper/big_sensitivity_2024/violin_plot.svg")


dens_b=[]; dens_n=[]; ppairs=[];
for pair in combinations(names(df1),2)
    col1, col2 = pair
    push!(dens_b, kde((df1[:,col1], df1[:,col2])))
    push!(dens_n, kde((df2[:,col1], df2[:,col2])))
    push!(ppairs, pair)
end

sps=[]
for i in range(1,length(dens_b))
    p = plot(dens_b[i], xticks=false, yticks=false, size=(400,400), color=:red,alpha=0.5)#, xlabel=ppairs[i][1], ylabel=ppairs[i][2], guidefontsize=9)#, alpha=0.8)#, xlabel=ppairs[i][1], ylabel=ppairs[i][2])
    if ppairs[i][1] == "kdeg" || ppairs[i][2] == "kdeg"
        plot!(p, dens_n[i], levels=140, xticks=false, yticks=false, size=(400,400), color=:blue,alpha=0.5)#,xlabel=ppairs[i][1], ylabel=ppairs[i][2], guidefontsize=9)#, alpha=0.8)
    else
        plot!(p, dens_n[i], xticks=false, yticks=false, size=(400,400), color=:blue,alpha=0.5)#,xlabel=ppairs[i][1], ylabel=ppairs[i][2], guidefontsize=9)#, alpha=0.8)
    end
    push!(sps, p)
end 

plot(sps[1])

e = plot(xticks=false, yticks=false, size=(50,50), framestyle=:none);

new_arr = [sps[1], e,       e,       e,       e,       e,       e,       e,
           sps[2], sps[9],  e,       e,       e,       e,       e,       e,
           sps[3], sps[10], sps[16], e,       e,       e,       e,       e,
           sps[4], sps[11], sps[17], sps[22], e,       e,       e,       e,
           sps[5], sps[12], sps[18], sps[23], sps[27], e,       e,       e,  
           sps[6], sps[13], sps[19], sps[24], sps[28], sps[31], e,       e,      
           sps[7], sps[14], sps[20], sps[25], sps[29], sps[32], sps[34], e,      
           sps[8], sps[15], sps[21], sps[26], sps[30], sps[33], sps[35], sps[36]]

l = @layout [grid(8,8)]

p = plot(new_arr..., layout=l, size=(500,500), margin=0mm)

savefig(p, "/home/hollie_hindley/Documents/paper/big_sensitivity_2024/contours.svg")
