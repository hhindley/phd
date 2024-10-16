using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes, KernelDensity, Arrow

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_switch_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

@load "/Users/s2257179/Desktop/saved_variables/high_kdam_stops.jld2" df_lengths df_stops 

df_rtca = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/high_kdam_rtca.arrow"))
df_times = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/high_kdam_times.arrow"))

kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

res = Dict()
times_res = Dict()
for i in eachindex(names(df_rtca))
    res[names(df_rtca)[i]] = collect(skipmissing(df_rtca[!,names(df_rtca)[i]]))
    times_res[names(df_times)[i]] = collect(skipmissing(df_times[!,names(df_times)[i]]))
end
res

# split into original simulations
sims = Dict(i=>Dict(kdam=>[] for kdam in kdams) for i in 1:20)
times = Dict(i=>Dict(kdam=>[] for kdam in kdams) for i in 1:20)
for kdam in eachindex(kdams)
    sims[1][kdams[kdam]] = res["$(kdams[kdam])"][1:df_stops[1,"$(kdams[kdam])"]]
    times[1][kdams[kdam]] = times_res["$(kdams[kdam])"][1:df_stops[1,"$(kdams[kdam])"]]
    for j in eachindex(kdams)
        if j < 20
            ind = j+1
        else 
            continue
        end
        sims[ind][kdams[kdam]] = res["$(kdams[kdam])"][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
        times[ind][kdams[kdam]] = times_res["$(kdams[kdam])"][df_stops[ind-1,"$(kdams[kdam])"]+1:df_stops[ind,"$(kdams[kdam])"]]
    end
end
sims[9]

lines(times[1][0.0], sims[1][0.0])
lines(times[2][0.0], sims[2][0.0])

f = Figure()
ax = Axis(f[1,1], yscale=log10)
[hist!(ax, res[i]) for i in eachindex(res)]

function waterfall_makie(x, y, z; zmin = minimum(z), lw = 1., colmap = :linear_bgy_10_95_c74_n256, colorband = (:white, 1.), xlab = "x", ylab = "y", zlab = "z")
    # Initialisation
    fig = Figure()
    ax = Axis3(fig[1,1], xlabel = xlab, ylabel = ylab, zlabel = zlab)
    for (j, yv) in enumerate(y)
        
        zj = z[j, :]
        println(length(zj))
        println(length(yv))
        lower = Point3f.(x, yv, zmin)
        upper = Point3f.(x, yv, zj)
        # edge_start = [Point3f(x[1], yv, zmin), Point3f(x[1], yv, zj[1])]
        # edge_end = [Point3f(x[end], yv, zmin), Point3f(x[end], yv, zj[end])]

        # # Surface
        band!(ax, lower, upper, color = colorband)

        # Line
        lines!(ax, upper, color = zj, colormap = colmap, linewidth = lw)
        
        # Edges
        # lines!(ax, edge_start, color = zj[1]*ones(2), colormap = colmap, linewidth = lw)
        # lines!(ax, edge_end, color = zj[end]*ones(2), colormap = colmap, linewidth = lw)
    end

    # Set axes limits
    xlims!(ax, minimum(x), maximum(x))
    ylims!(ax, minimum(y), maximum(y))
    zlims!(ax, zmin, maximum(z))

    fig
end

# Generate 5 different datasets
data1 = randn(100) .+ 2  # Normal distribution, mean = 2
data2 = randn(100) .+ 0  # Normal distribution, mean = 0
data3 = randn(100) .* 2  # Normal distribution, larger variance
data4 = 3 .+ rand(100) * 2  # Uniform distribution, between 3 and 5
data5 = exp.(randn(100) .* 0.5)  # Log-normal distribution

datasets = [data1, data2, data3, data4, data5]
kdes = [kde(data) for data in datasets]
x = kdes[1].x
y=1:length(datasets)
z = hcat([kde.density for kde in kdes]...)
z = permutedims(z, (2, 1)) 
nx=length(x)
ny=length(y)

fig = waterfall_makie(x, y, z, xlab="conc", ylab="freq", zlab="dam")


kdes_res = [kde(res[i]) for i in eachindex(res)]
kdes_res = [kde(log.(res[i] .+ 1)) for i in eachindex(res)]  # Apply log transformation

x = kdes_res[1].x
y=1:length(res)
z = hcat([kde.density for kde in kdes_res]...)
z = permutedims(z, (2, 1))
nx = length(x)
ny = length(y)

fig = waterfall_makie(x, y, z, xlab="conc", ylab="dam", zlab="freq")


hist(res[1], bins=20)
[lines!(kdes_res[i].x, kdes_res[i].density) for i in eachindex(kdes_res)]