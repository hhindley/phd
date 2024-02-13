using Statistics, PlotlyJS, DataInterpolations

include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

SF = 1e6/(6.022e23*1e-15)
1/SF
4000*SF

# ribosomes
ribo_cell =  [6800,13500,26300,45100,72000]
ribo_conc = []
for i in ribo_cell
    push!(ribo_conc, i*SF)
end
ribo_conc = convert(Vector{Float64}, ribo_conc)

#checking gr_c values from real data 
grs = [0.6,1,1.5,2,2.5]/60
print(grs)
cs = []
for (i,j) in zip(grs, ribo_conc)
    push!(cs, i/j)
end
print(cs)
growth_rate_constant = mean(cs)

plot(scatter(x=grs, y=ribo_conc))
int = LinearInterpolation(ribo_conc, grs)
r_c = int(0.0) # ribosome conc at growth rate of 0.013

function plot_int(fit, title, t1, gr) # must be plotted with plots.jl
    p1 = Plots.scatter(t1, gr, label="inut pdata",title=title, xlabel="Growth rate (/min)", ylabel="Ribosome conc (μM)")
    Plots.plot!(fit, label="fit")
    return display(p1)
end
plot_int(int, "Linear Interpolation", grs, ribo_conc)


# growth rate constant
lam_low = grs[1] # taken from lowest growth rate in table based off data from google colab notebook 
lam_high = grs[end]
rh_low =  ribo_conc[1] # taken from ribosome table at the lowest growth rate 
rh_high = ribo_conc[end]
gr_c_low = lam_low/rh_low # low
gr_c_high = lam_high/rh_high
gr_c_mean = mean(cs)

lam_data = 0.013
rh_data = r_c
gr_c_data = lam_data/rh_data

#initial value of rh 
rh_0_low = lam_low/gr_c_low
rh_0_high = lam_high/gr_c_high
# or using gr_c mean 
rh_0_low_mean = lam_low/gr_c_mean
rh_0_high_mean = lam_high/gr_c_mean

rh_0_data = lam_data/gr_c_data

# influx of healthy ribosomes
kin_low = lam_low*rh_low/g_max
kin_high = lam_high*rh_high/g_max

kin_data = lam_data*rh_data/g_max