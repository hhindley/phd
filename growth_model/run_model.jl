using DifferentialEquations, PlotlyJS, DataFrames, Measures, LabelledArrays, BenchmarkTools, ModelingToolkit, TickTock

PATH = "/home/holliehindley/phd"

include("$PATH/growth_model/model/model.jl")
include("$PATH/general_funcs/solving.jl")
include("$PATH/growth_model/parameters/parameters.jl")
include("$PATH/growth_model/parameters/uM_parameters.jl")

tspan = (0,1e9)

labels = ["cr" "em" "cp" "cq" "ct" "et" "cm" "mt" "mm" "q" "p" "si" "mq" "mp" "mr" "r" "a"]

# params_uM = [0, 0.1, 0.00166, 0.00166, 0, 15594.7, 0, 1661, 2.0923, 0, 160.01, 1.66, 1e8, 0.00687, 1.66, 5800, 300, 252.77, 0.299, 726, 1.54, 1.576, 0.00687, 4, 7459, 0.5, 255.7]

init_abx=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 1000, 0, 0, 0, 0]
params_abx = [dm_val, kb_val, ku_val, thetar_val, s0_val, gmax_val, thetax_val, Kt_val, M_val, we_val, Km_val, vm_val, nx_val, Kq_val, vt_val, wr_val, wq_val, nq_val, nr_val, ns_val, Kgamma_val, Cm_val, k_cm_val]

prob = ODEProblem(growth_model_incl_abx, init_abx, tspan, params_abx);
solu = solve(prob, Rodas4(), abstol=1e-12, reltol=1e-9);
df_gm = create_solu_df(solu, species_gm);
p = plot([scatter(x=df_gm.time, y=col, name="$(names(df_gm)[i])") for (col, i) in zip(eachcol(df_gm[:,2:end]), range(2,length(names(df_gm))))], Layout(xaxis_type="log"))

prob2 = ODEProblem(growth_model, init_gm, tspan, params_gm; jac=true);
solu2 = solve(prob2, Rodas4(), abstol=1e-12, reltol=1e-9);
df = create_solu_df(solu2, species_gm);
p = plot([scatter(x=df.time, y=col, name="$(names(df)[i])") for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log"))


# tick(); solve(prob, Rodas4(), abstol=1e-12, reltol=1e-9); tock() 

# tick(); solve(prob2, Rodas4(), abstol=1e-12, reltol=1e-9); tock() 
# solu_gm = simple_solve!(growth_model_incl_abx, init_abx, tspan, params_abx)

# gm_species = [:cr, :em, :cq, :ct, :et, :cm, :mt, :mm, :q, :si, :mq, :mr, :r, :a]
# gm_species_abx = [:cr, :em, :cq, :ct, :et, :cm, :mt, :mm, :q, :si, :mq, :mr, :r, :a, :zmr, :zmq, :zmt, :zmm]

# solu_gm = solve(prob, Rodas4(), abstol=1e-12, reltol=1e-9);
df_gm = create_solu_df(solu, species_gm);
p = plot([scatter(x=df_gm.time, y=col, name="$(names(df_gm)[i])") for (col, i) in zip(eachcol(df_gm[:,2:end]), range(2,length(names(df_gm))))], Layout(xaxis_type="log"))

df_gm_uM = create_solu_df(solu_gm_uM, species_gm);
p = plot([scatter(x=df_gm_uM.time, y=col, name="$(names(df_gm_uM)[i])") for (col, i) in zip(eachcol(df_gm_uM[:,2:end]), range(2,length(names(df_gm_uM))))], Layout(xaxis_type="log"))




gamma = @. gmax*df_gm.a/(Kgamma+df_gm.a) # aa min-1
ttrate = @. (df_gm.cq+df_gm.cr+df_gm.ct+df_gm.cm)*gamma 
lam_gm = @. ttrate/M

params1 = deepcopy(params)
nut = range(0,1e4, length=20)
gr = []
for i in nut
    params1[6] = i
    solu_gm = simple_solve!(growth_model, init_gm, tspan, params1)
    df_gm = create_solu_df(solu_gm, gm_species)
    gamma = gmax*df_gm.a[end]/(Kgamma+df_gm.a[end]) # aa min-1
    ttrate = (df_gm.cq[end]+df_gm.cr[end]+df_gm.ct[end]+df_gm.cm[end])*gamma 
    lam_gm = ttrate/M
    push!(gr, lam_gm)
end
# open("$PATHgrowth_model_2021/growth_model/result_plot.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end

plot(scatter(x=nut, y=gr))

solu_gm_uM = simple_solve!(growth_model, init_gm_uM, tspan, params_uM)
df_gm_uM = create_solu_df(solu_gm_uM, gm_species)

p = plot([scatter(x=df_gm_uM.time, y=col, name="$(names(df_gm_uM)[i])") for (col, i) in zip(eachcol(df_gm_uM[:,2:end]), range(2,length(names(df_gm_uM))))], Layout(xaxis_type="log"))

gamma_uM = @. gmax*df_gm_uM.a/(Kgamma_uM+df_gm_uM.a) # aa min-1
ttrate_uM = @. (df_gm_uM.cq+df_gm_uM.cr+df_gm_uM.ct+df_gm_uM.cm)*gamma_uM 
lam_uM = @. ttrate_uM/M_uM

plot([scatter(x=df_gm.time, y=lam_gm, name="molecs"),scatter(x=df_gm_uM.time, y=lam_uM, name="uM")], Layout(xaxis_type="log"))







species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
a = solDF[!, :a]
cq = solDF[!, :cq]
cr = solDF[!, :cr]
cp = solDF[!, :cp]
ct = solDF[!, :ct]
cm = solDF[!, :cm]
si = solDF[!, :si]
em = solDF[!, :em]
et = solDF[!, :et]
mt = solDF[!, :mt]
mm = solDF[!, :mm]
mq = solDF[!, :mq]
mp = solDF[!, :mp]
mr = solDF[!, :mr]
r = solDF[!, :r]
p = solDF[!, :p]
q = solDF[!, :q]

plot(sol.t, r)

plot(sol.t, cp)