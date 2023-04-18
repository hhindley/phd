using DifferentialEquations, PyCall, DataFrames, Measures, CSV, DiffEqCallbacks, DataInterpolations, LabelledArrays, Statistics
using PlotlyJS


include("model.jl")
include("parameters.jl")
include("functions.jl")


tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)


ssinit = get_ss_init(solu)
ns = 5.5; #0.85 
pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, sf]

solu_ss = simple_solve!(pop_model!, ssinit, tspan, pop_params)
plotly_plot_sol_pop(solu_ss, "log", "")

M = 10^8
cr = get_curve(solu_ss, :cr); cp = get_curve(solu_ss, :cp); cq = get_curve(solu_ss, :cq); ct = get_curve(solu_ss, :ct); cm = get_curve(solu_ss, :cm);
Rt = @. cr+cp+cq+ct+cm
g_max = 1260
k_g = 7*sf
a1 = get_curve(solu_ss, :a)
SF = 1e6/(6.022e23*1e-15)
a = a1*SF

lam = @. (Rt*(g_max*a/(k_g+a)))/(M)
# p_atp1 = plot(scatter(x=solu_ss.t, y=a1), Layout(xaxis_type="log"))
p_atp = plot(scatter(x=solu_ss.t, y=a), Layout(xaxis_type="log", title="ATP"))
p_lam = plot(scatter(x=solu_ss.t, y=lam), Layout(xaxis_type="log", title="λ"))


open("./atp_growth_model.html", "w") do io
    PlotlyBase.to_html(io, p_atp.plot)
end

# # species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :s, :N]
# function get_maxval(sol, species)
#     M = 10^8
#     cr = get_curve(sol, :cr); cp = get_curve(sol, :cp); cq = get_curve(sol, :cq); ct = get_curve(sol, :ct); cm = get_curve(sol, :cm);
#     Rt = @. cr+cp+cq+ct+cm
#     g_max = 1260
#     k_g = 7*sf
#     a1 = get_curve(sol, :a)
#     SF = 1e6/(6.022e23*1e-15)
#     a = a1*SF
#     lam = @. (Rt*(g_max*a/(k_g+a)))/(M)
#     lam_max = maximum(lam)
#     a_max = maximum(a)
#     if species == "lam"
#         return lam_max
#     else
#         return a_max
#     end
# end


# param_range1 = 10 .^(range(-3, stop=0.1, length=100))
# param_range2 = 10 .^(range(2,stop=10,length=100))

# function sweep_params(param_range1, param_range2)
#     res_lam = []
#     res_a = []
#     params = @LArray [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, sf] (:dm, :kb, :ku, :f, :thetar, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :wp, :nq, :nr, :ns, :kin, :d_s, :d_n, :sf)

#     for i in param_range1
#         params[:ns] = i
#         res1 = []
#         res2 = []
#         for val in param_range2 
#             ssinit[:s_0] = val
#             solu = simple_solve!(pop_model!, ssinit, tspan, params)
#             push!(res1, get_maxval(solu, "lam"))
#             push!(res2, get_maxval(solu, "a"))
#         end
#         push!(res_lam, res1)
#         push!(res_a, res2)
#     end

#     vec_lam = []
#     for i in (1:length(param_range1))
#         append!(vec_lam, values(res_lam[i]))
#     end

#     vec_a = []
#     for i in (1:length(param_range1))
#         append!(vec_a, values(res_a[i]))
#     end

#     vec_lam = reshape(vec_lam, (length(param_range1),length(param_range1)))
#     vec_a = reshape(vec_a, (length(param_range1),length(param_range1)))
#     p1 = plot(contour(x=param_range1, y=param_range2, z=vec_lam, colorbar=attr(title="λ", titleside="right", len=0.5, yanchor="bottom")), Layout(xaxis_title="ns", yaxis_title="s_0"))
#     p2 = plot(contour(x=param_range1, y=param_range2, z=vec_a, colorbar=attr(title="ATP", titleside="right", len=0.5, yanchor="top")), Layout(xaxis_title="ns", yaxis_title="s_0"))
#     return [p1; p2]
# end

# sweep_params(param_range1, param_range2)


# df = DataFrame(gr = lam, atp=a, time = solu_ss.t)
# CSV.write("interpolation_data.csv", df)

np = pyimport("numpy")
interp = pyimport("scipy.interpolate")

points = np.vstack((solu_ss.t, lam))'

atp_interp = interp.NearestNDInterpolator(points, a)

atp_test = atp_interp(solu_ss.t, lam)

plot(scatter(x=solu_ss.t, y=atp_test), Layout(xaxis_type="log"))




#load other growth rate 
csv_lam = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv"))
csv_lam = select!(csv_lam, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_lam."gr"[csv_lam."gr".< 0] .= 0 #zero(eltype(lam_colD))
plot(scatter(x=csv_lam."t", y=csv_lam."gr"))


function extend_gr_curve(csv)
    mean_gr = mean((csv."gr"[Int64((length(csv."t")*2/3)+1):end]))
    df = DataFrame(t=Float64[], gr=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, mean_gr])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."gr",new_df."t")
    return lam, new_df
end

lam_ext, new_df = extend_gr_curve(csv_lam)

lam_ext[lam_ext.< 0] .= 0 #zero(eltype(lam_colD))
plot(scatter(x=new_df."t", y=lam_ext), Layout(xaxis_range=(0,1000)))

plot(scatter(x=new_df."t", y=new_df."gr"), Layout(xaxis_type="log"))



atp_rtc = atp_interp(new_df."t", new_df."gr")

plot(scatter(x=new_df."t", y=atp_rtc), Layout(xaxis_type="log"))

# CSV.write("colD_gr_data.csv", new_df)

p1 = plot(scatter(x=new_df."t", y=new_df."gr"), Layout(xaxis_range=(0,1000), title="λ_data"))
p2 = plot(scatter(x=solu_ss.t, y=lam), Layout(xaxis_range=(0,1000), title="λ_growth_model"))

[p1;p2]


p = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
add_trace!(p, (scatter(x=new_df."t", y=new_df."gr")), row=1, col=1)
add_trace!(p, (scatter(x=solu_ss.t, y=lam)), row=2, col=1)
relayout!(p, showlegend=false, xaxis_type="", xaxis2_type="")
p

data_t = new_df."t"[2:end]
data_lam = new_df."gr"[2:end]
start_t = solu_ss.t[1:50]
start_lam = lam[1:50]
addon = data_t[1] - start_t[end] -1
start_t_moved = start_t .+ addon
beg_range = range(start_t[3], start_t_moved[3], length=50)
beg_lam = repeat([mean(start_lam[1:3])], outer = 50)

full_t = [start_t; data_t]
full_lam = [start_lam; data_lam]



plot(scatter(x=full_t, y=full_lam), Layout(xaxis_type="log"))

atp_rtc = atp_interp(full_t, full_lam)

p = make_subplots(rows=2, cols=2, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["Model λ" "Model ATP"; "Data λ" "Inferred ATP"])
add_trace!(p, (scatter(x=solu_ss.t, y=lam)), row=1, col=1)
add_trace!(p, (scatter(x=solu_ss.t, y=a)), row=2, col=1)
add_trace!(p, (scatter(x=full_t, y=full_lam)), row=1, col=2)
add_trace!(p, (scatter(x=full_t, y=atp_rtc)), row=2, col=2)
relayout!(p, showlegend=false, xaxis_type="log", xaxis2_type="log", xaxis3_type="log", xaxis4_type="log")
p

atp_lam_data = DataFrame(t = full_t, atp = atp_rtc, lam = full_lam)

CSV.write("/home/holliehindley/phd/data/atp_for_rtcmodel.csv", atp_lam_data)

