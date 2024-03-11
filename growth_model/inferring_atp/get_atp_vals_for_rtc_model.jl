using DifferentialEquations, DataFrames, Measures, CSV, DiffEqCallbacks, LabelledArrays
using PlotlyJS

include("model.jl")
include("parameters.jl")
include("functions.jl")


tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)

# plotly_plot_sol(solu, "log", "")



ssinit = get_ss_init(solu)
ns = 0.65;
pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, sf]

solu_ss = simple_solve!(pop_model!, ssinit, tspan, pop_params)
# plotly_plot_sol_pop(solu_ss, "log", "")

# a1 = get_curve(solu_ss, :a)
# SF = 1e6/(6.022e23*1e-15)
# a = a1*SF
# p_atp = plot(scatter(x=solu_ss.t, y=a), Layout(xaxis_type="log"))

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

# species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :s, :N]
function get_maxval(sol, species)
    M = 10^8
    cr = get_curve(sol, :cr); cp = get_curve(sol, :cp); cq = get_curve(sol, :cq); ct = get_curve(sol, :ct); cm = get_curve(sol, :cm);
    Rt = @. cr+cp+cq+ct+cm
    g_max = 1260
    k_g = 7*sf
    a1 = get_curve(sol, :a)
    SF = 1e6/(6.022e23*1e-15)
    a = a1*SF
    lam = @. (Rt*(g_max*a/(k_g+a)))/(M)
    lam_max = maximum(lam)
    a_max = maximum(a)
    if species == "lam"
        return lam_max
    else
        return a_max
    end
end


param_range1 = 10 .^(range(-3, stop=0.1, length=100))
param_range2 = 10 .^(range(2,stop=10,length=100))

function sweep_params(param_range1, param_range2)
    res_lam = []
    res_a = []
    params = @LArray [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, sf] (:dm, :kb, :ku, :f, :thetar, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :wp, :nq, :nr, :ns, :kin, :d_s, :d_n, :sf)

    for i in param_range1
        params[:ns] = i
        res1 = []
        res2 = []
        for val in param_range2 
            ssinit[:s_0] = val
            solu = simple_solve!(pop_model!, ssinit, tspan, params)
            push!(res1, get_maxval(solu, "lam"))
            push!(res2, get_maxval(solu, "a"))
        end
        push!(res_lam, res1)
        push!(res_a, res2)
    end

    vec_lam = []
    for i in (1:length(param_range1))
        append!(vec_lam, values(res_lam[i]))
    end

    vec_a = []
    for i in (1:length(param_range1))
        append!(vec_a, values(res_a[i]))
    end

    vec_lam = reshape(vec_lam, (length(param_range1),length(param_range1)))
    vec_a = reshape(vec_a, (length(param_range1),length(param_range1)))
    p1 = plot(contour(x=param_range1, y=param_range2, z=vec_lam, colorbar=attr(title="λ", titleside="right", len=0.5, yanchor="bottom")), Layout(xaxis_title="ns", yaxis_title="s_0"))
    p2 = plot(contour(x=param_range1, y=param_range2, z=vec_a, colorbar=attr(title="ATP", titleside="right", len=0.5, yanchor="top")), Layout(xaxis_title="ns", yaxis_title="s_0"))
    return [p1; p2]
end

sweep_params(param_range1, param_range2)

# plot(scatter(x=lam, y=a), Layout(yaxis_type="log", yaxis_title="ATP (μM)", xaxis_title="λ (/min)"))


# plot(scatter(x=lam, y=lam_vs_a), Layout(yaxis_type="log"))


df = DataFrame(gr = lam, atp=a, time = solu_ss.t)
CSV.write("interpolation_data.csv", df)


using Dierckx

interp_func = Spline2D(lam, solu_ss.t, a, kx=5, ky=5, s=1e10)
#check Interpolation
atp_check = []
for i in 1:length(solu_ss.t)
    push!(atp_check, interp_func(lam[i], solu_ss.t[i]))
end
atp_test = plot(scatter(x=solu_ss.t, y=atp_check), Layout(xaxis_type="log"))
# atp_check[atp_check.< 0] .= 1e-90 
# plot(scatter(x=solu_ss.t, y=atp_check), Layout(xaxis_type="log"))

[p_atp; atp_test]




#load other growth rate 
csv_lam = DataFrame(CSV.File("$PATHdata/results_colD_grfit.csv"))
csv_lam = select!(csv_lam, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_lam."gr"[csv_lam."gr".< 0] .= 0 #zero(eltype(lam_colD))
plot(scatter(x=csv_lam."t", y=csv_lam."gr"))

atp_from_gr = []
for i in 1:length(csv_lam."gr")
    push!(atp_from_gr, interp_func(csv_lam."gr"[i], csv_lam."t"[i]))
end

atp_from_gr[atp_from_gr.< 0] .= 1e-90 #zero(eltype(lam_colD))

plot(scatter(x=csv_lam."t", y=atp_from_gr))


# # for the actual data growth rate find the value of atp to use in rtc model 
# atp_gr = []
# for i in csv_lam."gr"
#     push!(atp_gr, lam_vs_a(i))
# end
# atp_gr

# csv_lam."gr"[1]
# lam_vs_a(csv_lam."gr"[1])

# plot(scatter(x=solu_ss.t, y=atp_gr))#, Layout(yaxis_type="log", xaxis_type="log"))

#make a dataframe with the atp values against time (for data) - now in rtc model script need to interpolate this so it can be used over time 
df_atp = DataFrame(t=csv_lam."t", atp=atp_gr)
CSV.write("$PATHdata/atp_curve_from_growth_model.csv", df_atp)







zeros(length(lam), length(solu_ss.t))



plot(contour(y=lam, x=solu_ss.t, z=a, colorbar=attr(title="ATP")), Layout(xaxis_type="log", yaxis_title="λ", xaxis_title="t"))


scatterplot = scatter3d(
    x=a,
    y=lam,
    z=solu_ss.t,
    mode="markers",
    marker_color=lam
)



layout=Layout(scene=attr(zaxis=attr( range=[z[1],650]), xaxis_title="ATP", yaxis_title="λ", zaxis_title="t"))
# Create the plot
plot(scatterplot, layout)


using Interpolations

itp_a = LinearInterpolation(solu_ss.t, a)

a_at_growth_rate = [itp_a(t) for t in solu_ss.t]

combined_data = hcat(lam, solu_ss.t, a_at_growth_rate)



