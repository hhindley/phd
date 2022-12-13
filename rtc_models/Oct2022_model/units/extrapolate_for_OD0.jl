# extrapolate to get value of OD_0
using Plots, Dierckx
time1 = dfe[:,1]
spl = Spline1D(time1, WT2; k=1, bc="extrapolate")
OD_0 = evaluate(spl, 0)

time2 = dff[:,1]
spl2 = Spline1D(time2, WT3; k=1, bc="extrapolate")
OD_02 = evaluate(spl2, 0)

time3 = df2[:,1]
spl3 = Spline1D(time3, WT4; k=1, bc="extrapolate")
OD_03 = evaluate(spl3, 0)

spl4 = Spline1D(time2, hpx_rtcon; k=1, bc="extrapolate")
OD_0_hpxon = evaluate(spl4, 0)

spl5 = Spline1D(time3, WT_colD; k=1, bc="extrapolate")
OD_0_wtcolD = evaluate(spl5, 0)

spl6 = Spline1D(time3, nA_colD; k=1, bc="extrapolate")
OD_0_nAcolD = evaluate(spl6, 0)

spl7 = Spline1D(time3, nB_colD; k=1, bc="extrapolate")
OD_0_nBcolD = evaluate(spl7, 0)

spl8 = Spline1D(time3, nB_B_colD; k=1, bc="extrapolate")
OD_0_nBBcolD = evaluate(spl8, 0)

spl9 = Spline1D(time3, nB_Bmut_colD; k=1, bc="extrapolate")
OD_0_nBBmutcolD = evaluate(spl9, 0)

spl10 = Spline1D(time3, nR_colD; k=1, bc="extrapolate")
OD_0_nRcolD = evaluate(spl10, 0)


function check_OD_0(time, data, OD)
    p1 = plot(time*60, data, markershape=:circle)
    scatter!((0,OD));
    time2 = [[0];time]
    data2 = [[OD];data]
    p2 = plot(time2*60, data2, markershape=:circle)
    return display(plot(p1, p2, layout=(2,1)))
end

check_OD_0(time1, WT2, OD_0)

check_OD_0(time2, WT3, OD_02)

check_OD_0(time3, WT4, OD_03)

check_OD_0(time2, hpx_rtcon, OD_0_hpxon)

check_OD_0(time3, WT_colD, OD_0_wtcolD)

check_OD_0(time3, nA_colD, OD_0_nAcolD)

check_OD_0(time3, nB_colD, OD_0_nBcolD)

check_OD_0(time3, nB_B_colD, OD_0_nBBcolD)

check_OD_0(time3, nB_Bmut_colD, OD_0_nBBmutcolD)

check_OD_0(time3, nR_colD, OD_0_nRcolD)
