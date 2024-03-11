using CSV, DataFrames, PlotlyJS


rpon = DataFrame(CSV.File("$PATHdata/rpon_means.csv")) # read csv to a datafram
rpon_stds = DataFrame(CSV.File("$PATHdata/rpon_stds.csv")) # read csv to a datafram


plot(rpon, x=:time, y=:cm1, error_y=attr(type="data",array=rpon_stds.cm1))


cm1 = scatter(x=rpon.time, y=rpon.cm1, error_y=(attr(type="data",array=rpon_stds.cm1)), name="Cm 1.25")
cm2 = scatter(x=rpon.time, y=rpon.cm2, error_y=(attr(type="data",array=rpon_stds.cm2)), name="Cm 2.5")
cm3 = scatter(x=rpon.time, y=rpon.cm3, error_y=(attr(type="data",array=rpon_stds.cm3)), name="Cm 5")
cm10 = scatter(x=rpon.time, y=rpon.cm10, error_y=(attr(type="data",array=rpon_stds.cm10)), name="Cm 10")

tet1 = scatter(x=rpon.time, y=rpon.tet1, error_y=(attr(type="data",array=rpon_stds.tet1)), name="Tet 0.75")
tet2 = scatter(x=rpon.time, y=rpon.tet2, error_y=(attr(type="data",array=rpon_stds.tet2)), name="Tet 1.5")
tet4 = scatter(x=rpon.time, y=rpon.tet4, error_y=(attr(type="data",array=rpon_stds.tet4)), name="Tet 4")

wt = scatter(x=rpon.time, y=rpon.wt, error_y=(attr(type="data",array=rpon_stds.wt)), name="WT")

p = plot([cm1,cm2,cm3,cm10,tet1,tet2,tet4,wt], Layout(title="rpoN expression", xaxis_title="Time (hours)", yaxis_title="Copy number per cell", yaxis_type="log"))

open("$PATHdata/rpon_log.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end






