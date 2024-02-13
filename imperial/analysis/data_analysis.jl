using CSV
using DataFrames
using PlotlyJS

csv = DataFrame(CSV.File("/home/holliehindley/BayesOpt/data.csv")) # read csv to a dataframe

names(csv) # get column names

data = select(csv, Not([:"Details",:"RtcA",:"RtcB",:"RtcR"])) # remove certain columns 

transform!(data, :Name => ByRow(x-> "Δ" *x[2:end]) => :Name) # remove first character in names column for each data entry 

# transform!(data, :SD => ByRow(x->x[2:end]) => :SD)
depends = ["RtcR, RtcA, RtcB", "RtcR, RtcA, RtcB", "RtcR, RtcA, RtcB", "RtcR, RtcB", "RtcR, RtcB", "No data"]
data.depends = depends
rename!(data, :"Effect (fold induction)" => :Effect) # rename column 

# data[!,:SD] = parse.(Float64, data[:, :SD]) # change column from string to Float64 so the SD can be plotted 
rep = data."Replicates"
effect = data."Effect"

# data
# color_vec = fill("SteelBlue", 6)
# color_vec[4] = "crimson"
# color_vec[5] = "crimson"
# color_vec[6] = "lightslategray"

fig = plot(data, x=:Name, y=:Effect, kind="bar", error_y=attr(type="data", array=:SD), color=:depends, labels=attr(depends="Dependent on:"),#marker_color=color_vec,
    Layout(yaxis_title_text="Effect on rtc expression (-fold induction)", xaxis_title_text="Condition", legend=attr(x=0.8, y=0.99),
    # annotations=[attr(x="ΔahpC",y=2.5, text=rep[1], showarrow=false, font=attr(color="white")), 
    # attr(x="ΔmazF", y=effect[2]/2, text=rep[2], showarrow=false, font=attr(color="white")),
    # attr(x="ΔybaK", y=effect[3]/2, text=rep[3], showarrow=false, font=attr(color="white")),
    # attr(x="Δgor", y=effect[4]/2, text=rep[4], showarrow=false, font=attr(color="white")),
    # attr(x="ΔyobF", y=effect[5]/2, text=rep[5], showarrow=false, font=attr(color="white")),
    # attr(x="ΔsrmB", y=effect[6]/2, text=rep[6], showarrow=false, font=attr(color="white"))],
    paper_bgcolor="white",
    plot_bgcolor="white",
    yaxis=attr(showline=true, linewidth=1, linecolor="black"),
    xaxis=attr(showline=true, linewidth=1, linecolor="black"),
    
    ))



savefig(fig, "different_conditions_rtc_expression_wo_replicates.jpeg")