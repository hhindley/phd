function plotBIG(x, y; xtitle="", ytitle="", title="")
    xvals = range(minimum(x), maximum(x), length=length(x))

    f = Figure()
    ax = Axis(f[1, 1], xlabel = "$xtitle", ylabel = "$ytitle", title="$title")
    ilines!(f[1,1], xvals, y)

    return f
end


function plot_subplots(data, species, threshold_vals)

    f = Figure()
    axess=[]

    for i in 1:5
        ax = Axis(f[i, 1], xlabel = "time", ylabel = "rm_a", title="threshold: $(threshold_vals[i])")
        push!(axess, ax)
    end

    [hidexdecorations!(ax) for ax in axess[1:end-1]]
    [ilines!(f[i,1], makiex(data[i].time), data[i][:,species]) for i in 1:5]

    return f
end

function makiex(x)
    return range(minimum(x), maximum(x), length=length(x))
end