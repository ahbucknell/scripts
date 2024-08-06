using Pandas
using CairoMakie

function makePAE(path::String, interval::Float64)
    pklDict = Pandas.read_pickle(path)
    arrlen = sqrt(length(pklDict["predicted_aligned_error"]))
    f = Figure()
    ax1 = Axis(f[1,1], yreversed = true, title = "Predicted Aligned Error (PAE)", aspect = 1, xticks = 0:interval:arrlen, yticks = 0:interval:arrlen)
    hm = heatmap!(ax1, pklDict["predicted_aligned_error"], colorrange = (0,30), colormap = Reverse(:BuGn))
    return f
end


