using CSV
using CairoMakie


ptmMx = CSV.read(joinpath(@__DIR__, "PTM.csv"),CSV.Tables.matrix, header=1, types=Float64)
iptmMx = CSV.read(joinpath(@__DIR__, "iPTM.csv"),CSV.Tables.matrix, header=1, types=Float64)

score = Float64[]
for x in eachindex(ptmMx)
    push!(score, 0.2*ptmMx[x] + 0.8*iptmMx[x])
end
scoreMx = reshape(score, (5,8))

ticklistSep = (1:8, ["Pil1","Pil2","Slm1/2","Eis1","GFA-1","Sep4","Sep5","Sep6"])
ticklistx = (1:5, ["Pil1","Pil2","Slm1/2","Eis1","GFA-1"])

doubleFig = Figure()
ax1 = Axis(doubleFig[1, 1], yreversed = true, title = "pTM", xticks = ticklistx, yticks = ticklistSep, aspect = 1, xticklabelsize = 12, yticklabelsize = 12)
ptmhm = heatmap!(ax1, ptmMx, colorrange = (0,0.75))
ax2 = Axis(doubleFig[1,2], yreversed = true, title = "ipTM", xticks = ticklistx, yticks = ticklistSep,  aspect = 1, xticklabelsize = 12, yticklabelsize = 12)
iptmhm = heatmap!(ax2, iptmMx, colorrange = (0,0.75))
Colorbar(doubleFig[1,3], ptmhm)
rowsize!(doubleFig.layout, 1, Aspect(2, 1))

singleFig = Figure()
ax3 = Axis(singleFig[1,1], yreversed = true, title = "Combined pTM/ipTM score", xticks = ticklistx, yticks = ticklistSep, aspect = 1)
combinedhm = heatmap!(ax3, scoreMx, colorrange = (0,0.75))
Colorbar(singleFig[1,2], combinedhm)
singleFig

save("double_heatmap.svg", doubleFig)