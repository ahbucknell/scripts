using Pandas
using CSV
using DataFrames

paths = ["/Users/doz23per/Documents/testPkl1.pkl","/Users/doz23per/Documents/testPkl2.pkl",]
outDir = ""


function extractPkl(path::String)
    pklDict = Pandas.read_pickle(path)
    outDict = Dict{String, Any}()
    outDict["ptm"] = pklDict["ptm"][1]
    outDict["avg_plddt"] = mean(pklDict["plddt"])
    return outDict, pklDict["predicted_aligned_error"], pklDict["plddt"]
end

csvData = DataFrames.DataFrame(binder =String[], pTM =Float64[], avg_pLDDT = Float64[])
for path in paths
    tmpName = String(split(basename(path), ".")[1])
    tmpOut, tmpPae, tmpPlddt = extractPkl(path)
    push!(csvData, [tmpName, tmpOut["ptm"], tmpOut["avg_plddt"]], promote=true)
    CSV.write(tmpName*"_PAE.csv", DataFrames.DataFrame(tmpPae, :auto), header = false)
    CSV.write(tmpName*"_pLDDT.csv", DataFrames.DataFrame(x=tmpPlddt), header = false)
end
CSV.write("pklData.csv", csvData)

