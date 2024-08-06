using BioStructures

function isSingleModel(struc)
    if length(models(inStructure)) > 1
        throw(DomainError(struc, "More than one model within the PDB!"))
    end
    return collect(keys(models(struc)))[1]
end

function isSingleChain(singleModel)
    if length(chains(singleModel)) > 1
        throw(DomainError(singleModel, "More than one chain within the PDB model!"))
    end
    return collect(keys(chains(struc[modelID])))[1]
end

function getPLDDTs(struc, modelNum, chainNum)
    plddtDict = Dict{Int, Float64}()
    for resNum in 1:length(residues(struc[modelNum][chainNum]))
        residueBfactor = tempfactor(struc[modelNum][chainNum][resNum]["N"])
        plddtDict[resNum] = residueBfactor
    end
    return plddtDict
end

function trimPLDDTs(resDict, minPLDDT)
    return filter(p -> (last(p) >= minPLDDT), resDict)
end

function findRanges(arr::Vector{Int64})
    tmp = [arr[1]]
    out = []
    for x in 2:length(arr)
        if arr[x] - arr[x-1] != 1
            push!(tmp, arr[x-1])
        else
            push!(out, tmp)
        end
    end
    return out
end

inStructure = read("MGG_00511T0.pdb", PDBFormat)
modelID = isSingleModel(inStructure)
chainID = isSingleChain(inStructure[modelID])
chainPLDDT = getPLDDTs(inStructure, modelID, chainID)
trimmedDict = trimPLDDTs(chainPLDDT, 20)
trimmedRes = sort(collect(keys(trimmedDict)))

list = [1,2,3,4,8,9,10,11]
findRanges(list)

