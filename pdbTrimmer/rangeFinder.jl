function findRanges(arr::Vector{Int64})
    tmp = [arr[1]]
    out = []
    for x in 2:length(arr)
        println("$(arr[x]) - $(arr[x-1]) = $(arr[x] - arr[x-1])")
        if arr[x] - arr[x-1] == 1
            push!(tmp, arr[x-1])
        else
            push!(out, tmp)
            tmp = [arr[x]]
        end
    end
    return out
end

arr = [1,2,3,4,5,7,8,9,10]

findRanges(arr)