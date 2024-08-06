function count_nucleotides(strand)
    out = Dict{Char, Int}('A' => 0, 'T' => 0, 'G' => 0, 'C' => 0)
    for i in unique(strand)
        if i in ['A','T','G','C']
            out[i] += 1
        else
            throw(DomainError(i, "Argument contains invalid nucleotide"))
        end
    end
    print(out)
end

count_nucleotides()




#=
    valid = true
    for i in unique(strand)
        if !(i in ['A','T','G','C']) || length(strand) < 1
            print("This is invalid!")
        else
            push!(out, i => count(==(i),strand))
        end
    end
    if valid == false
        print("INVALID -> error")
    else
        print(string, out)
    end
end
count_nucleotides("AGC")=#