# This script will convert all pickle files in the current directoy to JSON
# It assumes that all pickle files within the directory are generated from AlphaFold2
# Most importantly, it assumes that YOU trust all the pickle files within this directory!
ENV["JULIA_CONDAPKG_VERBOSITY"] = "-1"
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using PythonCall, JSON3, CondaPkg, ArgParse

# Manage dependencies for Python environment
pkl = pyimport("pickle")
CondaPkg.add("numpy")
CondaPkg.resolve()

function parseArgs()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "arg1"
            help = "path of directory to translate"
            required = true
            arg_type = String
    end
    return parse_args(settings)
end

function readPkl(fpath)
    file = pybuiltins.open(fpath, "rb")
    try
        data = pkl.load(file)
        pyconvert(Dict, data)
    finally
        file.close()
    end
end

# Converts Python objects to Julia counterparts
# Scalars to numbers, arrays to array, and matrices to nested arrays
function translateTypes(elem)
    if ndims(elem) == 0
        return elem[] 
    elseif ndims(elem) == 2
        [vec(elem[i, :]) for i in 1:size(elem, 1)]
    else
        return pyconvert(Array, elem)
    end
end

function convertFile(fPath)
    targetKeys = ["ptm", "plddt", "predicted_aligned_error", "max_predicted_aligned_error"]
    pklData = readPkl(fPath)
    outData = Dict(k => translateTypes(pklData[k]) for k in targetKeys)
    filename = split(basename(fPath), ".")[1]*".json"
    oPath = joinpath(dirname(fPath), filename)
    open("$oPath", "w") do io
        JSON3.write(io, outData)
    end
end


function main()
    args = parseArgs()
    targetDir = abspath(args["arg1"])
    files = readdir(targetDir)
    pkls = filter(f -> endswith(f, ".pkl"), files)
    for p in pkls
        path = joinpath(targetDir, p)
        convertFile(path)
    end
end

main()





