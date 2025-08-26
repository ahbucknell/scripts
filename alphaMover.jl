using JSON, ArgParse

function parseArgs()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "arg1"
            help = "position of AlphaFold output directory"
            required = true
            arg_type = String
    end
    return parse_args(settings)
end

function findBestPkl(path)
    pklContents = JSON.parsefile(joinpath(path, "ranking_debug.json"))
    topPikl = pklContents["order"][1]
    return "result_" * topPikl * ".pkl"
end

function mvTargetFiles(subdir)
end

function main()
    args = parseArgs()
    targetDir = joinpath(@__DIR__, args["arg1"])
    subDirs = readdir(targetDir)
    notHidden = filter(f -> !startswith(f, "."), subDirs)
    for subDir in joinpath.(targetDir, notHidden)
        fileName = split(subDir, "/")[end]
        try
            targetPikl = findBestPkl(subDir)
            for file in [["ranked_0.pdb", ".pdb"], [targetPikl, ".pkl"]]
                old = joinpath(subDir, file[1])
                new = joinpath(dirname(subDir), fileName * file[2])
                mv(old, new)
            end
        catch
            println("Missing at least one necessary file. Skipping: ", fileName)
        end
    end
end

main()





