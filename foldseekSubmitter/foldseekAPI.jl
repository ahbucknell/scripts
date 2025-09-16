using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using HTTP, JSON3


# Given a PDB/mmCIF file, returns FoldSeek ticket for the job
function submitJob(file)
    path = joinpath(@__DIR__, file)
    # Below is the HTTP form sent to the FoldSeek API
    form = HTTP.Form(
    Dict(
        "q" => open(path),
        "mode" => "3diaa",
        "database[]" => "pdb100"))
    # Below sends the job request via HTTP POST request
    submission = HTTP.post("https://search.foldseek.com/api/ticket", [], form)
    # Job requests creates a ticket we'll need - this is embedded in the HTML response in JSON format as hex
    bodyData = JSON3.read(String(submission.body))
    return bodyData[:id]
end

# Given a FoldSeek ticket, waits for the job to complete and DLs the output.
function downloadResults(ticket, filename)
    target = "https://search.foldseek.com/api/ticket/" * ticket
    complete = false
    # Below will keep checking the ticket via HTTP GET requests until the job is done
    while complete == false
        sleep(1)
        response = HTTP.get(target)
        responseData = JSON3.read(String(response.body))
        if responseData[:status] == "COMPLETE"
            complete = true
        end
    end
    dlTarget = "http://search.foldseek.com/api/result/download/" * ticket
    destination = joinpath(@__DIR__, filename*".tar.gz")
    HTTP.download(dlTarget, destination)
end




file = "541_MGG_13622.pdb"
ticket = submitJob(file)
outname = split(file, ".")[1]
downloadResults(ticket, outname)



