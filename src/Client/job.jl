function DFControl.DFJob(dir::String, s="localhost"; version::Int = -1)
    server = maybe_start_server(s)
    # server = Server(s)
    # dir = dir[1] == '/' ? dir[2:end] : dir
    resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    # Supplied dir was not a valid path, so we ask
    # previously registered jobs on the server that
    # contain dir.
    if resp.status == 204
        dir = request_job_dir(dir, server)
        dir === nothing && return
        resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    end 
    # @show String(resp.body)
    job = JSON3.read(resp.body, DFJob)
    job.server = server.name
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    return job
end

function request_job_dir(dir::String, server::Server)
    resp = HTTP.get(server, "/registered_jobs/" *  dir)
    matching_jobs = reverse(JSON3.read(resp.body, Vector{Tuple{String, DateTime}}))
    if length(matching_jobs) == 1
        return matching_jobs[1][1]
    elseif length(matching_jobs) == 0
        error("No jobs found matching $dir")
    elseif isdefined(Base, :active_repl)
        choices = ["$j -- $t" for (j, t) in matching_jobs]
        menu = RadioMenu(choices)
        choice = request("Multiple matching jobs were found, choose one:", menu)
        if choice != -1
            return matching_jobs[choice][1]
        else
            return nothing
        end
    else
        err_msg = "No concrete job for $dir found, closest matches are:"
        for m in choices
            err_msg = join([err_msg, m], "\n")
        end
        error(err_msg)
    end
end

function save(job::DFJob)
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    @assert !DFC.isarchived(job)
        "Not allowed to save a job in a archived directory, please specify a different directory with `set_dir!"
    server = maybe_start_server(job)
    @assert !isrunning(job) "Can't save a job in a directory where another is running."

    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    
    curver = job.version
    resp_job = JSON3.read(HTTP.post(server, "/jobs/" * job.dir, [], JSON3.write(job)).body, DFJob)
    @info "Job version: $(curver) => $(resp_job.version)."
    for f in fieldnames(DFJob)
        setfield!(job, f, getfield(resp_job, f))
    end
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    return job
end

"""
    isrunning(job::DFJob)

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_isrunning/" * job.dir).body, Bool)
end

"""
    versions(job::DFJob)

Returs the valid versions of `job`.
"""
function versions(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})
end

"""
    last_version(job::DFJob)

Returns the last version number of `job`.
"""
function last_version(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})[end]
end

"""
    last_running_calculation(job::DFJob)

Returns the last `DFCalculation` for which an output file was created.
"""
function last_running_calculation(job::DFJob)
    server = maybe_start_server(job)
    resp = HTTP.get(server, "/last_running_calculation", [], JSON3.write(job))
    if resp.status == 204
        return nothing
    else
        return job[JSON3.read(resp.body, Int)]
    end
end

function outputdata(job::DFJob)
    server = maybe_start_server(job)
    tmp_path = JSON3.read(HTTP.get(server, "/outputdata", [], JSON3.write(job)).body, String)
    local_temp = tempname() *".jld2"
    pull(server, tmp_path, local_temp)
    return DFC.JLD2.load(local_temp)["outputdata"]
end

"""
    pull(server::Server, server_file::String, local_file::String)

Pulls `server_file` from the server the `local_file`.
"""
function pull(server::Server, server_file::String, filename::String)
    if server.name == "localhost"
        cp(server_file, filename, force=true)
    else
        run(`scp $(ssh_string(server) * ":" * server_file) $filename`)
    end
end
