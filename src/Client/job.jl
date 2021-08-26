Servers.Server(j::Job) = Server(j.server)

function Jobs.Job(dir::AbstractString, s = "localhost"; version::Int = -1)
    server = Servers.maybe_start_server(s)
    if occursin(Jobs.VERSION_DIR_NAME, dir)
        error("It is not allowed to directly load a job version, please use `Job($dir, version=$(splitdir(dir)[end]))`")
    end
    resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    # Supplied dir was not a valid path, so we ask
    # previously registered jobs on the server that
    # contain dir.
    if resp.status == 204
        dir = request_job_dir(dir, server)
        dir === nothing && return
        resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    end
    job = JSON3.read(resp.body, Job)
    job.server = server.name
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    return job
end

function request_job_dir(dir::String, server::Server)
    resp = HTTP.get(server, "/registered_jobs/" * dir)
    matching_jobs = reverse(JSON3.read(resp.body, Vector{Tuple{String,DateTime}}))
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

function save(job::Job)
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    server = Servers.maybe_start_server(job)

    @assert !isrunning(job) "Can't save a job in a directory where another is running."

    @assert !Jobs.isarchived(job)
    "Not allowed to save a job in a archived directory, please specify a different directory with `set_dir!"

    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    Structures.sanitize!(job.structure)
    Calculations.sanitize_flags!(job.calculations, job.structure, job.name,
                                 joinpath(job, Jobs.TEMP_CALC_DIR))


    Jobs.sanitize_cutoffs!(job)
    files_to_copy = sanitize_pseudos!(job)

    curver = job.version
    resp_job = JSON3.read(HTTP.post(server, "/jobs/" * job.dir, [], JSON3.write(job)).body,
                          Job)
    for f in files_to_copy
        Servers.push(f, server, joinpath(job.dir, splitdir(f)[end]))
    end

    @info "Job version: $(curver) => $(resp_job.version)."
    for f in fieldnames(Job)
        if f == :server
            continue
        end
        setfield!(job, f, getfield(resp_job, f))
    end
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    Calculations.rm_tmp_flags!.(job.calculations)
    return job
end

function submit(job::Job)
    server = Servers.maybe_start_server(job)
    verify_execs(job, server)
    save(job)
    return HTTP.put(server, "/jobs/" * job.dir)
end

function sanitize_pseudos!(job::Job)
    all_pseudos = map(x -> x.pseudo, job.structure.atoms)
    uni_dirs = unique(map(x -> x.dir, all_pseudos))
    uni_pseudos = unique(all_pseudos)
    s = Server(job)
    pseudo_paths = Structures.path.(uni_pseudos)
    if !all(x -> JSON3.read(HTTP.get(s, "/get_ispath/" * x).body, Bool), pseudo_paths)
        if all(ispath, pseudo_paths)
            @info "Some Pseudos could not be found on Server $(s.name), pushing them from local storage to job dir."
            for p in all_pseudos
                p.dir = job.dir
            end
            return pseudo_paths
        else
            # Find if the pseudos are part of a pseudoset
            if length(uni_dirs)==1
                for set in list_pseudosets(s)
                    pseudos = pseudos(s, set)
                    dir = values(pseudos)[1].dir
                    pseudo_names = map(x -> x.name, values(pseudos))
                    if splitpath(uni_dirs[1])[end] == splitpath(dir)[end] && all(x -> x.name ∈ pseudo_names, uni_pseudos)
                        @warn "Matching pseudoset found on Server: $set."
                        set_pseudos!(job, set)
                        return String[]
                    end
                end
            else
                @warn "Some pseudos could not be found locally, and neither on the remote server."
            end
        end
    end
    return String[]
end

"""
    isrunning(job::Job)

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job::Job)
    server = Servers.maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_isrunning/" * job.dir).body, Bool)
end

"""
    versions(job::Job)

Returs the valid versions of `job`.
"""
versions(job::Job) = versions(job.dir, job.server) 
function versions(dir, s = "localhost")
    server=Servers.maybe_start_server(s) 
    return JSON3.read(HTTP.get(server, "/job_versions/" * Jobs.main_job_dir(dir)).body, Vector{Int})
end

"""
    last_version(job::Job)

Returns the last version number of `job`.
"""
last_version(job::Job) = last_version(job.dir, job.server) 
function last_version(dir, s="localhost")
    server = Servers.maybe_start_server(s)
    return JSON3.read(HTTP.get(server, "/job_versions/" * Jobs.main_job_dir(dir)).body, Vector{Int})[end]
end

"""
    last_running_calculation(job::Job)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(job::Job)
    server = Servers.maybe_start_server(job)
    resp = HTTP.get(server, "/last_running_calculation", [], JSON3.write(job))
    if resp.status == 204
        return nothing
    else
        return job[JSON3.read(resp.body, Int)]
    end
end

function outputdata(job::Job; extra_parse_funcs = nothing)
    server = Servers.maybe_start_server(job)
    tmp_path = JSON3.read(HTTP.get(server, "/outputdata", [], JSON3.write(job)).body,
                          String)
    local_temp = tempname() * ".jld2"
    Servers.pull(server, tmp_path, local_temp)
    dat = JLD2.load(local_temp)
    rm(local_temp)
    out = dat["outputdata"]
    if extra_parse_funcs !== nothing
        for c in job.calculations
            n = c.name
            if haskey(out, n)
                # try
                    f = Calculations.outpath(c)
                    local_f = tempname()
                    Servers.pull(server, f, local_f)
                    FileIO.parse_file(local_f, extra_parse_funcs, out = out[n])
                    rm(local_f)
                # catch
                #     nothing
                # end
            end
        end
    end
    return out     
end

function verify_execs(job::Job, server::Server)
    for e in unique(map(x->x.exec, job.calculations))
        if !JSON3.read(HTTP.get(server, "/verify_exec/", [], JSON3.write(e)).body, Bool)
            possibilities = JSON3.read(HTTP.get(server, "/known_execs/" * e.exec).body, Vector{Calculations.Exec})
            replacement = length(possibilities) == 1 ? possibilities[1] : getfirst(x -> x.dir == e.dir, possibilities)
            if replacement !== nothing
                @warn "Modules mismatched, but found a matching replacement executable on the server with the correct modules.\nUsing that one..."
                for e1 in map(x->x.exec, job.calculations)
                    if e1.exec == replacement.exec
                        e1.modules = replacement.modules
                        e1.dir = replacement.dir
                    end
                end
            else
                error("$e is not a valid executable on server $(server.name).\nReplace it with one of the following ones:\n$possibilities")
            end
        end
    end
end

#TODO: work this
"""
    abort(job::Job)

Will try to remove the job from the scheduler's queue.
If the last running calculation happened to be a `Calculation{QE}`, the correct abort file will be written.
For other codes the process is not smooth, and restarting is not guaranteed.
"""
function abort(job::Job)
    lastrunning = job.calculations[last_running_calculation(job)]
    if lastrunning == nothing
        error("Is this job running?")
    end
    if eltype(lastrunning) == QE
        length(filter(x -> eltype(x) == QE, job.calculations)) > 1 &&
            @warn "It's absolutely impossible to guarantee a graceful abort of a multi job script with QE."

        abortpath = writeabortfile(job, lastrunning)
        while ispath(abortpath)
            continue
        end
        qdel(job)
    else
        qdel(job)
    end
end

function registered_jobs(fuzzy="", s="localhost")
    server = Servers.maybe_start_server(s) 
    resp = HTTP.get(server, "/registered_jobs/" * fuzzy)
    return reverse(JSON3.read(resp.body, Vector{Tuple{String,DateTime}}))
end

function running_jobs(fuzzy="", s="localhost")
    server = Servers.maybe_start_server(s) 
    resp = HTTP.get(server, "/running_jobs/" * fuzzy)
    return reverse(JSON3.read(resp.body, Vector{String}))
end

function switch_version!(job::Job, version::Int)
    allvers = versions(job)
    if !(version in allvers)
        error("Version $version does not exist.")
    end
    tj = Job(Jobs.main_job_dir(job.dir), job.server, version=version)
    for f in fieldnames(Job)
        setfield!(job, f, getfield(tj, f))
    end
    return job
end

function rm_version!(job::Job, version::Int)
    server = Servers.maybe_start_server(job.server) 
    allvers = versions(job)
    if !(version in allvers)
        error("Version $version does not exist.")
    end

    HTTP.put(server, "/rm_version/$version", [], JSON3.write(job))    
    if version == job.version
        @warn "Job version is the same as the one to be removed, switching to last known version."
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        switch_version!(job, lv)
    end
end

function environment_from_jobscript(scriptpath, s="localhost")
    server = Servers.maybe_start_server(s)
    tmp = tempname()
    Servers.pull(server, scriptpath, tmp)
    return Jobs.environment_from_jobscript(tmp)
end
    
function get_environment(name, s="localhost")
    server = Servers.maybe_start_server(s)
    return JSON3.read(HTTP.get(server, "/environment/$name").body, Environment)    
end
    
function add_environment(env::Environment, name::String, s="localhost")
    server = Servers.maybe_start_server(s)
    return HTTP.post(server, "/environment/$name", [], JSON3.write(env))
end

function rm_environment!(name::String, s="localhost")
    server = Servers.maybe_start_server(s)
    return HTTP.put(server, "/environment/$name")
end
    
    
