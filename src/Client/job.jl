Servers.Server(j::Job) = Server(j.server)

function Jobs.Job(dir::AbstractString, s = "localhost"; version::Int = -1)
    server = Servers.maybe_start_server(s)
    if !ispath(server, dir)
        dir = request_job_dir(dir, server)
        dir === nothing && return
    end
    if version == last_version(dir; server= server)
        dir = Jobs.main_job_dir(dir)
    elseif !occursin(Jobs.VERSION_DIR_NAME, dir) && version != -1
        dir = Jobs.version_dir(Jobs.main_job_dir(dir), version)
    end
    resp = HTTP.get(server, "/jobs/" * abspath(server, dir))
    # Supplied dir was not a valid path, so we ask
    # previously registered jobs on the server that
    # contain dir.
    job = JSON3.read(resp.body, Job)
    job.server = server.name
    # if haskey(job.metadata, :timestamp)
    #     job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    # end
    return job
end

function request_job_dir(dir::String, server::Server)
    resp = HTTP.get(server, "/registered_jobs/" * dir)
    matching_jobs = reverse(JSON3.read(resp.body, Vector{Tuple{String,DateTime}}))
    choices = ["$j -- $t" for (j, t) in matching_jobs]
    if length(matching_jobs) == 1
        return matching_jobs[1][1]
    elseif length(matching_jobs) == 0
        error("No jobs found matching $dir")
    elseif isdefined(Base, :active_repl)
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

function save(job::Job, server::Server = Servers.maybe_start_server(job))
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    @assert !isrunning(job, server) "Can't save a job in a directory where another is running."

    @assert !Jobs.isarchived(job)
    "Not allowed to save a job in a archived directory, please specify a different directory."

    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    Structures.sanitize!(job.structure)
    Calculations.sanitize_flags!(job.calculations, job.structure, job.name,
                                 "./"*Jobs.TEMP_CALC_DIR)

    Jobs.sanitize_cutoffs!(job)

    curver = job.version

    tmpdir = mkpath(tempname())
    apath = abspath(job)

    if job.environment != ""
        environment = get_environment(job.environment; server= job.server)
        @assert environment !== nothing "Environment with name $(job.environment) not found!"
    else
        environment = nothing
    end
    job.dir = tmpdir
    FileIO.write_job_files(job, environment)
    files_to_send = Dict([f => read(joinpath(tmpdir, f), String) for f in readdir(tmpdir)])
    
    job.dir = apath
    rm(tmpdir, recursive=true)
  
    resp_job_version = JSON3.read(HTTP.post(server, "/jobs/" * apath, [], JSON3.write(files_to_send)).body,
                          Int)
    @info "Job version: $(curver) => $(resp_job_version)."
    job.version = resp_job_version
    Calculations.rm_tmp_flags!.(job.calculations)
    return job
end

function submit(job::Job)
    server = Servers.maybe_start_server(job)
    verify_execs(job, server)
    save(job)
    return HTTP.put(server, "/jobs/" * abspath(job))
end

function submit(jobs::Vector{Job}, run = true)
    # To verify all execs
    server_names = unique(map(x -> x.server, jobs))
    buckets = [jobs[findall(x -> x.server == s, jobs)] for s in server_names]
    outbuckets = [Vector{Job}(undef, length(b)) for b in buckets]
    Threads.@threads for i in 1:length(server_names)
        server = Servers.maybe_start_server(server_names[i])
        bucket = buckets[i]
        outbucket = outbuckets[i]
        execs = unique(vcat([map(x->x.exec, j.calculations) for j in bucket]...))
        
        replacements = verify_execs(execs, server)
        Threads.@threads for ij in 1:length(bucket)
            job = bucket[ij]
            for c in job.calculations
                for (e, rep) in replacements
                    if c.exec == e
                        c.exec.dir = rep.dir
                        c.exec.modules = rep.modules
                    end
                end
            end
            outbucket[ij] = save(job, server)
            run && HTTP.put(server, "/jobs/" * abspath(job))
        end
    end
    return outbuckets
end

"""
    state(job::Job, [server::Server])

Returns the state of a job.
"""
function state(job::Job, server::Server = Servers.maybe_start_server(job))
    return JSON3.read(HTTP.get(server, "/job_state/" * abspath(job)).body, Jobs.JobState)
end
    
"""
    isrunning(job::Job, [server::Server])

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job::Job, server::Server = Servers.maybe_start_server(job))
    s = state(job, server)
    return s == Jobs.Running || s == Jobs.Pending
end

"""
    versions(job::Job)

Returs the valid versions of `job`.
"""
versions(job::Job) = versions(abspath(job); server = job.server) 
function versions(dir::AbstractString; server = "localhost")
    server=Servers.maybe_start_server(server) 
    return JSON3.read(HTTP.get(server, "/job_versions/" * Jobs.main_job_dir(dir)).body, Vector{Int})
end

"""
    last_version(job::Job)

Returns the last version number of `job`.
"""
last_version(job::Job) = last_version(abspath(job), server = job.server) 
function last_version(dir::AbstractString; server="localhost")
    t = versions(dir, server=server)
    return isempty(t) ? 0 : t[end]
end

function submission_time(job::Job)
    server = Servers.maybe_start_server(job)
    resp = HTTP.get(server, "/job_submission_time/" * abspath(job))
    return JSON3.read(resp.body, Float64)
end

"""
    last_running_calculation(job::Job)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(job::Job)
    server = Servers.maybe_start_server(job)
    resp = HTTP.get(server, "/last_running_calculation/" * abspath(job))
    if resp.status == 204
        return nothing
    else
        return job[JSON3.read(resp.body, Int)]
    end
end

outputdata(job::Job; kwargs...) =
    outputdata(job.dir; server=job.server, kwargs...)
    
function outputdata(jobdir::String; server = "localhost", calcs::Vector{String}=String[], extra_parse_funcs = nothing)
    server = Servers.maybe_start_server(server)
    jobdir = isabspath(jobdir) ? jobdir : joinpath(server, jobdir)
    resp = HTTP.get(server, "/outputdata/" * jobdir, [], JSON3.write(calcs))
    if resp.status == 204
        error("No outputdata found yet. Is the job running?")
    end
    tmp_path = JSON3.read(resp.body,
                          String)
                          
    local_temp = tempname() * ".jld2"
    Servers.pull(server, tmp_path, local_temp)
    dat = JLD2.load(local_temp)
    rm(local_temp)
    out = dat["outputdata"]
    if extra_parse_funcs !== nothing
        for k in keys(out)
            if !isempty(calcs) && k ∈ calcs
                n = c.name
                try
                    f = joinpath(jobdir, c.outfile)
                    local_f = tempname()
                    Servers.pull(server, f, local_f)
                    FileIO.parse_file(local_f, extra_parse_funcs, out = out[n])
                    rm(local_f)
                catch
                    nothing
                end
            end
        end
    end
    return out     
end

function known_execs(e::String, dir::String = ""; server = Server("localhost"))
    s = Servers.maybe_start_server(server)
    return JSON3.read(HTTP.get(s, "/known_execs/", [], JSON3.write(Dict("exec" => e, "dir" => dir))).body, Vector{Calculations.Exec})
end
known_execs(e::Calculations.Exec; kwargs...) = known_execs(e.exec,e.dir; kwargs...)

function get_exec(name::String; server="localhost")
    s = Servers.maybe_start_server(server)
    return JSON3.read(HTTP.get(s, "/exec/$name").body, Calculations.Exec)
end
function save(e::Exec; server="localhost")
    s = Servers.maybe_start_server(server)
    return JSON3.read(HTTP.post(s, "/exec/", [], JSON3.write(e)).body, Calculations.Exec)
end

function verify_execs(job::Job, server::Server)
    replacements = verify_execs(unique(map(x->x.exec, filter(x->x.run, job.calculations))), server)
    for (e, rep) in replacements
        for c in job.calculations
            if c.exec == e
                c.exec.dir = rep.dir
                c.exec.modules = rep.modules
            end
        end
    end
end

function verify_execs(execs::Vector{Exec}, server::Server)
    replacements = Dict{Exec, Exec}()
    for e in execs 
        if !JSON3.read(HTTP.get(server, "/verify_exec/", [], JSON3.write(e)).body, Bool)
            possibilities = known_execs(e, server=server)
            curn = 0
            while length(possibilities) > 1
                possibilities = filter(x -> all(splitpath(x.dir)[end-curn:end] .== splitpath(e.dir)[end-curn:end]), possibilities)
                curn += 1
            end
            if !isempty(possibilities)
                replacement = possibilities[1]
                @warn "Executable ($(e.exec)) in dir ($(e.dir)) not runnable,\n but found a matching replacement executable in dir ($(replacement.dir)).\nUsing that one..."
                replacements[e] = replacement 
            else
                error("$e is not a valid executable on server $(server.name).\nReplace it.")
            end
        end
    end
    return replacements
end

#TODO: work this
"""
    abort(job::Job)

Will try to remove the job from the scheduler's queue.
If the last running calculation happened to be a `Calculation{QE}`, the correct abort file will be written.
For other codes the process is not smooth, and restarting is not guaranteed.
"""
function abort(job::Job)
    @assert isrunning(job) "Is this job running?"
    id = JSON3.read(HTTP.get(Server(job.server), "/abort/" * abspath(job)).body, Int)
    @info "Aborted job $id"
end

function registered_jobs(fuzzy::String=""; server=nothing)
    if server === nothing
        all_jobs = Dict{String,Vector{Tuple{String, DateTime}}}()
        for s in Servers.known_servers()
            jobs = registered_jobs(fuzzy, s)
            if !isempty(jobs)
                all_jobs[s.name] = jobs
            end
        end
        return all_jobs
    else
        server = Servers.maybe_start_server(server) 
        resp = HTTP.get(server, "/registered_jobs/" * fuzzy)
        return reverse(JSON3.read(resp.body, Vector{Tuple{String,DateTime}}))
    end
end
      
function running_jobs(fuzzy=""; server="localhost")
    server = Servers.maybe_start_server(server) 
    resp = HTTP.get(server, "/running_jobs/" * fuzzy)
    return reverse(JSON3.read(resp.body, Vector{String}))
end

function switch_version!(job::Job, version::Int)
    allvers = versions(job)
    if !(version in allvers)
        error("Version $version does not exist.")
    end
    tj = Job(Jobs.main_job_dir(job), job.server, version=version)
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

    HTTP.put(server, "/rm_version/" * abspath(job), [], JSON3.write(version))    
    if version == job.version
        @warn "Job version is the same as the one to be removed, switching to last known version."
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        switch_version!(job, lv)
    end
end

function environment_from_jobscript(scriptpath; server="localhost")
    server = Servers.maybe_start_server(server)
    tmp = tempname()
    Servers.pull(server, scriptpath, tmp)
    return Jobs.environment_from_jobscript(tmp)
end
    
function get_environment(name; server="localhost")
    server = Servers.maybe_start_server(server)
    return JSON3.read(HTTP.get(server, "/environment/$name").body, Environment)    
end
    
function add_environment(env::Environment, name::String; server="localhost")
    server = Servers.maybe_start_server(server)
    return HTTP.post(server, "/environment/$name", [], JSON3.write(env))
end

function rm_environment!(name::String; server="localhost")
    server = Servers.maybe_start_server(server)
    return HTTP.put(server, "/environment/$name")
end
 
