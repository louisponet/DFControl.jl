function read_job_info(job::Job)
    job_dir = job.dir
    orig_dir = job_dir
    scriptpath = joinpath(job_dir, "job.tt")
    if ispath(scriptpath)
        version = Jobs.version(job_dir)
    else
        error("No valid job found in $job_dir.")
    end
    metadata::Dict{Symbol, Any} = ispath(joinpath(job_dir, ".metadata.jld2")) ? JLD2.load(joinpath(job_dir, ".metadata.jld2"))["metadata"] : Dict{Symbol, Any}()

    name, tcalcs, header, environment = FileIO.read_job_script(scriptpath)
    cs = map(tcalcs) do x
        infile = splitdir(x.infile)[2]
        outfile = splitdir(x.outfile)[2]
        return (name = splitext(infile)[1], infile = infile, outfile = outfile, contents = read(x.infile, String), exec=x.exec, run=x.run)
    end
    
    # TODO Define a better way of working with pseudos
    upf_files = filter(x->occursin(".UPF",x), readdir(job_dir))
    pseudos = Dict([n => read(joinpath(job_dir, n), String) for n in upf_files])
    
    return Dict(:name => name,
                :version => version,
                :header => header,
                :environment => environment,
                :calculations => cs,
                :pseudos => pseudos) 
end

function load(job::Job)
    inf = read_job_info(job)
    name        = inf[:name]
    header      = inf[:header]
    environment = inf[:environment]
    calculations, structure = FileIO.parse_calculations(inf[:calculations])
    for a in structure.atoms
        a.pseudo = get(inf[:pseudos], a.pseudo, "")
    end
    return Job(name, structure, calculations, job.dir, header, inf[:version], job.copy_temp_folders, job.server, environment)
    
end

function save(jobdir::String, files; kwargs...)

    if jobdir[end] == '/'
        jobdir = jobdir[1:end-1]
    end
    #Since at this stage we know the job will belong to the current localhost we change the server
    # Here we find the main directory, needed for if a job's local dir is a .versions one
    dir = Jobs.main_job_dir(jobdir)
    version = Jobs.last_job_version(dir) + 1
    if ispath(joinpath(dir, "job.tt"))
        tj = load(Job(dir))
        cp(tj, joinpath(tj, Jobs.VERSION_DIR_NAME, "$(tj.version)"); force = true)
    end
    if jobdir != dir
        # We know for sure it was a previously saved job
        # Now that we have safely stored it we can clean out the directory to then fill
        # it with the files from the job.version
        clean_dir!(dir)
        for f in readdir(jobdir)
            cp(f, dir; force = true)
        end
    end

    
    # Needs to be done so the inputs `dir` also changes.
    mkpath(dir)

    for (name, f) in files
        d = splitdir(name)[1]
        mkpath(joinpath(dir, d))
        write(joinpath(dir, name), f)
    end

    JOB_QUEUE[].full_queue[dir] = (-1, Jobs.Saved)
    return version
end

job_versions(args...) = Jobs.job_versions(args...)

function registered_jobs(jobdir::String)
    queue!(JOB_QUEUE[], Servers.local_server().scheduler, false)
    dirs = Tuple{String,DateTime}[]
    queue = JOB_QUEUE[]
    for qu in (queue.current_queue, queue.full_queue)
        for dir in keys(qu)
            if occursin(jobdir, dir)
                push!(dirs, (dir, Jobs.timestamp(dir)))
            end
        end
    end
    return sort(dirs, by = x -> x[2], rev=true)
end

function running_jobs(fuzzy)
    out = Tuple{String, Int}[]
    for (j, info) in JOB_QUEUE[].current_queue
        if occursin(fuzzy, j)
            push!(out, (j, info[1]))
        end
    end
    return out
end

function workflow_logger(job::Job)
    return TeeLogger(MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                               append = true), Logging.Info),
                     MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                               append = true), Logging.Warn),
                     MinLevelLogger(FileLogger(joinpath(job, ".workflow", "error.log");
                                               append = true), Logging.Error))
end

queued_dir(job::Job, args...) = joinpath(job, ".workflow/queued", args...)
finished_dir(job::Job, args...) = joinpath(job, ".workflow/finished", args...)

function clear_queue!(job::Job)
    qd = queued_dir(job)
    if !ispath(qd)
        return
    end
    for f in readdir(qd)
        rm(joinpath(qd, f))
    end
end

queued(job::Job) = readdir(queued_dir(job))
finished(job::Job) = readdir(finished_dir(job))

# function submit_workflow(job::Job, funcs, d::Daemon = init_daemon())
#     queue_steps(job, funcs)
#     write_workflow_files(job)
#     while !is_started(d)
#         sleep(1)
#     end
#     return runexpr("""
#                DFControl.spawn_worker(DAEMON, Job("$(job.dir)"))
#                DFControl.save(DAEMON)
#                """; port = d.port)
# end

"""
    submit(dir::String, workflow::Bool)

Writes the directory to either pending workflows or pending jobs file.
"""
function submit(job_dir::String, workflow::Bool)
    if workflow
        open(PENDING_WORKFLOWS_FILE, "a", lock=true) do f
            return write(f, job_dir * "\n")
        end
    else
        open(PENDING_JOBS_FILE, "a", lock=true) do f
            return write(f, job_dir * "\n")
        end
    end
end

"""
    last_running_calculation(path::String)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(path::String)
    scrpath = joinpath(path, "job.tt")
    job = load(Job(path))
    t = mtime(Jobs.scriptpath(job))
    times = map(x -> (o = joinpath(job, x.outfile); ispath(o) ? mtime(o) : 0.0), job.calculations)
    return isempty(times) ? nothing : findmax(times)[2]
end

"""
    state(job_dir::String)

Returns the job state of the job in `job_dir`.
"""
function state(job_dir::String)
    if haskey(JOB_QUEUE[].current_queue, job_dir)
        return JOB_QUEUE[].current_queue[job_dir][2]
    elseif haskey(JOB_QUEUE[].full_queue, job_dir)
        return JOB_QUEUE[].full_queue[job_dir][2]
    else
        return Jobs.Unknown
    end
end

function dirsize(path::String)
    totsize = 0.0
    for (root, dirs, files) in walkdir(path)
        for file in files
            totsize += filesize(root, file)
        end
    end
    return totsize
end

"""
    filesize(job::Job)

Total filesize on disk for a job and all its versions.
"""
Base.filesize(job::Job) = dirsize(job.dir)

timestamp(job::Job) = timestamp(job.dir)
has_timestamp(job) = haskey(job.metadata, :timestamp)

function clean_dir!(dir::AbstractString)
    for f in readdir(dir)
        if f == Jobs.TEMP_CALC_DIR ||
           f == Jobs.VERSION_DIR_NAME ||
           splitext(f)[end] == ".jl"
            continue
        end
        rm(joinpath(dir, f); recursive = true)
    end
end

exists_job(d::AbstractString) = ispath(d) && ispath(joinpath(d, "job.tt"))

function outputdata(jobdir::String, calculations::Vector{String})
    job = load(Job(jobdir))
    calculations = isempty(calculations) ? map(x->x.name, job.calculations) : calculations
    respath = joinpath(job, "results.jld2")
    if ispath(respath)
        datadict = JLD2.load(respath, "outputdata")
    else
        datadict = Dict{String,Dict{Symbol,Any}}()
    end
    stime = isempty(datadict) ? 0.0 : mtime(respath)
    new_data = false
    for c in calculations
        calculation = job[c]
        p = joinpath(job, calculation.outfile)
        if mtime(p) > stime
            try 
                tout = outputdata(calculation, p)
                if !isempty(tout)
                    datadict[calculation.name] = tout
                    new_data = true
                end
            catch e
                @warn "Something went wrong reading output for calculation $c."
                @warn e
            end
        end
    end
    if new_data
        JLD2.jldsave(respath; outputdata=datadict)
        return respath
    elseif ispath(respath)
        return respath
    else
        return nothing
    end
end

rm_version!(jobdir::String, version::Int) = Jobs.rm_version!(load(Job(jobdir)), version)

function Servers.abort(job_dir::String)
    @assert haskey(JOB_QUEUE[].current_queue, job_dir) "No job exists in dir: $job_dir!"
    s = Servers.local_server()
    id = JOB_QUEUE[].current_queue[job_dir][1]
    Servers.abort(s.scheduler, id)
    delete!(JOB_QUEUE[].current_queue, job_dir)
    JOB_QUEUE[].full_queue[job_dir] = (id, Jobs.Cancelled)
    return id
end
