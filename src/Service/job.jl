function load_job(job_dir::AbstractString, version::Int=-1)
    if ispath(job_dir)
        if occursin(Jobs.VERSION_DIR_NAME, job_dir)
            error("It is not allowed to directly load a job version, please use `Job(dir, version=$(splitdir(job_dir)[end]))`")
        end
        if version != -1
            real_path = version_dir(job_dir, version)
            real_version = version
        elseif ispath(joinpath(job_dir, "job.tt"))
            real_path = job_dir
            real_version = Jobs.main_job_version(job_dir)
        else
            error("No valid job found in $job_dir.")
        end
        job = Job(; merge((dir = real_path, version = real_version),
                    Jobs.read_job_calculations(joinpath(real_path, "job.tt")))...)
        Jobs.maybe_register_job(job)
        return job
    else
        return nothing
    end
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

function queue_steps(job::Job, funcs)
    qd = queued_dir(job)
    if !ispath(qd)
        mkpath(qd)
    end
    prev_steps = readdir(qd)
    last_step = !isempty(prev_steps) ? parse(Int, splitext(prev_steps[end])[1]) : 0
    for (i, f) in enumerate(funcs)
        write(joinpath(qd, "$(i + last_step).jl"), @code_string f(job, Dict{Symbol,Any}()))
    end
end

function write_workflow_files(job::Job)
    all = string.(values(Base.loaded_modules))
    valid = String[]
    ks = keys(Pkg.project().dependencies)
    for p in all
        if p ∈ ks
            push!(valid, p)
        end
    end
    DFControl.JLD2.save(joinpath(job, ".workflow/environment.jld2"), "modules", valid, "project",
              Base.current_project())
    return DFControl.JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", Dict{Symbol,Any}())
end


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

mods_test() = Base.loaded_modules

"""
    save(job::Job)

Saves the job's calculations and `job.tt` submission script in `job.dir`.
Some sanity checks will be performed on the validity of flags, execs, pseudopotentials, etc.
The job will also be registered for easy retrieval at a later stage.

If a previous job is present in the job directory (indicated by a valid job script),
it will be copied to the `.versions` sub directory as the previous version of `job`,
and the version of `job` will be incremented. 
"""
function save(job::Job; kwargs...)

    #Since at this stage we know the job will belong to the current localhost we change the server
    job.server = "localhost"
    # Here we find the main directory, needed for if a job's local dir is a .versions one
    dir = Jobs.main_job_dir(job)
    if dir != job.dir
        # We know for sure it was a previously saved job
        # Now that we have safely stored it we can clean out the directory to then fill
        # it with the files from the job.version
        clean_dir!(dir)
        cp(job, dir; force = true)
    end
    if ispath(joinpath(dir, "job.tt"))
        tj = load_job(dir, -1)
        cp(tj, joinpath(tj, Jobs.VERSION_DIR_NAME, "$(tj.version)"); force = true)
    end
    
    set_dir!(job, dir) # Needs to be done so the inputs `dir` also changes.
    mkpath(dir)

    job.version = last_version(job) + 1
    timestamp!(job, now())
    save_metadata(job)
    writejobfiles(job; kwargs...)
    maybe_register_job(job)
    return job
end

"""
    submit(job::Job; kwargs...)

First saves the job, then tries to submit the job script through `sbatch job.tt` or `bash job.tt` if the former command fails.
`kwargs...` get passed to `save(job; kwargs...)`.
"""
function submit(job_dir::String)
    open(PENDING_JOBS_FILE, "a") do f
        write(f, job_dir * "\n")
    end
end


"""
    last_running_calculation(job::Job)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(job::Job)
    t = mtime(DFC.scriptpath(job))
    for (i, c) in enumerate(reverse(job.calculations))
        p = DFC.outpath(c)
        if ispath(p) && mtime(p) > t
            return i
        end
    end
end

for (f, strs) in zip((:cp, :mv), (("copy", "Copies"), ("move", "Moves")))
    @eval begin
        """
            $($f)(job::Job, dest::AbstractString; all=false, temp=false, kwargs...)

        $($(strs[2])) the contents of `job.dir` to `dest`. If `all=true`, it will also $($(strs[1])) the
        `.version` directory with all previous versions. If `temp=true` it will override
        `job.copy_temp_folders` and $($(strs[1])) also the temporary calculation directories.
        The `kwargs...` are passed to `Base.$($f)`.
        """
        function Base.$f(job::Job, dest::AbstractString; all = false, temp = false,
                         kwargs...)
            if !ispath(dest)
                mkpath(dest)
            end
            for file in readdir(job.dir)
                if !all
                    if file == VERSION_DIR_NAME
                        continue
                    elseif file == Jobs.TEMP_CALC_DIR && !(temp || job.copy_temp_folders)
                        continue
                    end
                end
                if joinpath(job, file) == abspath(dest)
                    continue
                end
                $f(joinpath(job, file), joinpath(dest, file); kwargs...)
            end
        end
    end
end

is_slurm_job(job::Job) = haskey(job.metadata, :slurmid)

"""
    isrunning(job_dir::String)

Returns whether a job exists in the `job_dir` and if it is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job_dir::String)
    !ispath(joinpath(job_dir, "job.tt")) && return false
    job = load_job(job_dir, -1)
    server = DFC.Server(job)
    n = now()
    if server.scheduler == DFC.Slurm
        return slurm_isrunning(job)
    else
        u = username()
        i = last_running_calculation(job)
        i === nothing && return false
        l = job[i] 
        codeexec = execs(l)[end].exec
        try
            pids = parse.(Int, split(read(`pgrep $codeexec`, String)))
            if isempty(pids)
                return false
            end
            pwd = split(strip(read(`pwdx $(pids[end])`, String)))[end]
            return abspath(pwd) == job.dir
        catch
            return false
        end
    end
    return false
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

"""
    cleanup(job::Job)
    
Cleanup `job.dir` interactively.
"""
function cleanup(job::Job)
    labels = String[]
    paths = String[]
    for v in versions(job)
        vpath = version_dir(job, v)
        s = round(dirsize(vpath) / 1e6; digits = 3)
        push!(labels, "Version $v:  $s Mb")
        push!(paths, vpath)
        opath = joinpath(vpath, Jobs.TEMP_CALC_DIR)
        if ispath(opath)
            s_out = round(dirsize(opath) / 1e6; digits = 3)
            push!(labels, "Version $v/outputs:  $s_out Mb")
            push!(paths, opath)
        end
    end
    menu = MultiSelectMenu(labels)
    choices = request("Select job files to delete:", menu)
    for i in choices
        if ispath(paths[i]) # Could be that outputs was already deleted
            @info "Deleting $(paths[i])"
            rm(paths[i]; recursive = true)
        end
    end
end

function save_metadata(job)
    return jldsave(joinpath(job, ".metadata.jld2"); metadata = job.metadata,
                   version = job.version)
end

timestamp(job::Job) = job.metadata[:timestamp]
timestamp!(job, time) = job.metadata[:timestamp] = time
has_timestamp(job) = haskey(job.metadata, :timestamp)

function clean_dir!(dir::AbstractString)
    for f in readdir(dir)
        if f == Jobs.TEMP_CALC_DIR || f == Jobs.VERSION_DIR_NAME || splitext(f)[end] == ".jl"
            continue
        end
        rm(joinpath(dir, f); recursive = true)
    end
end

exists_job(d::AbstractString) = ispath(d) && ispath(joinpath(d, "job.tt"))

"Finds the output files for each of the calculations of a job, and groups all found data into a dictionary."
function outputdata(job::Job, calculations::Vector{Calculation}; print = true,
                    onlynew = false)
    if DFC.isarchived(job) && ispath(joinpath(job, "results.jld2"))
        return DFC.JLD2.load(joinpath(job, "results.jld2"))["outputdata"]
    end
    datadict = Dict{String, Dict{Symbol,Any}}()
    stime = DFC.starttime(job)
    #TODO Think about storing results.jld2
    for calculation in calculations
        #TODO Would I ever not want to do onlynew ?
        newout = mtime(outpath(calculation)) > stime
        if onlynew && !newout
            continue
        end
        tout = outputdata(calculation; print = print, overwrite = newout)
        if !isempty(tout)
            datadict[calculation.name] = tout
        end
    end
    tmp = tempname() * ".jld2" 
    out = DFC.JLD2.save(tmp, "outputdata", datadict)
    return tmp
end
outputdata(job::Job; kwargs...) = outputdata(job, job.calculations; kwargs...)
