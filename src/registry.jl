const JOB_REGISTRY = (in_progress = String[], archived = String[])

all_known_jobs() = vcat(JOB_REGISTRY[1], JOB_REGISTRY[2]) 

function init_job_registry()
    append!(JOB_REGISTRY.in_progress, readlines(config_path("in_progress_jobs.txt")))
    append!(JOB_REGISTRY.archived, readlines(config_path("archived_jobs.txt")))
end

function write_job_registry()
    open(config_path("in_progress_jobs.txt"), "w") do f
        for j in JOB_REGISTRY.in_progress
            write(f, "$j\n")
        end
    end
    open(config_path("archived_jobs.txt"), "w") do f
        for j in JOB_REGISTRY.archived
            write(f, "$j\n")
        end
    end
end

function cleanup_job_registry!(; print = true)
    for REG in JOB_REGISTRY
        stale_ids = findall(x -> !ispath(x) || !ispath(joinpath(x, "job.tt")), REG)
        if stale_ids !== nothing
            jobs_to_remove = REG[stale_ids]
            message = "Removing $(length(jobs_to_remove)) stale jobs (job folder removed) from the registry:\n"
            for j in jobs_to_remove
                message *= "\t$j\n"
            end
            print && @warn message
            deleteat!(REG, stale_ids)
        end
    end
    return write_job_registry()
end

function maybe_register_job(abspath::String)
    if ispath(abspath)
        REG = occursin(".archived", abspath) ? JOB_REGISTRY.archived : JOB_REGISTRY.in_progress
        jid = findfirst(isequal(abspath), REG)
        if jid === nothing
            push!(REG, abspath)
            write_job_registry()
        end
    end
end
maybe_register_job(job::DFJob) = maybe_register_job(job.local_dir)

"""
    registered_jobs(fuzzy::AbstractString = "")

Lists all the known [`DFJobs`](@ref DFJob) directories that contain `fuzzy`.
Intended to be used as:
```julia
job_dirs = registered_jobs("NiO")
job = DFJob(job_dirs[1])
```
"""
function registered_jobs(fuzzy::AbstractString = "")
    cleanup_job_registry!(; print = false)
    choices = filter(x -> occursin(fuzzy, x), all_known_jobs()) 
    timestamps = timestamp.(choices)
    sort_ids = sortperm(timestamps)
    return choices[sort_ids]
end

function timestamp(jobdir::AbstractString)
    if ispath(joinpath(jobdir, ".metadata.jld2"))
        md = load(joinpath(jobdir, ".metadata.jld2"))["metadata"]
        return get(md, :timestamp, DateTime(0))
    else
        return DateTime(0)
    end
end

"""
    archived_jobs(fuzzy::AbstractString = "")

Returns a `Vector` of pairs with all archived [`DFJobs`](@ref DFJob) whose directory contains the fuzzy as the first, and their description as the second item of the pair.
"""
function archived_jobs(fuzzy::AbstractString = "")
    cleanup_job_registry!(; print = false)
    jobs = filter(x -> occursin(fuzzy, x), JOB_REGISTRY.archived)
    return Dict([j => ispath(joinpath(j, "description.txt")) ? read(joinpath(j, "description.txt")) : "" for j in jobs])
end

"""
    load_jobs(fuzzy::AbstractString)

Loads all the known [`DFJobs`](@ref DFJob) whose `local_dir` contains `fuzzy`.
"""
load_jobs(fuzzy::AbstractString) = DFJob.(registered_jobs(fuzzy))

function request_job(job_dir::String)
    matching_jobs = reverse(registered_jobs(job_dir))
    choices = ["$j -- $t" for (j, t) in zip(matching_jobs, timestamp.(matching_jobs))]

    if length(matching_jobs) == 1
        return matching_jobs[1]
    elseif isdefined(Base, :active_repl)
        menu = RadioMenu(choices)
        choice = request("Multiple matching jobs were found, choose one:", menu)
        if choice != -1
            return matching_jobs[choice]
        else
            return nothing
        end
    else
        err_msg = "No concrete job for $job_dir found, closest matches are:"
        for m in choices
            err_msg = join([err_msg, m], "\n")
        end
        @error err_msg
    end
end

"""
    load_running_jobs(fuzzy::AbstractString = "")

Loads all [`DFJobs`](@ref DFJob) that are currently running.
"""
function load_running_jobs(fuzzy::AbstractString = "")
    jobs = load_jobs(fuzzy)
    return filter(isrunning, jobs)
end

