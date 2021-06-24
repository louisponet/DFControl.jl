const JOB_REGISTRY = String[]

init_job_registry() = append!(JOB_REGISTRY, readlines(config_path("job_registry.txt")))

function write_job_registry()
    open(config_path("job_registry.txt"), "w") do f
        for j in JOB_REGISTRY
            write(f, "$j\n")
        end
    end
end

function cleanup_job_registry(; print = true)
    stale_ids = findall(x -> !ispath(x) || !ispath(joinpath(x, "job.tt")), JOB_REGISTRY)
    if stale_ids !== nothing
        jobs_to_remove = JOB_REGISTRY[stale_ids]
        message = "Removing $(length(jobs_to_remove)) stale jobs (job folder removed) from the registry:\n"
        for j in jobs_to_remove
            message *= "\t$j\n"
        end
        print && @warn message
        deleteat!(JOB_REGISTRY, stale_ids)
    end
    return write_job_registry()
end

function maybe_register_job(abspath::String)
    if ispath(abspath)
        jid = findfirst(isequal(abspath), JOB_REGISTRY)
        if jid === nothing
            push!(JOB_REGISTRY, abspath)
            write_job_registry()
        end
    end
end
maybe_register_job(job::DFJob) = maybe_register_job(job.local_dir)

"""
    registered_jobs(fuzzy::AbstractString)

Lists all the known [`DFJob](@ref) directories that contain `fuzzy`.
Intended to be used as:
```julia
job_dirs = registered_jobs("NiO")
job = DFJob(job_dirs[1])
```
"""
registered_jobs(fuzzy::AbstractString) = filter(x -> occursin(fuzzy, x), JOB_REGISTRY)

"""
    load_jobs(fuzzy::AbstractString)

Loads all the known [`DFJobs`](@ref DFJob) whose `local_dir` contains `fuzzy`.
"""
load_jobs(fuzzy::AbstractString) = DFJob.(registered_jobs(fuzzy))

function request_job(job_dir::String)
    cleanup_job_registry(; print = false)
    function timestamp(jobdir)
        if ispath(joinpath(jobdir, ".metadata.jld2"))
            md = load(joinpath(jobdir, ".metadata.jld2"))["metadata"]
            return get(md, :timestamp, DateTime(0))
        else
            return DateTime(0)
        end
    end
    matching_jobs = registered_jobs(job_dir)
    timestamps = timestamp.(matching_jobs)
    sort_ids = sortperm(timestamps, rev=true)

    choices = ["$j -- $t" for (j, t) in zip(matching_jobs[sort_ids], timestamps[sort_ids])]
    
    if length(matching_jobs) == 1
        return matching_jobs[1]
    elseif isdefined(Base, :active_repl)
        menu = RadioMenu(choices)
        choice = request("Multiple matching jobs were found, choose one:", menu)
        if choice != -1
            return matching_jobs[sort_ids[choice]]
        else
            return nothing
        end
    else
        err_msg = "No concrete job for $job_dir found, closest matches are:"
        for m in choices
            err_msg = join(err_msg, m, "\n")
        end
        @error err_msg
    end
end
