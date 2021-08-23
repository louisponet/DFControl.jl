const JOB_REGISTRY = (in_progress = DFC.config_path("in_progress_jobs.txt"), archived = DFC.config_path("archived_jobs.txt"))

all_known_jobs() = vcat(readlines(JOB_REGISTRY[1]), readlines(JOB_REGISTRY[2]))

function cleanup_job_registry!()
    in_progress, archived = readlines(JOB_REGISTRY[1]), readlines(JOB_REGISTRY[2])
    for (ir, REG) in enumerate((in_progress, archived))
        stale_ids = findall(x -> !ispath(x) || !ispath(joinpath(x, "job.tt")), REG)
        if stale_ids !== nothing
            jobs_to_remove = REG[stale_ids]
            deleteat!(REG, stale_ids)
            unique!(REG)
        end
        write(JOB_REGISTRY[ir], join(REG, "\n"))
    end
    return (; in_progress, archived)
end

function maybe_register_job(abspath::AbstractString)
    if ispath(abspath)
        n = occursin(".archived", abspath) ? JOB_REGISTRY.archived :
              JOB_REGISTRY.in_progress
        REG = readlines(n)
        jid = findfirst(isequal(abspath), REG)
        if jid === nothing
            push!(REG, abspath)
            write(n, join(REG, "\n"))
            # cleanup_job_registry!()
        end
    end
end
maybe_register_job(job::Job) = maybe_register_job(main_job_dir(job.dir))

"""
    registered_jobs(fuzzy::AbstractString = "")

Lists all the known [`Jobs`](@ref Job) directories that contain `fuzzy`.
Intended to be used as:
```julia
job_dirs = registered_jobs("NiO")
job = Job(job_dirs[1])
```
"""
function registered_jobs(fuzzy::AbstractString = "")
    choices = filter(x -> occursin(fuzzy, x), vcat(cleanup_job_registry!()...))
    timestamps = timestamp.(choices)
    sort_ids = sortperm(timestamps)
    return [(choices[i], timestamps[i]) for i in sort_ids]
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

Returns a `Vector` of pairs with all archived [`Jobs`](@ref Job) whose directory contains the fuzzy as the first, and their description as the second item of the pair.
"""
function archived_jobs(fuzzy::AbstractString = "")
    jobs = filter(x -> occursin(fuzzy, x), cleanup_job_registry!().archived)
    return [j => ispath(joinpath(j, "description.txt")) ?
                 read(joinpath(j, "description.txt"), String) : "" for j in jobs]
end

"""
    load_jobs(fuzzy::AbstractString)

Loads all the known [`Jobs`](@ref Job) whose `dir` contains `fuzzy`.
"""
load_jobs(fuzzy::AbstractString) = map(x -> Job(x[1]), registered_jobs(fuzzy))

"""
    load_running_jobs(fuzzy::AbstractString = "")

Loads all [`Jobs`](@ref Job) that are currently running.
"""
function load_running_jobs(fuzzy::AbstractString = "")
    jobs = load_jobs(fuzzy)
    return filter(isrunning, jobs)
end

"""
    running_jobs(fuzzy::AbstractString = "")

Returns all the paths to running ['Jobs'](@ref Job).
"""
running_jobs(fuzzy::AbstractString = "") = getfield.(load_running_jobs(fuzzy), :dir)
