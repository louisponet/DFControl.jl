const JOB_REGISTRY = String[]

init_job_registry() =
    append!(JOB_REGISTRY, readlines(config_path("job_registry.txt")))
    
function write_job_registry()
    open(config_path("job_registry.txt"), "w") do f
        for j in JOB_REGISTRY
            write(f, "$j\n")
        end
    end
end

function cleanup_job_registry()
    stale_ids = findall(x -> !ispath(x) || !ispath(joinpath(x, "job.tt")), JOB_REGISTRY)
    if stale_ids !== nothing
        jobs_to_remove = JOB_REGISTRY[stale_ids] 
        message = "Removing $(length(jobs_to_remove)) stale jobs (job folder removed) from the registry:\n"
        for j in jobs_to_remove
            message *= "\t$j\n"
        end
        @warn message
        deleteat!(JOB_REGISTRY, stale_ids)
    end
    write_job_registry()
end

function maybe_register_job(abspath::String)
    if !ispath(abspath)
        push!(JOB_REGISTRY, abspath)
    else
        jid = findfirst(isequal(abspath), JOB_REGISTRY)
        if jid === nothing
            push!(JOB_REGISTRY, abspath)
        end
    end
    write_job_registry()
end
