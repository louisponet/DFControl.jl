function job_path(req)
    id = findnext(isequal('/'), req.target, 2)
    if length(req.target) < id + 1
        return ""
    else
        return req.target[id+1:end]
    end
end

save_job(req) = Service.save(job_path(req), JSON3.read(req.body, Dict{String, String}))
HTTP.@register(ROUTER, "POST", "/jobs", save_job)

submit_job(req) = Service.submit(job_path(req))
HTTP.@register(ROUTER, "PUT", "/jobs/*", submit_job)

get_job(req) = Service.load_job(job_path(req))
HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)
HTTP.@register(ROUTER, "GET", "/jobs/", get_job)

job_state(req) = Service.state(job_path(req))
HTTP.@register(ROUTER, "GET", "/job_state/*", job_state)

function job_submission_time(req)
    scriptpath = joinpath(job_path(req), "job.tt")
    if ispath(scriptpath)
        return mtime(scriptpath)
    else
        return 0
    end
end
HTTP.@register(ROUTER, "GET", "/job_submission_time/*", job_submission_time)

registered_jobs(req) = Service.registered_jobs(job_path(req))
HTTP.@register(ROUTER, "GET", "/registered_jobs/*", registered_jobs)
HTTP.@register(ROUTER, "GET", "/registered_jobs/", registered_jobs)

job_versions(req) = Service.job_versions(job_path(req))
HTTP.@register(ROUTER, "GET", "/job_versions/*", job_versions)

rm_version!(req) = Service.rm_version!(job_path(req), JSON3.read(req.body, Int))
HTTP.@register(ROUTER, "PUT", "/rm_version/*", rm_version!)

last_running_calculation(req) = Service.last_running_calculation(job_path(req))
HTTP.@register(ROUTER, "GET", "/last_running_calculation/*", last_running_calculation)

outputdata(req) = Service.outputdata(job_path(req), JSON3.read(req.body, Vector{String}))
HTTP.@register(ROUTER, "GET", "/outputdata/*", outputdata)

running_jobs(req) = Service.running_jobs(job_path(req))
HTTP.@register(ROUTER, "GET", "/running_jobs/", running_jobs)
HTTP.@register(ROUTER, "GET", "/running_jobs/*", running_jobs)
