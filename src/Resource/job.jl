function job_path(req)
    id = findnext(isequal('/'), req.target, 2)
    if length(req.target) < id + 1
        return ""
    else
        return req.target[id+1:end]
    end
end

save_job(req) = Service.save(JSON3.read(req.body, DFJob))
HTTP.@register(ROUTER, "POST", "/jobs", save_job)

submit_job(req) = Service.submit(job_path(req))
HTTP.@register(ROUTER, "PUT", "/jobs/*", submit_job)

get_job(req) = Service.load_job(job_path(req), JSON3.read(req.body, Int))
HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)
HTTP.@register(ROUTER, "GET", "/jobs/", get_job)

job_isrunning(req) = Service.isrunning(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/job_isrunning/*", job_isrunning)

registered_jobs(req) = Service.registered_jobs(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/registered_jobs/*", registered_jobs)
HTTP.@register(ROUTER, "GET", "/registered_jobs/", registered_jobs)

job_versions(req) = Service.job_versions(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/job_versions/*", job_versions)

last_running_calculation(req) = Service.last_running_calculation(JSON3.read(req.body, DFJob)) 
HTTP.@register(ROUTER, "GET", "/last_running_calculation", last_running_calculation)

outputdata(req) = Service.outputdata(JSON3.read(req.body, DFJob)) 
HTTP.@register(ROUTER, "GET", "/outputdata", outputdata)

