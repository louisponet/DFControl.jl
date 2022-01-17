save_job(req) = Service.save(path(req), JSON3.read(req.body, Dict{String, String}))
HTTP.@register(ROUTER, "POST", "/jobs", save_job)

submit_job(req) = Service.submit(path(req))
HTTP.@register(ROUTER, "PUT", "/jobs/*", submit_job)

get_job(req) = Service.load_job(path(req))
HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)
HTTP.@register(ROUTER, "GET", "/jobs/", get_job)

job_state(req) = Service.state(path(req))
HTTP.@register(ROUTER, "GET", "/job_state/*", job_state)

registered_jobs(req) = Service.registered_jobs(path(req))
HTTP.@register(ROUTER, "GET", "/registered_jobs/*", registered_jobs)
HTTP.@register(ROUTER, "GET", "/registered_jobs/", registered_jobs)

job_versions(req) = Service.job_versions(path(req))
HTTP.@register(ROUTER, "GET", "/job_versions/*", job_versions)

rm_version!(req) = Service.rm_version!(path(req), JSON3.read(req.body, Int))
HTTP.@register(ROUTER, "PUT", "/rm_version/*", rm_version!)

last_running_calculation(req) = Service.last_running_calculation(path(req))
HTTP.@register(ROUTER, "GET", "/last_running_calculation/*", last_running_calculation)

outputdata(req) = Service.outputdata(path(req), JSON3.read(req.body, Vector{String}))
HTTP.@register(ROUTER, "GET", "/outputdata/*", outputdata)

running_jobs(req) = Service.running_jobs(path(req))
HTTP.@register(ROUTER, "GET", "/running_jobs/", running_jobs)
HTTP.@register(ROUTER, "GET", "/running_jobs/*", running_jobs)

abort_job(req) = Service.abort(path(req))
HTTP.@register(ROUTER, "GET", "/abort/*", abort_job)
