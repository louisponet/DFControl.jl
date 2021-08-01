job_path(req) = req.target[findnext(isequal('/'), req.target, 2)+1:end]

save_job(req) = Service.save_job(JSON3.read(req.body, DFJob))
HTTP.@register(ROUTER, "POST", "/jobs", save_job)

get_job(req) = (@show req; Service.load_job(job_path(req), JSON3.read(req.body, Int)))
HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)

job_isrunning(req) = Service.isrunning(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/job_isrunning/*", job_isrunning)

registered_jobs(req) = Service.registered_jobs(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/registered_jobs/*", registered_jobs)

job_versions(req) = Service.job_versions(job_path(req)) 
HTTP.@register(ROUTER, "GET", "/job_versions/*", job_versions)

last_running_calculation(req) = Service.last_running_calculation(JSON3.read(req.body, DFJob)) 
HTTP.@register(ROUTER, "GET", "/last_running_calculation", last_running_calculation)
