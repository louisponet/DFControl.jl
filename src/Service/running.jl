function main_loop()
    mkpath(PENDING_JOBS_DIR)
    running_jobs = ispath(RUNNING_JOBS_FILE) ? readlines(RUNNING_JOBS_FILE) : String[]
    job_dirs_procs = Dict{String, Tuple{Int, Future}}()
    for j in running_jobs
        if DFControl.exists_job(j)
            tjob = DFJob(j)
            job_dirs_procs[j] = spawn_worker(tjob)
        end
    end
    while true
        handle_workflow_runners!(job_dirs_procs)
        handle_job_submission!(job_dirs_procs)
        sleep(SLEEP_TIME)
    end
end

function spawn_worker(job::DFJob)
    if ispath(joinpath(job, ".workflow/environment.jld2"))
        env_dat = DFControl.JLD2.load(joinpath(job, ".workflow/environment.jld2"))
        proc = addprocs(1; exeflags = "--project=$(env_dat["project"])")[1]
        using_expr = Base.Meta.parse("""using $(join(env_dat["modules"], ", "))""")
        ctx_path = joinpath(job, ".workflow/ctx.jld2")

        Distributed.remotecall_eval(Distributed.Main, proc, using_expr)
        f = remotecall(DFControl.run_queue, proc, job, DFControl.JLD2.load(ctx_path)["ctx"];
                       sleep_time = SLEEP_TIME)
    else
        proc = addprocs(1; exeflags = "--project=$(Base.current_project())")[1]
        Distributed.remotecall_eval(Distributed.Main, proc, :(using DFControl))
        Distributed.remotecall_eval(Distributed.Main, proc, :(DFControl.global_logger(DFControl.FileLogger(joinpath($(job.dir), "submission.log"); append = true))))
        f = Distributed.remotecall(DFControl.submit, proc, job)
    end
    return proc, f
end

function handle_workflow_runners!(job_dirs_procs)
    to_rm = String[]
    for (k, t) in job_dirs_procs
        if isready(t[2])
            try
                t_ = fetch(t[2])
                @info """Workflow for job directory $(k) done."""
                push!(to_rm, k)
            catch e
                @warn """Workflow in job directory $(k) failed.
                See $(joinpath(k, ".workflow/error.log")) for more info."""
                push!(to_rm, k)
            end
        end
    end
    if !isempty(to_rm)
        for r in to_rm
            rmprocs(job_dirs_procs[r][1])
        end
        pop!.((job_dirs_procs,), to_rm)
        save_running_jobs(job_dirs_procs)
    end
end

save_running_jobs(job_dirs_procs) = DFControl.writelines(RUNNING_JOBS_FILE, keys(job_dirs_procs))

# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(job_dirs_procs)
    pending_job_submissions = readdir(PENDING_JOBS_DIR)
    for j in pending_job_submissions
        dat = DFControl.JLD2.load(joinpath(PENDING_JOBS_DIR, j))
        job = dat["job"]
        mkpath(job.dir)
        job.server = "localhost"
        for (fname, contents) in dat["files"]
            write(joinpath(job, fname), contents)
        end
        for a in atoms(job)
            a.pseudo.dir = job.dir
        end
        job_dirs_procs[job.dir] = spawn_worker(job)
        rm(joinpath(PENDING_JOBS_DIR, j))
    end
end

function run_queue(job::DFJob, ctx::Dict; sleep_time = 10)
    logger = workflow_logger(job)
    with_logger(logger) do
        qd = queued_dir(job)
        queued_steps = readdir(qd)
        order = sortperm(parse.(Int, getindex.(splitext.(queued_steps), 1)))
        fd = finished_dir(job)
        if !ispath(qd) || isempty(queued_steps)
            return
        end
        if !ispath(fd)
            mkpath(fd)
        else
            for f in readdir(fd)
                rm(joinpath(fd, f))
            end
        end

        for f in queued_steps[order]
            step_file = joinpath(qd, f)
            t = include(step_file)
            @eval $(t)($(job), $(ctx))
            DFControl.JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", ctx)
            while isrunning(job)
                sleep(sleep_time)
            end
            mv(step_file, joinpath(fd, f); force = true)
        end
    end
    return true
end

