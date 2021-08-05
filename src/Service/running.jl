function main_loop()
    running_jobs = ispath(RUNNING_JOBS_FILE) ? readlines(RUNNING_JOBS_FILE) : String[]
    job_dirs_procs = Dict{String,Tuple{Int,Future}}()
    for j in running_jobs
        if DFControl.exists_job(j)
            tjob = load_job(j)
            job_dirs_procs[j] = spawn_worker(tjob)
        end
    end
    while true
        handle_workflow_runners!(job_dirs_procs)
        handle_job_submission!(job_dirs_procs)
        sleep(SLEEP_TIME)
    end
end

function spawn_worker(job::Job)
    if ispath(joinpath(job, ".workflow/environment.jld2"))
        env_dat = JLD2.load(joinpath(job, ".workflow/environment.jld2"))
        proc = addprocs(1; exeflags = "--project=$(env_dat["project"])")[1]
        using_expr = Base.Meta.parse("""using $(join(env_dat["modules"], ", "))""")
        ctx_path = joinpath(job, ".workflow/ctx.jld2")

        Distributed.remotecall_eval(Distributed.Main, proc, using_expr)
        f = remotecall(DFControl.run_queue, proc, job, JLD2.load(ctx_path)["ctx"];
                       sleep_time = SLEEP_TIME)
        return proc, f
    else
        s = DFC.Server("localhost")
        proc = addprocs(1)[1]
        remotecall_fetch(cd, proc, job.dir)
        if s.scheduler == Servers.Bash
            cmd = Cmd(`bash job.tt`)
        else
            cmd = Cmd(`sbatch job.tt`)
        end
        f = remotecall(run, proc, cmd)
        return proc, f
    end
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

save_running_jobs(job_dirs_procs) = write(RUNNING_JOBS_FILE, join(keys(job_dirs_procs), "\n"))

# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(job_dirs_procs)
    if ispath(PENDING_JOBS_FILE)
        pending_job_submissions = readlines(PENDING_JOBS_FILE)
        for j in pending_job_submissions
            job = load_job(j)
            job_dirs_procs[job.dir] = spawn_worker(job)
        end
        write(PENDING_JOBS_FILE, "")
    end
end

function run_queue(job::Job, ctx::Dict; sleep_time = 10)
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
            JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", ctx)
            while isrunning(job)
                sleep(sleep_time)
            end
            mv(step_file, joinpath(fd, f); force = true)
        end
    end
    return true
end
