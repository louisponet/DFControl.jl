module Service
    using ..DFControl
    using Distributed, Pkg, CodeTracking
    using ..DFControl: config_path
    # const DAEMON_CONFIG_PATH = config_path("daemon.jld2")
    # delete_daemon_config!() = rm(DAEMON_CONFIG_PATH)

    const RUNNING_JOBS_FILE = config_path("running_jobs.txt")
    const PENDING_JOBS_DIR = config_path("pending_jobs")
    const SERVICE_LOG = config_path("daemon.log")
    const SLEEP_TIME = 10.0
    
    function launch_service(s::DFControl.Server)
        julia_exec   = joinpath(Sys.BINDIR, "julia")
        project_path = Pkg.project().path
        cmd = Cmd(`$julia_exec -t auto --project=$project_path -e "using DFControl; DFControl.Resource.run($s.port)"`; detach = true)
        proc = DFControl.establish_connection(s)
        
        p   = remotecall(run, cmd; wait = false)
        @info "Daemon on Server $(s.name) started, listening on port $(s.port)."
    end

    daemon_logger() = FileLogger(SERVICE_LOG; append = true)
    function workflow_logger(job::DFJob)
        return TeeLogger(MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                                   append = true), Logging.Info),
                         MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                                   append = true), Logging.Warn),
                         MinLevelLogger(FileLogger(joinpath(job, ".workflow", "error.log");
                                                   append = true), Logging.Error))
    end

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
            Distributed.remotecall_eval(Distributed.Main, proc, :(DFControl.global_logger(DFControl.FileLogger(joinpath($(job.local_dir), "submission.log"); append = true))))
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
    load_running_jobs() = readlines(RUNNING_JOBS_FILE)

    # Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
    # This means we need to turn job.server_dir into job.local_dir
    # Additional files are packaged with the job
    function handle_job_submission!(job_dirs_procs)
        pending_job_submissions = readdir(PENDING_JOBS_DIR)
        for j in pending_job_submissions
            dat = DFControl.JLD2.load(joinpath(PENDING_JOBS_DIR, j))
            job = dat["job"]
            set_localdir!(job, job.server_dir)
            mkpath(job.local_dir)
            job.server_dir = ""
            job.server = "localhost"
            for (fname, contents) in dat["files"]
                write(joinpath(job, fname), contents)
            end
            for a in atoms(job)
                a.pseudo.dir = job.local_dir
            end
            job_dirs_procs[job.local_dir] = spawn_worker(job)
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

    queued_dir(job::DFJob, args...) = joinpath(job, ".workflow/queued", args...)
    finished_dir(job::DFJob, args...) = joinpath(job, ".workflow/finished", args...)

    function queue_steps(job::DFJob, funcs)
        qd = queued_dir(job)
        if !ispath(qd)
            mkpath(qd)
        end
        prev_steps = readdir(qd)
        last_step = !isempty(prev_steps) ? parse(Int, splitext(prev_steps[end])[1]) : 0
        for (i, f) in enumerate(funcs)
            write(joinpath(qd, "$(i + last_step).jl"), @code_string f(job, Dict{Symbol,Any}()))
        end
    end

    function write_workflow_files(job::DFJob)
        all = string.(values(Base.loaded_modules))
        valid = String[]
        ks = keys(Pkg.project().dependencies)
        for p in all
            if p ∈ ks
                push!(valid, p)
            end
        end
        DFControl.JLD2.save(joinpath(job, ".workflow/environment.jld2"), "modules", valid, "project",
                  Base.current_project())
        return DFControl.JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", Dict{Symbol,Any}())
    end


    function clear_queue!(job::DFJob)
        qd = queued_dir(job)
        if !ispath(qd)
            return
        end
        for f in readdir(qd)
            rm(joinpath(qd, f))
        end
    end

    queued(job::DFJob) = readdir(queued_dir(job))
    finished(job::DFJob) = readdir(finished_dir(job))

    # function submit_workflow(job::DFJob, funcs, d::Daemon = init_daemon())
    #     queue_steps(job, funcs)
    #     write_workflow_files(job)
    #     while !is_started(d)
    #         sleep(1)
    #     end
    #     return runexpr("""
    #                DFControl.spawn_worker(DAEMON, DFJob("$(job.local_dir)"))
    #                DFControl.save(DAEMON)
    #                """; port = d.port)
    # end

    mods_test() = Base.loaded_modules
    save_job(j::DFJob) = DFControl.save(j)
    get_job(p::AbstractString) = DFJob(p)
end
