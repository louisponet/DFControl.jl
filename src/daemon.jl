const DAEMON_CONFIG_PATH = config_path("daemon.jld2")
delete_daemon_config!() = rm(DAEMON_CONFIG_PATH)

@with_kw mutable struct Daemon
    port::Int = 9584
    pid::Int = -1
    job_dirs_procs::Dict = Dict{String,Tuple{Int,Future}}()
    query_time::Int = 10
    started::Bool = false
    main_loop = nothing
end

function isalive(d::Daemon)
    d.pid == -1 && return false
    try
        run(pipeline(`ps -p $(d.pid)`, devnull))
        return true
    catch
        return false
    end
end

function save(d::Daemon)
    return JLD2.save(DAEMON_CONFIG_PATH, "port", d.port, "pid", d.pid, "query_time",
                     d.query_time, "job_dirs_procs", d.job_dirs_procs, "started", d.started)
end

function start(d::Daemon)
    julia_exec   = joinpath(Sys.BINDIR, "julia")
    project_path = Pkg.project().path
    p            = run(Cmd(`$julia_exec -t auto --project=$project_path -e "using DFControl; DAEMON = DFControl.server_start($(d.port)); using DFControl.DaemonMode; serve($(d.port), true)"`))#; detach = true); wait = false)
    d.pid        = getpid(p)
    @info "Daemon started, listening on port $(d.port), with PID $(d.pid).#"
    d.started = false
    save(d)
    return d
end

is_started(d::Daemon) = JLD2.load(DAEMON_CONFIG_PATH)["started"]

function server_start(port::Int)
    with_logger(daemon_logger()) do
        d = Daemon(; port = port, pid = getpid())
        if ispath(DAEMON_CONFIG_PATH)
            prev_daemon = Daemon(DAEMON_CONFIG_PATH)
            d.query_time = prev_daemon.query_time
            for j in keys(prev_daemon.job_dirs_procs)
                if exists_job(j)
                    tjob = DFJob(j)
                    spawn_worker(d, tjob)
                end
            end
        end
        d.main_loop = @async main_loop(d)
        d.started = true
        save(d)
        return d
    end
end

daemon_logger() = FileLogger(config_path("daemon.log"); append = true)
function workflow_logger(job::DFJob)
    return TeeLogger(MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                               append = true), Logging.Info),
                     MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log");
                                               append = true), Logging.Warn),
                     MinLevelLogger(FileLogger(joinpath(job, ".workflow", "error.log");
                                               append = true), Logging.Error))
end

function init_daemon()
    if ispath(DAEMON_CONFIG_PATH)
        data = JLD2.load(DAEMON_CONFIG_PATH)
        d = Daemon(; port = data["port"], pid = data["pid"], query_time = data["query_time"],
                  job_dirs_procs = data["job_dirs_procs"], started = data["started"])
    else
        d = Daemon()
    end
    # create necessary daemon dirs (if they didn't exist yet)
    mkpath(config_path("pending_jobs"))
    if isalive(d)
        @info "Daemon found running at port $(d.port)."
        return d
    else
        return start(d)
    end
end

kill_daemon(d::Daemon) = run(`kill $(d.pid)`)

function main_loop(d::Daemon)
    while true
        handle_workflow_runners!(d)
        handle_job_submission(d)
        sleep(d.query_time)
    end
end

function handle_workflow_runners!(d::Daemon)
    to_rm = String[]
    for (k, t) in d.job_dirs_procs
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
            rmprocs(d.job_dirs_procs[r][1])
        end
        pop!.((d.job_dirs_procs,), to_rm)
        save(d)
    end
end 

# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# This means we need to turn job.server_dir into job.local_dir
# Additional files are packaged with the job
function handle_job_submission(d::Daemon)
    pending_job_submissions = readdir(config_dir("pending_jobs")) 
    for j in pending_job_submissions
        dat = JLD2.load(j)
        job = dat["job"]
        set_localdir!(job, job.server_dir)
        job.server_dir = ""
        job.server = "localhost"
        for (fname, contents) in dat["files"]
            write(fname, contents)
        end
        for a in atoms(job)
            a.pseudo.dir = ""
        end
        # TODO For now just submit
        submit(job)
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
            JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", ctx)
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
    JLD2.save(joinpath(job, ".workflow/environment.jld2"), "modules", valid, "project",
              Base.current_project())
    return JLD2.save(joinpath(job, ".workflow/ctx.jld2"), "ctx", Dict{Symbol,Any}())
end

function spawn_worker(d::Daemon, job::DFJob)
    env_dat = JLD2.load(joinpath(job, ".workflow/environment.jld2"))
    proc = addprocs(1; exeflags = "--project=$(env_dat["project"])")[1]
    using_expr = Base.Meta.parse("""using $(join(env_dat["modules"], ", "))""")
    ctx_path = joinpath(job, ".workflow/ctx.jld2")

    Distributed.remotecall_eval(Distributed.Main, proc, using_expr)
    f = remotecall(DFControl.run_queue, proc, job, DFControl.JLD2.load(ctx_path)["ctx"];
                   sleep_time = d.query_time)
    d.job_dirs_procs[job.local_dir] = (proc, f)
    return save(d)
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

function submit_workflow(job::DFJob, funcs, d::Daemon = init_daemon())
    queue_steps(job, funcs)
    write_workflow_files(job)
    while !is_started(d)
        sleep(1)
    end
    return runexpr("""
               DFControl.spawn_worker(DAEMON, DFJob("$(job.local_dir)"))
               DFControl.save(DAEMON)
               """; port = d.port)
end

mods_test() = Base.loaded_modules
