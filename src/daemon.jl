const DAEMON_CONFIG_PATH = config_path("daemon.jld2") 

@with_kw mutable struct Daemon
    port::Int = 9584
    pid::Int  = -1
    jobs_tasks::Dict{String, Task} = Dict{String, Task}()
    query_time::Int = 10
    main_loop = nothing
end

function isalive(d::Daemon)
    d.pid == -1 && return false
    try
        run(pipeline(`ps -p $(d.pid)`, devnull));
        return true
    catch
        return false
    end
end

save(d::Daemon) = JLD2.save(DAEMON_CONFIG_PATH, "daemon", d)

function start(d::Daemon)
    julia_exec   = joinpath(Sys.BINDIR, "julia")
    project_path = Pkg.project().path
    d.pid = getpid(run(`$julia_exec -e "using DFControl; DAEMON = DFControl.server_start($(d.port)); using DaemonMode; serve($(d.port), true)"`, wait=false))
    @info "Daemon started, listening on port $(d.port), with PID $(d.pid)."
    return d
end

function server_start(port::Int)
    d = Daemon(port = port, pid = getpid())
    if ispath(DAEMON_CONFIG_PATH)
        prev_daemon = JLD2.load(DAEMON_CONFIG_PATH)["daemon"]
        for k in keys(prev_daemon)
            if exists_job(k)
                tjob = DFJob(k)
                d.jobs_tasks[k] = run_queue(tjob; sleep_time = d.query_time)
            end
        end
    end
    save(d)
    d.main_loop = @async main_loop(d)
    return d
end

function main_loop(d::Daemon)
    while true
        to_rm = []
        for (k, t) in d.jobs_tasks
            if istaskfailed(t)
                @warn """Workflow in job directory $(k) failed.
                See $(joinpath(k, ".workflow/error.log")) for more info."""
                push!(to_rm, k)
            elseif istaskdone(t)
                @info """Workflow for job directory $(k) done."""
                push!(to_rm, k)
            end
        end
        if !isempty(to_rm)
            for r in to_rm
                pop!(d.jobs_tasks, r, nothing)
            end
            save(d)
        end
        sleep(d.query_time)
    end
end
            
function init_daemon()
    if ispath(DAEMON_CONFIG_PATH)
        d = JLD2.load(DAEMON_CONFIG_PATH)["daemon"]
    else
        d = Daemon()
    end
    if isalive(d)
        return d
    else
        return start(d)
    end
end

kill_daemon() =
    run(`kill $(DAEMON_PID[])`)
        
function run_queue(job::DFJob; sleep_time = 10)
    Threads.@spawn begin
        logger = TeeLogger(MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log"), append=true), Logging.Info),
                           MinLevelLogger(FileLogger(joinpath(job, ".workflow", "info.log"), append=true), Logging.Warn),
                           MinLevelLogger(FileLogger(joinpath(job, ".workflow", "error.log"), append=true), Logging.Error))
        with_logger(logger) do 
            qd = queued_dir(job)  
            queued_steps = readdir(qd)
            fd = finished_dir(job)
            if !ispath(qd) || isempty(queued_steps)
                return
            end
            for f in queued_steps
                step_file = joinpath(qd, f)
                t = include(step_file)
                @eval $(t)($(job))
                while isrunning(job)
                    sleep(sleep_time)
                    yield()
                end
                if !ispath(fd)
                    mkpath(fd)
                end
                mv(step_file, joinpath(fd, f), force=true)
            end
        end
    end
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
        write(joinpath(qd, "$(i + last_step).jl"), @code_string f(job))
    end
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

function submit_queue(job::DFJob, d::Daemon)
    runexpr("""
        DAEMON.jobs_tasks["$(job.local_dir)"] = DFControl.run_queue(DFJob("$(job.local_dir)"))
        DFControl.save(DAEMON)
        """, port=d.port)
end

    
    
    
