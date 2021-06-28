const DAEMON_PORT = Ref(9584)
const DAEMON_PID  = Ref(0)
const DAEMON_CONFIG_PATH = config_path("daemon.jld2") 

function init_daemon()
    try
        c = JLD2.load(DAEMON_CONFIG_PATH)
        DAEMON_PORT[] = c["PORT"]
        DAEMON_PID[] = c["PID"]
        runexpr("5+1", port=DAEMON_PORT[])
    catch
        julia_exec = joinpath(Sys.BINDIR, "julia")
        project_path = Pkg.project().path
        DAEMON_PID[] = getpid(run(`$julia_exec -t auto --project=$project_path -e "using DaemonMode; serve($(DAEMON_PORT[]), async=true, threaded=true, shared=true)"`, wait=false))
        @info "Daemon created, listening on port $(DAEMON_PORT[]), with PID $(DAEMON_PID[])."
        JLD2.save(DAEMON_CONFIG_PATH, "PID", DAEMON_PID[], "PORT", DAEMON_PORT[])
    end
    runexpr("using DFControl; jobs_tasks = Pair{String,Task}[]", port=DAEMON_PORT[])
end
        
function run_queue(job::DFJob; sleep_time = 10)
    qd = queued_dir(job)  
    queued_steps = readdir(qd)
    fd = finished_dir(job)
    if !ispath(qd) || isempty(queued_steps)
        return
    end
    Threads.@spawn begin
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

