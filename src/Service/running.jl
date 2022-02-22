const JOB_QUEUE = Ref(Dict{String, Tuple{Int, Jobs.JobState}}())

function main_loop(s::Server)
    job_dirs_procs = Dict()
    queue!(JOB_QUEUE[], s.scheduler, true)
    @info (timestamp = Dates.now(), username = ENV["USER"], host = gethostname(), pid=getpid())
    for (j, info) in JOB_QUEUE[]
        if info[2] == Jobs.Running 
            job_dirs_procs[j] = spawn_task(s, j)
        end
    end

    # Used to identify if multiple servers are running in order to selfdestruct 
    log_mtimes = mtime.(joinpath.((config_path("logs/runtimes/"),), readdir(config_path("logs/runtimes/"))))
    
    while true
        queue!(JOB_QUEUE[], s.scheduler, false)
        handle_workflow_runners!(job_dirs_procs)
        handle_job_submission!(s, job_dirs_procs)
        # handle_workflow_submission!(s, job_dirs_procs)
        monitor_issues(log_mtimes)
        print_log(s, job_dirs_procs)
        if ispath(config_path("self_destruct"))
            @info (timestamp = Dates.now(), message = "self_destruct found, self destructing...")
            for (k, t) in job_dirs_procs
                if t isa Task
                    continue
                else
                    run(`pkill $t`)
                end
            end
            @info (timestamp = Dates.now(), message = "self destructing complete")
            exit()
        end
            
        sleep(SLEEP_TIME)
    end
end

function print_log(s, job_dirs_procs)
    if length(job_dirs_procs) > 0
        ntasks = sum(x-> x isa Task ? 1 : 0, values(job_dirs_procs))
        nprocs = length(job_dirs_procs) - ntasks
    else
        ntasks = 0
        nprocs = 0
    end
    @info (timestamp = Dates.now(), ntasks = ntasks, nprocs = nprocs)
end

function monitor_issues(log_mtimes)
    new_mtimes = mtime.(joinpath.((config_path("logs/runtimes"),), readdir(config_path("logs/runtimes"))))
    if length(new_mtimes) != length(log_mtimes)
        @error "More Server logs got created signalling a server was started while a previous was running."
        touch(config_path("self_destruct"))
    end
    ndiff = length(filter(x -> log_mtimes[x] != new_mtimes[x], 1:length(log_mtimes)))
    if ndiff > 1
        @error "More Server logs modification times differed than 1."
        touch(config_path("self_destruct"))
    end
end
   
function handle_workflow_runners!(job_dirs_procs)
    to_rm = String[]
    for (k, t) in job_dirs_procs
        if t isa Task
            if istaskdone(t)
                if istaskfailed(t)
                    @error "Task failed with $(t.result)"
                else
                    t_ = fetch(t)
                    @info (timestamp = Dates.now(), jobdir = k, state=Jobs.Completed)
                end
                push!(to_rm, k)
            end
        else
            if t.exitcode >= 0
                push!(to_rm, k)
            end
        end
    end
    if !isempty(to_rm)
        pop!.((job_dirs_procs,), to_rm)
    end
end

function spawn_task(s::Server, jobdir)
    id = Servers.submit(s.scheduler, jobdir)
    @info (timestamp = Dates.now(), jobdir = jobdir, jobid = id, state = Jobs.Submitted)
    JOB_QUEUE[][jobdir] = (id, Jobs.Submitted)
    return Threads.@spawn begin
        with_logger(job_logger(id)) do
            @info (jobdir = jobdir,)
            sleep(3*SLEEP_TIME)
            info = state(jobdir)
            @info (timestamp = Dates.now(), state = info)
            while info == Jobs.Pending || info == Jobs.Running || info == Jobs.Submitted
                sleep(SLEEP_TIME)
                ninfo = state(jobdir)
                ninfo != info && @info (timestamp = Dates.now(), state = info) 
                info = ninfo
            end
            JOB_QUEUE[][jobdir] = (id, Jobs.Running)
            @info (timestamp = Dates.now(), state = Jobs.PostProcessing)
            Service.outputdata(jobdir, map(x -> splitext(splitpath(x.infile)[end])[1], FileIO.read_job_script(joinpath(jobdir, "job.tt"))[2]))
            JOB_QUEUE[][jobdir] = (JOB_QUEUE[][jobdir][1], Jobs.Completed)
        end
    end
end

# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(s::Server, job_dirs_procs)
    lines = filter(!isempty, readlines(PENDING_JOBS_FILE))
    write(PENDING_JOBS_FILE, "")
    njobs = length(filter(x -> x isa Task, collect(values(job_dirs_procs))))
    if length(lines) + njobs > s.max_concurrent_jobs
        
        to_submit = lines[1:s.max_concurrent_jobs - njobs]
        open(PENDING_JOBS_FILE, "a") do f
            for l in lines[s.max_concurrent_jobs - njobs + 1:end]
                write(f, l * "\n")
            end
        end
    else
        to_submit = lines
    end
        
    if !isempty(to_submit)
        while !isempty(to_submit)
            j = pop!(to_submit)
            try
                job_dirs_procs[j] = spawn_task(s, j)
            catch
                sleep(SLEEP_TIME)
                push!(to_submit, j)
            end
        end
    end
end

function spawn_worker(s::Server, jobdir)
    wdir = joinpath(jobdir, ".workflow")
        
    julia_exec = joinpath(Sys.BINDIR, "julia")
    envpath = joinpath(jobdir, ".workflow/run.jl")
    ppath = joinpath(wdir, "Project.toml")
    proc = run(pipeline(Cmd(`$julia_exec --project=$ppath --startup-file=no $envpath $jobdir`, detach=true), stderr = joinpath(wdir, "errors.log"), stdout = joinpath(wdir, "out.log")), wait=false) 
    return proc
end

function handle_workflow_submission!(s::Server, job_dirs_procs)
    lines = filter(!isempty, readlines(PENDING_WORKFLOWS_FILE))
    write(PENDING_WORKFLOWS_FILE, "")
    for l in lines
        job_dirs_procs[l] = spawn_worker(s, l)
    end
end
