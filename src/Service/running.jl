using StructTypes
mutable struct QueueInfo
    full_queue::Dict{String, Tuple{Int, Jobs.JobState}}
    current_queue::Dict{String, Tuple{Int, Jobs.JobState}}
end
StructTypes.StructType(::Type{QueueInfo}) = StructTypes.Mutable()

const JOB_QUEUE = Ref{QueueInfo}()

function main_loop(s::Server, queue)

    JOB_QUEUE[] = QueueInfo(Dict{String, Tuple{Int, Jobs.JobState}}(), Dict{String, Tuple{Int, Jobs.JobState}}())
    queue!(JOB_QUEUE[], s.scheduler, true)
    @info (timestamp = Dates.now(), username = ENV["USER"], host = gethostname(), pid=getpid())

    # Used to identify if multiple servers are running in order to selfdestruct 
    log_mtimes = mtime.(joinpath.((config_path("logs/runtimes/"),), readdir(config_path("logs/runtimes/"))))
    
    while true
        queue!(JOB_QUEUE[], s.scheduler, false)
        handle_job_submission!(JOB_QUEUE[], s)
        # handle_workflow_submission!(s, job_dirs_procs)
        monitor_issues(log_mtimes)
        print_log(JOB_QUEUE[])
        if ispath(config_path("self_destruct"))
            @info (timestamp = Dates.now(), message = "self_destruct found, self destructing...")
            exit()
        end
            
        sleep(SLEEP_TIME)
    end
end

function print_log(queue)
    @info (timestamp = Dates.now(), njobs = length(queue.current_queue), nprocs = nprocs)
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
   
# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(queue, s::Server)
    lines = filter(!isempty, readlines(PENDING_JOBS_FILE))
    write(PENDING_JOBS_FILE, "")
    njobs = length(queue.current_queue)
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
                id = Servers.submit(s.scheduler, j)
                @info (timestamp = Dates.now(), jobdir = j, jobid = id, state = Jobs.Submitted)
                queue.current_queue[j] = (id, Jobs.Submitted)
            catch
                sleep(SLEEP_TIME)
                push!(to_submit, j)
            end
        end
    end
end

### OBSOLETE ###
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
