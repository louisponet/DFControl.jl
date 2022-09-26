using StructTypes
Base.@kwdef mutable struct QueueInfo
    full_queue::Dict{String, Tuple{Int, Jobs.JobState}} = Dict{String, Tuple{Int, Jobs.JobState}}()
    current_queue::Dict{String, Tuple{Int, Jobs.JobState}} = Dict{String, Tuple{Int, Jobs.JobState}}()
    submit_queue::Vector{String} = String[]
end
StructTypes.StructType(::Type{QueueInfo}) = StructTypes.Mutable()

const JOB_QUEUE = Ref{QueueInfo}()

function main_loop(s::Server, submit_channel)

    # while !Servers.is_reachable(s.scheduler)
    #     @warn "Can't reach scheduler..."
    #     sleep(SLEEP_TIME)
    # end
    JOB_QUEUE[] = QueueInfo()
    queue!(JOB_QUEUE[], s.scheduler, true)
    @info (timestamp = Dates.now(), username = ENV["USER"], host = gethostname(), pid=getpid())

    # Used to identify if multiple servers are running in order to selfdestruct 
    log_mtimes = mtime.(joinpath.((config_path("logs/runtimes/"),), readdir(config_path("logs/runtimes/"))))

    queuelock = ReentrantLock()
    t = Threads.@spawn while true
        lock(queuelock)
        try
            queue!(JOB_QUEUE[], s.scheduler, false)
        catch e
            @error "queue error" e
        end
        unlock(queuelock)
        sleep(SLEEP_TIME)
    end
    Threads.@spawn while true
        try
            handle_job_submission!(JOB_QUEUE[], s, queuelock, submit_channel)
        catch e
            @error "job submission error" e
        end
        sleep(SLEEP_TIME)
    end
    Threads.@spawn while true
        monitor_issues(log_mtimes)

        try  
            print_log(JOB_QUEUE[])
        catch
            @error "logging error" e
        end
        if ispath(config_path("self_destruct"))
            @info (timestamp = Dates.now(), message = "self_destruct found, self destructing...")
            exit()
        end
            
        sleep(SLEEP_TIME)
    end
    fetch(t)
end

function print_log(queue)
    # @info (timestamp = string(Dates.now()), njobs = length(queue.current_queue), nprocs = nprocs)
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
    daemon_log = config_path("logs/daemon/restapi.log")
    if filesize(daemon_log) > 1e9
        open(daemon_log, "w") do f
            write(f, "")
        end
    end
end
   
# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(queue, s::Server, queuelock, submit_channel)
    lines = queue.submit_queue
    njobs = length(queue.current_queue)
    while !isempty(submit_channel)
        push!(lines, take!(submit_channel))
    end
    n_submit = min(s.max_concurrent_jobs - njobs, length(lines))
    for i in 1:n_submit
        j = lines[i]
        if ispath(j)
            curtries = 0
            while -1 < curtries < 3
                try
                    id = Servers.submit(s.scheduler, j)
                    @info (timestamp = Dates.now(), jobdir = j, jobid = id, state = Jobs.Pending)
                    lock(queuelock)
                    queue.current_queue[j] = (id, Jobs.Pending)
                    unlock(queuelock)
                    curtries = -1
                catch e
                    curtries += 1
                    sleep(SLEEP_TIME)
                    @error e
                end
            end
            if curtries != -1
                push!(lines, j)
            end
        else
            @warn "Submission job at dir: $j is not a directory."
        end
    end
    deleteat!(lines, 1:n_submit)
end

