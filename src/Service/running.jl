const JOB_QUEUE = Ref(Dict{String, Tuple{Int, Jobs.JobState}}())

function main_loop(s::Server)
    job_dirs_procs = Dict()
    queue!(JOB_QUEUE[], s, true)
    for (j, info) in JOB_QUEUE[]
        if info[2] == Jobs.Running 
            job_dirs_procs[j] = spawn_worker(j)
        end
    end
    while true
        queue!(JOB_QUEUE[], s)
        handle_workflow_runners!(job_dirs_procs)
        handle_job_submission!(s, job_dirs_procs)
        handle_workflow_submission!(s, job_dirs_procs)
        sleep(SLEEP_TIME)
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
                    @info """Workflow for job directory $(k) done."""
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
    # TODO: implement workflows again

    id = next_jobid()
    info = submit(s, jobdir)
    JOB_QUEUE[][jobdir] = info
    return Threads.@spawn begin
        with_logger(job_logger(id)) do
            @info (timestamp = Dates.now(), jobdir = jobdir) 
            sleep(3*SLEEP_TIME)
            info = state(jobdir)
            while info == Jobs.Pending || info == Jobs.Running || info == Jobs.Submitted
                sleep(SLEEP_TIME)
                info = state(jobdir)
                @info (timestamp = Dates.now(), jobdir = jobdir, state = info) 
            end
            # outputdata(jobdir, map(x -> splitext(splitpath(x.infile)[end])[1], FileIO.read_job_script(joinpath(jobdir, "job.tt"))[2]))
            JOB_QUEUE[][jobdir] = (info[1], Jobs.Completed)
        end
    end
end

function spawn_worker(s::Server, jobdir)
    wdir = joinpath(jobdir, ".workflow")
        
    julia_exec = joinpath(Sys.BINDIR, "julia")
    envpath = joinpath(jobdir, ".workflow/run.jl")
    ppath = joinpath(wdir, "Project.toml")
    proc = run(pipeline(Cmd(`$julia_exec --project=$ppath --startup-file=no $envpath $jobdir`, detach=true), stderr = joinpath(wdir, "errors.log"), stdout = joinpath(wdir, "out.log")), wait=false) 
    # proc = addprocs(1; exeflags = "--project=$ppath")[1]
    

    # f = Distributed.remotecall(include, proc, envpath)
    # f = remotecall(Main.run_queue, proc, jobdir)
    return proc
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

function handle_workflow_submission!(s::Server, job_dirs_procs)
    lines = filter(!isempty, readlines(PENDING_WORKFLOWS_FILE))
    write(PENDING_WORKFLOWS_FILE, "")
    for l in lines
        job_dirs_procs[l] = spawn_worker(s, l)
    end
end
# function run_queue(jobdir::String)
#     job = load_job(jobdir)
#     logger = workflow_logger(job)
#     with_logger(logger) do
#         qd = queued_dir(job)
#         queued_steps = readdir(qd)
#         order = sortperm(parse.(Int, getindex.(splitext.(queued_steps), 1)))
#         fd = finished_dir(job)
#         if !ispath(qd) || isempty(queued_steps)
#             return
#         end
#         if !ispath(fd)
#             mkpath(fd)
#         else
#             for f in readdir(fd)
#                 rm(joinpath(fd, f))
#             end
#         end
#         for f in queued_steps[order]
#             step_file = joinpath(qd, f)
#             t = include(step_file)
#             results = @eval $(t)($(job))
#             JLD2.jldopen(joinpath(job, ".workflow/results.jld2"), "w+") do j
#                 j[splitext(f)[1]] = results
#             end
#             mv(step_file, joinpath(fd, f); force = true)
#             info = state(jobdir)
#             while info == Jobs.Pending || info == Jobs.Running || info == Jobs.Submitted
#                 sleep(SLEEP_TIME)
#             end
#         end
#     end
#     return true
# end
