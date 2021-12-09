const MAX_CONCURRENT_JOBS = 1000

const JOB_QUEUE = Ref(Dict{String, Tuple{Int, Jobs.JobState}}())

function main_loop(s::Server)
    job_dirs_procs = Dict{String,Task}()
    JOB_QUEUE[] = queue(s, true)
    for (j, info) in JOB_QUEUE[]
        if info[2] == Jobs.Pending || info[2] == Jobs.Running
            job_dirs_procs[j] = spawn_worker(j)
        end
    end
    while true
        handle_workflow_runners!(job_dirs_procs)
        if length(job_dirs_procs) < MAX_CONCURRENT_JOBS
            handle_job_submission!(s, job_dirs_procs)
        end
        q = queue(s)
        JOB_QUEUE[] = merge(JOB_QUEUE[], q) 
        sleep(SLEEP_TIME)
    end
end

function handle_workflow_runners!(job_dirs_procs)
    to_rm = String[]
    for (k, t) in job_dirs_procs
        if istaskdone(t)
            if istaskfailed(t)
                @error "Task failed with $(t.result)"
            else
                t_ = fetch(t)
                @info """Workflow for job directory $(k) done."""
            end
            push!(to_rm, k)
        end
    end
    if !isempty(to_rm)
        pop!.((job_dirs_procs,), to_rm)
    end
end

function spawn_worker(jobdir)
    # TODO: implement workflows again
    # if ispath(joinpath(jobdir, ".workflow/environment.jld2"))
    #     env_dat = JLD2.load(joinpath(jobdir, ".workflow/environment.jld2"))
    #     proc = addprocs(1; exeflags = "--project=$(env_dat["project"])")[1]
    #     using_expr = Base.Meta.parse("""using $(join(env_dat["modules"], ", "))""")
    #     ctx_path = joinpath(jobdir, ".workflow/ctx.jld2")

    #     Distributed.remotecall_eval(Distributed.Main, proc, using_expr)
        # f = remotecall(DFControl.run_queue, proc, job, JLD2.load(ctx_path)["ctx"];
                       # sleep_time = SLEEP_TIME)
        # return proc, f
    # else
    return Threads.@spawn begin
        sleep(3*SLEEP_TIME)
        while isrunning(jobdir)
            sleep(SLEEP_TIME)
        end
        outputdata(jobdir, map(x -> splitext(splitpath(x.infile)[end])[1], FileIO.read_job_script(joinpath(jobdir, "job.tt"))[2]))
    end
end


# Jobs are submitted by the daemon, using supplied job jld2 from the caller (i.e. another machine)
# Additional files are packaged with the job
function handle_job_submission!(s::Server, job_dirs_procs)
    s = DFC.Server("localhost")
    if ispath(PENDING_JOBS_FILE)
        lines = filter(!isempty, readlines(PENDING_JOBS_FILE))
        write(PENDING_JOBS_FILE, "")
        if length(lines) + length(job_dirs_procs) > MAX_CONCURRENT_JOBS
            
            to_submit = lines[1:MAX_CONCURRENT_JOBS - length(job_dirs_procs)]
            open(PENDING_JOBS_FILE, "a") do f
                for l in lines[MAX_CONCURRENT_JOBS - length(job_dirs_procs) + 1:end]
                    write(f, l * "\n")
                end
            end
        else
            to_submit = lines
        end
        
        if !isempty(to_submit)
            curdir = pwd()
            while !isempty(to_submit)
                j = pop!(to_submit)
                try
                    info = submit(s, j)
                    JOB_QUEUE[][j] = info
                    job_dirs_procs[j] = spawn_worker(j)
                catch
                    sleep(SLEEP_TIME)
                    push!(to_submit, j)
                end
            end
            cd(curdir)
        end
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
            Service.outputdata(job)
            mv(step_file, joinpath(fd, f); force = true)
        end
    end
    return true
end
