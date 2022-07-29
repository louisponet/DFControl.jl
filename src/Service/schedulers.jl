using ..Servers: Bash, Slurm, HQ, Scheduler, jobstate

function queue!(q, s::Scheduler, init)

    if init
        if ispath(QUEUE_FILE())
            t = read(QUEUE_FILE())
            if !isempty(t)
                tq = JSON3.read(t, QueueInfo)
                copy!(q.full_queue, tq.full_queue)
                copy!(q.current_queue, tq.current_queue)
            end
        end
    end
    for qu in (q.full_queue, q.current_queue)
        for (dir, info) in qu
            if !isdir(dir)
                delete!(qu, dir)
            end
        end
    end
    # Here we check whether the scheduler died while the server was running and try to restart and resubmit   
    if maybe_scheduler_restart(s)
        for (d, i) in q.current_queue
            if isdir(d)
                q.full_queue[d] = (-1, Jobs.Saved)
                submit(d, false)
            end
            pop!(q.current_queue, d)
        end
    else
        squeue = queue(s)
        
        for (d, i) in q.current_queue
            if haskey(squeue, d)
                state = pop!(squeue, d)[2]
            else
                state = jobstate(s, i[1])
            end
            if in_queue(state)
                q.current_queue[d] = (i[1], state)
            else
                delete!(q.current_queue, d)
                q.full_queue[d] = (i[1], state)
            end
        end
        for (k, v) in squeue
            q.current_queue[k] = v
        end
            
    end
    
    JSON3.write(QUEUE_FILE(), q)
    return q
end

## BASH ##
maybe_scheduler_restart(::Bash) = false
queue(::Bash) = Dict()

function Servers.jobstate(::Bash, id::Int)
    out = Pipe()
    err = Pipe()
    p = run(pipeline(ignorestatus(`ps -p $id`), stderr = err, stdout=out))
    close(out.in)
    close(err.in)
    return p.exitcode == 1 ? Jobs.Completed : Jobs.Running
end

Servers.submit(::Bash, j::String) = 
    Int(getpid(run(Cmd(`bash job.sh`, detach=true, dir=j), wait=false)))

function Servers.abort(::Bash, id::Int)
    pids = [parse(Int, split(s)[1]) for s in readlines(`ps -s $id`)[2:end]]
    for p in pids
        run(ignorestatus(`kill $p`))
    end
end

in_queue(s::Jobs.JobState) =
    s in (Jobs.Submitted, Jobs.Pending, Jobs.Running, Jobs.Configuring, Jobs.Completing, Jobs.Suspended)
    
## SLURM ##
function maybe_scheduler_restart(::Slurm)
    if occursin("error", read(run(`squeue --me`), String))
        if occursin("(null)", read(run(Cmd(`slurmd`)), String))
            error("Can not start slurmctld automatically...")
        else
            return true
        end
    else
        return false
    end
end

function queue(sc::Slurm)
    qlines = readlines(`squeue -u $(ENV["USER"]) --format="%Z %i %T"`)[2:end]
    return Dict([(s = split(x); s[1] => (parse(Int, s[2]), jobstate(sc, s[3]))) for x in qlines])
end

function Servers.jobstate(s::Slurm, id::Int)
    cmd = `sacct -u $(ENV["USER"]) --format=State -j $id -P`
    try
        lines = readlines(cmd)
        if length(lines) <= 1
            return Jobs.Unknown
        end
        return jobstate(s, lines[2])
    catch
        return Jobs.Unknown
    end
end

function Servers.jobstate(::Slurm, state::AbstractString)
    if state == "PENDING"
        return Jobs.Pending
    elseif state == "RUNNING"
        return Jobs.Running
    elseif state == "COMPLETED"
        return Jobs.Completed
    elseif state == "CONFIGURING"
        return Jobs.Configuring
    elseif state == "COMPLETING"
        return Jobs.Completing
    elseif state == "CANCELLED"
        return Jobs.Cancelled
    elseif state == "BOOT_FAIL"
        return Jobs.BootFail
    elseif state == "DEADLINE"
        return Jobs.Deadline
    elseif state == "FAILED"
        return Jobs.Failed
    elseif state == "NODE_FAIL"
        return Jobs.NodeFail
    elseif state == "OUT_OF_MEMORY"
        return Jobs.OutOfMemory
    elseif state == "PREEMTED"
        return Jobs.Preempted
    elseif state == "REQUEUED"
        return Jobs.Requeued
    elseif state == "RESIZING"
        return Jobs.Resizing
    elseif state == "REVOKED"
        return Jobs.Revoked
    elseif state == "SUSPENDED"
        return Jobs.Suspended
    elseif state == "TIMEOUT"
        return Jobs.Timeout
    end
    return Jobs.Unknown
end

Servers.submit(::Slurm, j::String) =
    parse(Int, split(read(Cmd(`sbatch job.sh`, dir=j), String))[end])

Servers.abort(s::Slurm, id::Int) = 
    run(`scancel $id`)

#### HQ
function maybe_scheduler_restart(sc::HQ)
    function readinfo()
        out = Pipe()
        err = Pipe()
        run(pipeline(Cmd(Cmd(string.([split(sc.server_command)..., "server", "info"])), ignorestatus=true), stdout=out, stderr=err))
        close(out.in)
        close(err.in)
        return read(err, String)
    end
    
    if occursin("No online", readinfo())
        run(Cmd(Cmd(string.([split(sc.server_command)..., "server", "start"])), detach=true), wait=false)
        sleep(0.01)
        tries = 0
        while occursin("No online", readinfo()) && tries < 10
            sleep(0.01)
            tries += 1
        end
        if tries == 10
            error("HQ server not reachable")
        end

        maybe_restart_allocs(sc)
        return true
    else
        maybe_restart_allocs(sc)
        return false
    end
end

function maybe_restart_allocs(sc::HQ)
    alloc_lines = readlines(Cmd(Cmd(string.([split(sc.server_command)..., "alloc", "list"]))))
    if length(alloc_lines) == 3 # no allocs -> add all
        allocs_to_add = sc.allocs
    else
        alloc_args = map(a -> replace(strip(split(a, "|")[end]), "," => " "), alloc_lines[4:end-1])
        allocs_to_add = filter(a -> !any(x -> x == strip(split(a, "--")[end]), alloc_args), sc.allocs)
    end
        
    for ac in allocs_to_add
        run(Cmd(string.([split(sc.server_command)..., "alloc", "add", split(ac)...])))
    end
end
        
function queue(sc::HQ)
    all_lines = readlines(Cmd(string.([split(sc.server_command)..., "job", "list"])))

    start_id = findnext(x -> x[1] == '+', all_lines, 2)
    endid = findnext(x -> x[1] == '+', all_lines, start_id + 1)
    if endid === nothing
        return Dict()
    end
    qlines = all_lines[start_id+1:endid-1]

    jobinfos = [(s = split(x); (parse(Int, s[2]), jobstate(sc, s[6]))) for x in qlines]

    workdir_line_id =
        findfirst(x-> occursin("Working directory", x), readlines(Cmd(string.([split(sc.server_command)..., "job", "info", "$(jobinfos[1][1])"]))))

    workdir(id) = split(readlines(Cmd(string.([split(sc.server_command)..., "job", "info", "$id"])))[workdir_line_id])[end-1]
    
    return Dict([workdir(x[1]) => x for x in jobinfos])
end

function Servers.jobstate(s::HQ, id::Int)
    lines = readlines(Cmd(string.([split(s.server_command)..., "job", "info", "$id"])))
    
    if length(lines) <= 1
        return Jobs.Unknown
    end
    return jobstate(s, split(lines[findfirst(x->occursin("State", x), lines)])[4])
end

function Servers.jobstate(::HQ, state::AbstractString)
    if state == "WAITING"
        return Jobs.Pending
    elseif state == "RUNNING"
        return Jobs.Running
    elseif state == "FINISHED"
        return Jobs.Completed
    elseif state == "CANCELED"
        return Jobs.Cancelled
    elseif state == "FAILED"
        return Jobs.Failed
    end
    return Jobs.Unknown
end

function Servers.submit(h::HQ, j::String)
    chmod(joinpath(j, "job.sh"), 0o777)

    out = read(Cmd(Cmd(string.([split(h.server_command)..., "submit", "./job.sh"])), dir=j), String)
    if !occursin("successfully", out)
        error("Submission error for job in dir $j.")
    end
    return parse(Int, split(out)[end])
end

Servers.abort(h::HQ, id::Int) = 
    run(Cmd(string.([split(h.server_command)..., "job", "cancel", "$id"])))
