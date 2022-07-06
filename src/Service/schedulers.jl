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
        for qu in (q.full_queue, q.current_queue)
            for (dir, info) in qu
                if !isdir(dir)
                    delete!(qu, dir)
                    continue
                end
                id = info[1]
                if in_queue(info[2])
                    st = jobstate(s, id)
                    if in_queue(st)
                        delete!(q.full_queue, dir)
                        q.current_queue[dir] = (id, st)
                    else
                        q.full_queue[dir] = (id, st)
                    end
                end
            end
        end
    end
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
    for d in keys(q.full_queue)
        if !isdir(d)
            delete!(q.full_queue, d)
        end
    end
    for (k, v) in squeue
        q.current_queue[k] = v
    end
        
    JSON3.write(QUEUE_FILE(), q)
    return q
end

## BASH ##
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
function queue(sc::Slurm)
    qlines = readlines(`squeue -u $(ENV["USER"]) --format="%Z %i %T"`)[2:end]
    return Dict([(s = split(x); s[1] => (parse(Int, s[2]), jobstate(sc, s[3]))) for x in qlines])
end

function Servers.jobstate(s::Slurm, id::Int)
    cmd = `sacct -u $(ENV["USER"]) --format=State -j $id -P`
    lines = readlines(cmd)
    if length(lines) <= 1
        return Jobs.Unknown
    end
    return jobstate(s, lines[2])
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
function queue(sc::HQ)
    all_lines = readlines(`hq job list`)
    start_id = findnext(x -> x[1] == '+', all_lines, 2)
    endid = findnext(x -> x[1] == '+', all_lines, start_id + 1)
    if endid === nothing
        return Dict()
    end
    qlines = all_lines[start_id+1:endid-1]

    jobinfos = [(s = split(x); (parse(Int, s[2]), jobstate(sc, s[6]))) for x in qlines]

    workdir_line_id =
        findfirst(x-> occursin("Working directory", x), readlines(`hq job info $(jobinfos[1][1])`))

    workdir(id) = split(readlines(`hq job info $id`)[workdir_line_id])[end-1]
    
    return Dict([workdir(x[1]) => x for x in jobinfos])
end

function Servers.jobstate(s::HQ, id::Int)
    lines = readlines(`hq job info $id`)
    
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

function Servers.submit(::HQ, j::String)
    chmod(joinpath(j, "job.sh"), 0o777)

    out = read(Cmd(`hq submit ./job.sh`, dir=j), String)
    if !occursin("successfully", out)
        error("Submission error for job in dir $j.")
    end
    return parse(Int, split(out)[end])
end

Servers.abort(::HQ, id::Int) = 
    run(`hq job cancel $id`)
