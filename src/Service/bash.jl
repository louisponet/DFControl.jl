function bash_jobid(job::Job)
    i = last_running_calculation(job.dir)
    l = job[i]
    codeexec = l.exec.exec
    pids = parse.(Int, split(read(`pgrep $codeexec`, String)))
    if isempty(pids)
        return (-1, Jobs.Completed)
    end
    pwds = map(x -> abspath(split(strip(read(`pwdx $x`, String)))[end]), pids)
    id = findfirst(isequal(abspath(job)), pwds)
    if id !== nothing
        return pids[id]
    end
end

function bash_jobstate(jobdir::String)
    job = load(Job(jobdir))
    n = now()
    u = username()
    try
        id = bash_jobid(job)
        if id === nothing
            return (-1, Jobs.Completed)
        else
            return (id, Jobs.Running)
        end
    catch
        return (-1, Jobs.Completed)
    end
end

function bash_queue!(q, init)
    if init
        jobdirs = filter(exists_job, readlines(RUNNING_JOBS_FILE))
    else
        jobdirs = keys(q)
    end
    to_running = String[]
    for j in jobdirs
        info = bash_jobstate(j)
        q[j] = info
        if info[2] == Jobs.Running
            push!(to_running, j)
        end
    end
    write(RUNNING_JOBS_FILE, join(to_running, "\n"))
    return q
end

function bash_submit(j::String)
    run(Cmd(`bash job.tt`, detach=true, dir=j), wait=false)
    open(RUNNING_JOBS_FILE, "a") do f
        write(f, j * "\n")
    end
    return (-1, Jobs.Submitted)
end

function bash_abort(jobdir::String)
    job = load(Job(jobdir))
    id = bash_jobid(job)
    if id !== nothing
        run(`pkill $id`)
    end
    return id
end


