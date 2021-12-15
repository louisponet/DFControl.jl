function bash_jobstate(jobdir::String)
    job = load_job(jobdir)
    n = now()
    u = username()
    i = last_running_calculation(jobdir)
    i === nothing && return false
    l = job[i]
    codeexec = l.exec.exec
    try
        pids = parse.(Int, split(read(`pgrep $codeexec`, String)))
        if isempty(pids)
            return (-1, Jobs.Completed)
        end
        pwd = split(strip(read(`pwdx $(pids[end])`, String)))[end]
        return (pids[end], abspath(pwd) == abspath(job) ? Jobs.Running : Jobs.Completed)
    catch
        return (-1, Jobs.Completed)
    end
end

function bash_queue!(q, init)
    for j in keys(q)
        if exists_job(j)
            q[j] = bash_jobstate(j)
        end
    end
    return q
end

function bash_submit(j::String)
    cd(j)
    run(Cmd(`bash job.tt`, detach=true), wait=false)
    open(RUNNING_JOBS_FILE, "a") do f
        write(f, j * "\n")
    end
    return (-1, Jobs.Submitted)
end
        
