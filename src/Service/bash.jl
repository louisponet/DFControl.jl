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

function bash_queue(init)
    if ispath(RUNNING_JOBS_FILE)
        lines = filter(exists_job, readlines(RUNNING_JOBS_FILE))
        d = Dict([l => bash_jobstate(l) for l in lines])
        to_pop = String[]
        for (dir, state) in d
            if state[2] != Jobs.Running
                push!(to_pop, dir)
            end
        end
        for dir in to_pop
            pop!(d, dir)
        end
        open(RUNNING_JOBS_FILE, "w") do f
            for k in keys(d)
                write(f, k * "\n")
            end
        end
        return d
    else
        if init
            touch(RUNNING_JOBS_FILE)
        end
        return Dict{String, Tuple{Int, Jobs.JobState}}()
    end
end

function bash_submit(j::String)
    cd(j)
    run(`bash job.tt`)
    open(RUNNING_JOBS_FILE, "a") do f
        write(f, j * "\n")
    end
    return (-1, Jobs.Submitted)
end
        
