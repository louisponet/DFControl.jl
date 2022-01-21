slurm_process_command(cmd) = strip.(reverse(readlines(cmd)))[1:end-2]
"""
    slurm_history_jobdir(stardate=yesterday())

Returns the unique job directories of the jobs that ran since the `startdate`.
Startdate should be printed in following format: yyyy-mm-dd.
"""
function slurm_history_jobdir(startdate = yesterday()) #format of startdate = yyyy-mm-dd
    history = slurm_process_command(`sacct --starttime $startdate --format=Workdir%100`)
    output = String[]
    for h in history
        if h ∉ output && ispath(h)
            push!(output, h)
        end
    end
    return output
end

function slurm_state(state)
    if state == "PENDING"
        return Jobs.Pending
    elseif state == "RUNNING"
        return Jobs.Running
    elseif state == "COMPLETED"
        return Jobs.Completed
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
end

function slurm_queue!(q, all=false)
    if all # Should only be done veeeery sporadically
        cmd = `sacct -u $(ENV["USER"]) --format=Workdir%100,JobID%20,State%30 -S 2000-01-01`
    else
        cmd = `sacct -u $(ENV["USER"]) --format=Workdir%100,JobID%20,State%30`
    end
    all = map(x -> strip.(filter(!isempty, split(x, "  "))), readlines(cmd)[3:end])
    lines = filter(x->length(x) == 3, all)
    return merge!(q, Dict([l[1] => (parse(Int, l[2]), slurm_state(occursin("by", l[3]) ? split(l[3])[1] : l[3])) for l in lines]))
end

"""
    slurm_jobid(job::Job)

Looks through the jobs since the `startdate` and returns the job ID if found.
Returns -1 if the jobID was not found in the list of jobs since `startdate`.
"""
function slurm_jobid(job::Job)
    if haskey(JOB_QUEUE[], job.dir)
        return JOB_QUEUE[][job.dir][1]
    else
        @info "Job in directory $(job.dir) was not found in the slurm jobs."
        return -1
    end
end

"""
    slurm_jobstate(job::Job)

Returns the current state of a job.
"""
function slurm_job_state(job::Job)
    if haskey(JOB_QUEUE[], job.dir)
        return JOB_QUEUE[][job.dir][2]
    else
        error("Job in directory $(job.dir) was not found in the slurm jobs.")
        return -1
    end
end

"""
    slurm_mostrecent(index=1, jobfile="job.tt", startdate=lastmonth(), args...; kwargs...)

Returns whether the `index`th most recent job with job script `jobfile`.
Extra args and kwargs will be passed to the `Job` constructor.
"""
function slurm_mostrecent(index = 1, jobfile = "job.tt", startdate = lastmonth(), args...;
                          kwargs...)
    dirs = slurm_history_jobdir(startdate)
    @assert length(dirs) >= index "There are less recent jobs since startdate $startdate than the index $index."
    jobpath = joinpath(dirs[index], jobfile)
    @assert ispath(jobpath) "The directory $(dirs[index]) does not have a job file with filename $jobfile."
    return Job(dirs[index], args...; job_fuzzy = jobfile, kwargs...)
end

function slurm_submit(j::String)
    id = parse(Int, split(read(Cmd(`sbatch job.tt`, dir=j), String))[end])
    return (id, Jobs.Submitted)
end

slurm_abort(id::Int) = (run(`scancel $id`); return id)
