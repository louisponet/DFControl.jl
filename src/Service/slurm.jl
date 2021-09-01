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

"""
    slurm_jobid(job::Job)

Looks through the jobs since the `startdate` and returns the job ID if found.
Returns -1 if the jobID was not found in the list of jobs since `startdate`.
"""
function slurm_jobid(job::Job)
    if haskey(job.metadata, :slurmid)
        return job.metadata[:slurmid]
    end
    id_dir = filter(x -> length(x) == 2,
                    split.(slurm_process_command(`sacct --format=JobID,Workdir%100`)))
    id_ = -1
    for (id, dir) in id_dir
        if occursin(job.dir, dir) || occursin(dir, job.dir)
            id_ = parse(Int, id)
            break
        end
    end
    if id_ == -1
        @info "Job in directory $(job.dir) was not found in the slurm jobs since
        $startdate"
    else
        job.metadata[:slurmid] = id_
    end
    return id_
end

"""
    slurm_isrunning(job::Job)

Returns whether the job is running.
"""
function slurm_isrunning(job::Job)
    id = slurm_jobid(job)
    try
        if id != -1
            result = readlines(`squeue -j $id`)
            st_id = findfirst(x -> x=="ST", split(result[1]))
            return split(result[2])[st_id] ∈ ("R", "PD", "CF", "CG")
        else
            @warn "No jobid found. Was your job submitted through slurm?"
            return false
        end
    catch
        return false
    end
end

"""
    slurm_isqueued(job::Job)

Returns whether the job is queued.
"""
function slurm_isqueued(job::Job)
    id = slurm_jobid(job)
    if id != -1
        if slurm_process_command(`sacct -j $id --format=state`)[1] == "PENDING"
            return true
        else
            return false
        end
    else
        return false
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
