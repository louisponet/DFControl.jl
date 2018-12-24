"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
qstat(server) = server=="localhost" ? run(`qstat`) : sshcmd(server, "qstat")
qstat()       = qstat(getdefault_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
watch_qstat(server) = server=="localhost" ? run(`watch qstat`) : sshcmd(server, "watch qstat")
watch_qstat()       = watch_qstat(getdefault_server())

#--------------- Slurm Interactions -------------------#
"""
    slurm_history_jobdir(stardate=yesterday())

Returns the unique job directories of the jobs that ran since the `startdate`.
Startdate should be printed in following format: yyyy-mm-dd.
"""
function slurm_history_jobdir(startdate=yesterday()) #format of startdate = yyyy-mm-dd
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
    slurm_jobid(job::DFJob, startdate=yesterday())

Looks through the jobs since the `startdate` and returns the job ID if found.
Returns -1 if the jobID was not found in the list of jobs since `startdate`.
"""
function slurm_jobid(job::DFJob, startdate=yesterday())
    runslocal_assert(job)
    if haskey(job.metadata, :slurmid)
        return job.metadata[:slurmid]
    end
    id_dir = filter(x->length(x) ==2, split.(slurm_process_command(`sacct --starttime $startdate --format=JobID,Workdir%100`)))
    id_ = -1
    for (id, dir) in id_dir
        if dir == job.local_dir[1:end-1]
            id_ = parse(Int, id)
            break
        end
    end
    if id_ == -1
        @info "Job in directory $(job.local_dir) was not found in the slurm jobs since
        $startdate"
    else
        job.metadata[:slurmid] = id_
    end
    return id_
end

"""
    slurm_isrunning(job::DFJob)

Returns whether the job is running.
"""
function slurm_isrunning(job::DFJob)
    runslocal_assert(job)
    id = slurm_jobid(job)
    if id != -1
        if slurm_process_command(`sacct -j $id --format=state`)[1] == "RUNNING"
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
Extra args and kwargs will be passed to the `DFJob` constructor.
"""
function slurm_mostrecent(index=1, jobfile="job.tt", startdate=lastmonth(), args...; kwargs...)
    dirs = slurm_history_jobdir(startdate)
    @assert length(dirs) >= index "There are less recent jobs since startdate $startdate than the index $index."
    jobpath = joinpath(dirs[index], jobfile)
    @assert ispath(jobpath) "The directory $(dirs[index]) does not have a job file with filename $jobfile."
    return DFJob(dirs[index], args...; job_fuzzy=jobfile, kwargs...)
end
