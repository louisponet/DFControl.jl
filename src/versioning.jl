const VERSION_DIR_NAME = ".versions"

function job_versions(dir::String)
    verdir = joinpath(dir, VERSION_DIR_NAME)
    if ispath(verdir)
        return parse.(Int, readdir(joinpath(dir, VERSION_DIR_NAME)))
    else
        return Int[]
    end
end
versions(job::DFJob) = job_versions(job.local_dir)
version(job::DFJob) = job.version

function last_job_version(dir::String)
    versions = job_versions(dir)
    return isempty(versions) ? 0 : versions[end]
end
last_version(job::DFJob) = last_job_version(job.local_dir)

version_path(dir::String, version::Int) = joinpath(dir, VERSION_DIR_NAME, "$version")
version_path(job::DFJob) = version_path(job.local_dir, job.version)

exists_version(dir::String, version::Int) = version ∈ job_versions(dir)

function maybe_increment_version(job::DFJob)
    versions_path = joinpath(job, VERSION_DIR_NAME)
    if !ispath(versions_path)
        mkpath(versions_path)
    end
    if ispath(joinpath(job, "job.tt"))
        tjob = DFJob(job.local_dir, version = last_version(job) + 1)
        vpath = version_path(tjob)
        mkpath(vpath)
        cp(job, vpath)
        job.version = last_version(job) + 1
    end
end

function switch_version(job::DFJob, version)
    @assert !isrunning(job; print=false) "Can't switch job versions on a running job."
    cur_version = job.version
    if version != cur_version
        verpath = version_path(job.local_dir, version)
        if !ispath(verpath)
            err_str = "Requested job version ($version) is invalid, please choose from:\n"
            for jv in versions(job)
                err_str *= "\t$jv\n"
            end
            @error err_str
        else
            maybe_increment_version(job)
            out = DFJob(verpath)
            curdir = job.local_dir
            mv(out, job.local_dir, force=true)
            rm(verpath, recursive=true)
            for f in fieldnames(DFJob)
                setfield!(job, f, getfield(out,f))
            end
            set_localdir!(job, curdir)
            job.version = version
        end
    end
    return job
end
